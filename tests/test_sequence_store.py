"""Focused regression tests for PyamilySeq's SQLite sequence store and provenance."""

import sqlite3
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from families import (
    FamilyInputError,
    SequenceStore,
    build_family_ledger,
)


class SequenceStoreTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="pyamilyseq_sqlite_test_")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_persistent_database_round_trip(self):
        fasta = self.root / "input.fasta"
        fasta.write_text(
            ">genomeA|gene1\nATGAAATAG\n"
            ">genomeB|gene2\nATGCCCTAA\n",
            encoding="utf-8",
        )
        database = self.root / "checkpoint.sqlite3"

        store = SequenceStore.from_fasta(
            str(fasta), database_path=str(database), memory_limit_mb=64,
            buffer_cap_mb=8,
        )
        self.assertTrue(database.is_file())
        self.assertEqual(len(store), 2)
        self.assertEqual(store.sequence_length("genomeA|gene1"), 9)
        self.assertEqual(store["genomeB|gene2"], "ATGCCCTAA")
        store.close()

        self.assertTrue(database.is_file(), "persistent database must survive close()")
        reopened = SequenceStore.from_existing_db(str(database))
        self.assertEqual(len(reopened), 2)
        self.assertEqual(reopened.sequence_length("genomeA|gene1"), 9)
        reopened.close()
        self.assertTrue(database.is_file())

    def test_temporary_database_is_removed_on_close(self):
        fasta = self.root / "input.fasta"
        fasta.write_text(">g|x\nATGC\n", encoding="utf-8")
        store = SequenceStore.from_fasta(str(fasta), temp_dir=str(self.root))
        path = Path(store.path)
        self.assertTrue(path.exists())
        store.close()
        self.assertFalse(path.exists())
        self.assertFalse(Path(str(path) + "-wal").exists())
        self.assertFalse(Path(str(path) + "-shm").exists())

    def test_existing_non_pyamilyseq_sqlite_file_is_not_modified(self):
        database = self.root / "other.sqlite3"
        conn = sqlite3.connect(database)
        conn.execute("CREATE TABLE unrelated (value TEXT)")
        conn.commit()
        conn.close()

        with self.assertRaises(FamilyInputError):
            SequenceStore.from_existing_db(str(database))

        conn = sqlite3.connect(database)
        tables = {
            row[0] for row in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            )
        }
        conn.close()
        self.assertEqual(tables, {"unrelated"})

    def test_incomplete_checkpoint_is_rejected(self):
        database = self.root / "incomplete.sqlite3"
        store = SequenceStore(path=str(database))
        store._write_metadata({"schema_version": store.DB_SCHEMA_VERSION, "status": "building"})
        store.add("g|x", "ATGC")
        store.close()
        with self.assertRaisesRegex(FamilyInputError, "not a completed checkpoint"):
            SequenceStore.from_existing_db(str(database))

    def test_byte_safe_retrieval_reports_non_ascii_payload(self):
        database = self.root / "bad_payload.sqlite3"
        conn = sqlite3.connect(database)
        conn.execute(
            "CREATE TABLE sequences (sequence_id TEXT PRIMARY KEY, sequence TEXT NOT NULL, length INTEGER NOT NULL)"
        )
        conn.execute(
            "INSERT INTO sequences VALUES (?, CAST(? AS TEXT), ?)",
            ("g|bad", sqlite3.Binary(b"ATG\xff"), 4),
        )
        conn.commit()
        conn.close()

        store = SequenceStore.from_existing_db(str(database))
        self.assertEqual(store.sequence_length("g|bad"), 4)
        with self.assertRaisesRegex(FamilyInputError, "non-ASCII/corrupt bytes"):
            _ = store["g|bad"]
        store.close()


class SupplementationProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="pyamilyseq_provenance_test_")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_original_backfilled_singleton_is_not_supplement_only(self):
        fasta = self.root / "combined.fasta"
        fasta.write_text(
            ">original|clustered\nAAAAAA\n"
            ">original|backfilled\nCCCCCC\n"
            ">supplement|extended\nGGGGGG\n"
            ">supplement|only\nTTTTTT\n",
            encoding="utf-8",
        )
        original_clusters = self.root / "original.clstr"
        original_clusters.write_text(
            ">Cluster 0\n"
            "0 6aa, >original|clustered... *\n",
            encoding="utf-8",
        )
        reclustered = self.root / "reclustered.clstr"
        reclustered.write_text(
            ">Cluster 0\n"
            "0 6aa, >original|clustered... *\n"
            "1 6aa, >supplement|extended... at 99.00%\n"
            ">Cluster 1\n"
            "0 6aa, >supplement|only... *\n",
            encoding="utf-8",
        )

        ledger = build_family_ledger(
            fasta_path=str(fasta),
            cluster_path=str(original_clusters),
            cluster_format="CD-HIT",
            recluster_path=str(reclustered),
            supplemental_sequence_ids={"supplement|extended", "supplement|only"},
        )
        try:
            stages_by_member = {
                member: family.stage
                for family in ledger.families.values()
                for member in family.members
            }
            self.assertEqual(stages_by_member["original|clustered"], "extended")
            self.assertEqual(stages_by_member["original|backfilled"], "First")
            self.assertEqual(stages_by_member["supplement|extended"], "extended")
            self.assertEqual(stages_by_member["supplement|only"], "only_Second")
            self.assertEqual(len(ledger.families), 3)
        finally:
            ledger.sequences.close()

    def test_failed_ledger_build_removes_temporary_database(self):
        fasta = self.root / "input.fasta"
        fasta.write_text(">g|x\nATGC\n", encoding="utf-8")
        missing_clusters = self.root / "missing.clstr"
        options = SimpleNamespace(
            tmp_dir=str(self.root),
            mem=64,
            sequence_buffer_cap_mb=8,
        )
        with self.assertRaises(FileNotFoundError):
            build_family_ledger(
                fasta_path=str(fasta),
                cluster_path=str(missing_clusters),
                cluster_format="CD-HIT",
                options=options,
            )
        leftovers = list(self.root.glob("pyamilyseq_sequences_*.sqlite3*"))
        self.assertEqual(leftovers, [], "temporary database/sidecars should be cleaned on ledger errors")

    def test_supplement_id_cannot_appear_in_original_clusters(self):
        fasta = self.root / "combined.fasta"
        fasta.write_text(">g|new\nAAAAAA\n", encoding="utf-8")
        original_clusters = self.root / "original.clstr"
        original_clusters.write_text(
            ">Cluster 0\n0 6aa, >g|new... *\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(FamilyInputError, "already present in the original cluster file"):
            build_family_ledger(
                fasta_path=str(fasta),
                cluster_path=str(original_clusters),
                cluster_format="CD-HIT",
                recluster_path=str(original_clusters),
                supplemental_sequence_ids={"g|new"},
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
