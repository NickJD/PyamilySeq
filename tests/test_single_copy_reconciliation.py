"""Regression tests for single-copy filtering after First/Second reconciliation.

Run from the PyamilySeq source root with:
    python -m unittest discover -s tests -p 'test_single_copy_reconciliation.py'
"""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

try:
    from families import FamilyInputError, build_family_ledger
    from workflow import run_pangenome_workflow
except ImportError:  # pragma: no cover - package-style test discovery
    from PyamilySeq.families import FamilyInputError, build_family_ledger
    from PyamilySeq.workflow import run_pangenome_workflow


def write_cdhit_clusters(path: Path, clusters: list[list[str]]) -> None:
    """Write a minimal valid CD-HIT .clstr file.

    The first member of each list is marked as the representative.
    """
    rows = []
    for cluster_index, members in enumerate(clusters):
        rows.append(f">Cluster {cluster_index}")
        for member_index, member in enumerate(members):
            marker = "*" if member_index == 0 else "at 99%"
            rows.append(f"{member_index} 100aa, >{member}... {marker}")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


class SingleCopyReconciliationTests(unittest.TestCase):
    def build_ledger(
        self,
        temp_dir: Path,
        first_clusters: list[list[str]],
        second_clusters: list[list[str]],
        all_ids: list[str],
        *,
        mode: str = "Species",
    ):
        fasta_path = temp_dir / "combined.faa"
        first_path = temp_dir / "first.clstr"
        second_path = temp_dir / "second.clstr"

        with fasta_path.open("w", encoding="utf-8") as handle:
            for sequence_id in all_ids:
                handle.write(f">{sequence_id}\nM" + ("A" * 98) + "K\n")

        write_cdhit_clusters(first_path, first_clusters)
        write_cdhit_clusters(second_path, second_clusters)

        first_ids = {sequence_id for cluster in first_clusters for sequence_id in cluster}
        supplemental_ids = set(all_ids) - first_ids
        self.assertTrue(
            supplemental_ids,
            "Each regression scenario should contain at least one supplementary sequence",
        )

        return build_family_ledger(
            fasta_path=str(fasta_path),
            cluster_path=str(first_path),
            cluster_format="CD-HIT",
            mode=mode,
            recluster_path=str(second_path),
            sequence_tag="StORF",
            sequence_delimiter="|",
            supplemental_sequence_ids=supplemental_ids,
        )

    def test_first_and_second_copy_in_same_genome_is_excluded(self):
        first = "GCA_1|geneA"
        other = "GCA_2|geneB"
        supplement = "GCA_1|StORF_001"

        with tempfile.TemporaryDirectory() as directory:
            ledger = self.build_ledger(
                Path(directory),
                [[first, other]],
                [[first, supplement]],
                [first, other, supplement],
            )
            try:
                family = next(iter(ledger.families.values()))
                self.assertEqual(family.stage, "extended")
                self.assertEqual(ledger.family_copy_number(family.family_id), {
                    "GCA_1": 2,
                    "GCA_2": 1,
                })
                self.assertNotIn(
                    family.family_id,
                    ledger.selected_single_copy(0),
                    "A First and Second member from the same genome must count as two copies",
                )
            finally:
                ledger.sequences.close()

    def test_first_and_second_copies_in_different_genomes_are_retained(self):
        first = "GCA_1|geneA"
        other = "GCA_3|geneB"
        supplement = "GCA_2|StORF_001"

        with tempfile.TemporaryDirectory() as directory:
            ledger = self.build_ledger(
                Path(directory),
                [[first, other]],
                [[first, supplement]],
                [first, other, supplement],
            )
            try:
                family = next(iter(ledger.families.values()))
                self.assertEqual(
                    set(ledger.family_copy_number(family.family_id).values()),
                    {1},
                )
                self.assertIn(family.family_id, ledger.selected_single_copy(0))
            finally:
                ledger.sequences.close()

    def test_second_copy_is_counted_after_two_original_families_merge(self):
        first_a = "GCA_1|geneA"
        first_b = "GCA_2|geneB"
        second_a = "GCA_3|geneC"
        second_b = "GCA_4|geneD"
        supplement = "GCA_1|StORF_001"

        with tempfile.TemporaryDirectory() as directory:
            ledger = self.build_ledger(
                Path(directory),
                [[first_a, first_b], [second_a, second_b]],
                [[first_a, second_a, supplement]],
                [first_a, first_b, second_a, second_b, supplement],
            )
            try:
                self.assertEqual(len(ledger.families), 1)
                family = next(iter(ledger.families.values()))
                self.assertEqual(family.stage, "combined")
                self.assertEqual(
                    sum(not label.startswith("singleton:") for label in family.source_clusters),
                    2,
                    "The final family must contain both original source clusters; "
                    "the supplementary singleton is a third source record",
                )
                self.assertEqual(
                    ledger.family_copy_number(family.family_id),
                    {"GCA_1": 2, "GCA_2": 1, "GCA_3": 1, "GCA_4": 1},
                )
                self.assertEqual(ledger.selected_single_copy(0), [])
            finally:
                ledger.sequences.close()

    def test_genus_mode_counts_distinct_genomes_not_only_the_shared_genus(self):
        first = "Genus_001|geneA"
        supplement = "Genus_002|StORF_001"

        with tempfile.TemporaryDirectory() as directory:
            ledger = self.build_ledger(
                Path(directory),
                [[first]],
                [[first, supplement]],
                [first, supplement],
                mode="Genus",
            )
            try:
                family = next(iter(ledger.families.values()))
                self.assertEqual(family.prevalence, 1)
                self.assertEqual(
                    ledger.family_copy_number(family.family_id),
                    {"Genus_001": 1, "Genus_002": 1},
                )
                self.assertIn(family.family_id, ledger.selected_single_copy(0))
            finally:
                ledger.sequences.close()

    def test_tolerance_is_applied_to_final_occupied_genomes(self):
        first_members = [
            "GCA_1|geneA",
            "GCA_2|geneB",
            "GCA_3|geneC",
            "GCA_4|geneD",
        ]
        supplement = "GCA_1|StORF_001"

        with tempfile.TemporaryDirectory() as directory:
            ledger = self.build_ledger(
                Path(directory),
                [first_members],
                [[first_members[0], supplement]],
                first_members + [supplement],
            )
            try:
                family_id = next(iter(ledger.families))
                self.assertEqual(
                    ledger.selected_single_copy(0),
                    [],
                    "Zero tolerance must exclude one multi-copy genome",
                )
                self.assertIn(
                    family_id,
                    ledger.selected_single_copy(25),
                    "One multi-copy genome is allowed at 25% of four occupied genomes",
                )
                self.assertEqual(
                    ledger.selected_single_copy(24.99),
                    [],
                    "The permitted number of multi-copy genomes is rounded down",
                )
            finally:
                ledger.sequences.close()

    def test_workflow_filters_analysis_outputs_but_keeps_full_assignment_ledger(self):
        first = "GCA_1|geneA"
        other = "GCA_2|geneB"
        retained = "GCA_3|geneC"
        supplement = "GCA_1|StORF_001"

        with tempfile.TemporaryDirectory() as directory:
            temp_dir = Path(directory)
            fasta_path = temp_dir / "combined.faa"
            first_path = temp_dir / "first.clstr"
            second_path = temp_dir / "second.clstr"
            output_dir = temp_dir / "output"

            all_ids = [first, other, retained, supplement]
            with fasta_path.open("w", encoding="utf-8") as handle:
                for sequence_id in all_ids:
                    handle.write(f">{sequence_id}\nM" + ("A" * 98) + "K\n")
            write_cdhit_clusters(first_path, [[first, other], [retained]])
            write_cdhit_clusters(second_path, [[first, supplement]])

            options = SimpleNamespace(
                fasta=str(fasta_path),
                clusters=str(first_path),
                cluster_format="CD-HIT",
                output_dir=str(output_dir),
                reclustered=str(second_path),
                sequence_tag="StORF",
                sequence_delimiter="|",
                supplemental_sequence_ids=[supplement],
                taxon_map=None,
                sequence_db=None,
                sequence_db_out=None,
                mem=32,
                sequence_buffer_cap_mb=32,
                tmp_dir=None,
                species_groups="99,95,15",
                core_groups="1-10",
                single_copy_only=True,
                single_copy_tolerance=0.0,
                gene_presence_absence_out=True,
                gene_annotation_files=None,
                baseline_output_dir=None,
                gene_group_table=True,
                core_group_table=True,
                write_groups=None,
                write_individual_groups=False,
                align_core=False,
                visualise=False,
                added_genomes=None,
                threads=1,
            )

            ledger = run_pangenome_workflow(options, "Species")
            self.assertEqual(len(ledger.families), 2)

            with (output_dir / "family_summary.tsv").open(encoding="utf-8") as handle:
                family_rows = list(csv.DictReader(handle, delimiter="\t"))
            with (output_dir / "gene_presence_absence.csv").open(
                encoding="utf-8", newline=""
            ) as handle:
                gpa_rows = list(csv.DictReader(handle))
            with (output_dir / "gene_family_assignments.tsv").open(
                encoding="utf-8"
            ) as handle:
                assignment_rows = list(csv.DictReader(handle, delimiter="\t"))
            manifest = json.loads((output_dir / "run_manifest.json").read_text())

            self.assertEqual(len(family_rows), 1)
            self.assertEqual(len(gpa_rows), 1)
            self.assertEqual(len(assignment_rows), len(all_ids))
            self.assertEqual(manifest["selected_families"], 1)
            self.assertEqual(manifest["single_copy_filter"]["scope"], "final_reconciled_families")
            self.assertEqual(manifest["single_copy_filter"]["excluded_families"], 1)
            statistics = (output_dir / "summary_statistics.txt").read_text()
            self.assertIn("Scope: final reconciled families", statistics)
            self.assertIn("Families retained: 1", statistics)
            self.assertIn("Families excluded: 1", statistics)

    def test_non_finite_or_out_of_range_tolerance_is_rejected(self):
        first = "GCA_1|geneA"
        supplement = "GCA_1|StORF_001"

        with tempfile.TemporaryDirectory() as directory:
            ledger = self.build_ledger(
                Path(directory),
                [[first]],
                [[first, supplement]],
                [first, supplement],
            )
            try:
                for tolerance in (-0.1, 100.1, float("nan"), float("inf")):
                    with self.subTest(tolerance=tolerance):
                        with self.assertRaises(FamilyInputError):
                            ledger.selected_single_copy(tolerance)
            finally:
                ledger.sequences.close()


if __name__ == "__main__":
    unittest.main()
