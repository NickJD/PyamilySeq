import csv
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from utils import read_combined_files, read_separate_files, iter_gene_annotations
from incremental import prepare_incremental_inputs
from workflow import (_family_metadata_consensus, _write_gpa, run_pangenome_workflow,
                      GPA_METADATA_COLUMNS)


class AnnotationExtractionTests(unittest.TestCase):
    def _write_combined_gff(self, directory):
        gff = Path(directory) / "GCA_000001.1.gff3"
        gff.write_text(
            "##gff-version 3\n"
            "NZ_1\tBakta\tCDS\t1\t9\t.\t+\t0\t"
            "ID=geneA;Name=Alpha enzyme;gene=abcA;product=Alpha enzyme%2C subunit A\n"
            "NZ_1\tBakta\tCDS\t10\t18\t.\t+\t0\t"
            "ID=geneB;Name=Hypothetical protein;product=Hypothetical protein\n"
            "##FASTA\n"
            ">NZ_1\nATGAAACCCATGAAACCC\n",
            encoding="utf-8",
        )
        return gff

    def test_combined_gff_writes_gene_and_product_metadata(self):
        with tempfile.TemporaryDirectory() as td:
            self._write_combined_gff(td)
            dna = os.path.join(td, "genes_dna.fasta")
            annotations = os.path.join(td, "gene_annotations.tsv.gz")
            read_combined_files(
                td, ".gff3", "CDS", dna, False, False,
                annotation_out=annotations,
            )
            records = list(iter_gene_annotations(annotations))
            self.assertEqual(len(records), 2)
            by_id = {row["sequence_id"]: row for row in records}
            first = by_id["GCA_000001.1|geneA"]
            self.assertEqual(first["gene_name"], "abcA")
            self.assertEqual(first["annotation"], "Alpha enzyme, subunit A")
            second = by_id["GCA_000001.1|geneB"]
            self.assertEqual(second["gene_name"], "")
            self.assertEqual(second["annotation"], "Hypothetical protein")
            with open(dna, encoding="utf-8") as handle:
                self.assertEqual(sum(line.startswith(">") for line in handle), 2)

    def test_invalid_cds_does_not_replace_existing_fasta_or_annotation(self):
        with tempfile.TemporaryDirectory() as td:
            gff = Path(td) / "bad.gff3"
            gff.write_text(
                "##gff-version 3\n"
                "ctg\tBakta\tCDS\t1\t3\t.\t+\t0\tID=bad;gene=bad;product=bad\n"
                "##FASTA\n>ctg\nATX\n",
                encoding="utf-8",
            )
            dna = Path(td) / "genes_dna.fasta"
            annotation = Path(td) / "gene_annotations.tsv.gz"
            dna.write_text("PREVIOUS_DNA\n", encoding="utf-8")
            annotation.write_bytes(b"PREVIOUS_ANNOTATION")
            with self.assertRaises(ValueError):
                read_combined_files(
                    td, ".gff3", "CDS", str(dna), False, False,
                    annotation_out=str(annotation),
                )
            self.assertEqual(dna.read_text(encoding="utf-8"), "PREVIOUS_DNA\n")
            self.assertEqual(annotation.read_bytes(), b"PREVIOUS_ANNOTATION")

    def test_combined_gff_file_order_is_deterministic(self):
        with tempfile.TemporaryDirectory() as td:
            # Create in reverse lexical order; extraction should still sort by basename.
            for genome in ("GCA_000002.1", "GCA_000001.1"):
                Path(td, genome + ".gff3").write_text(
                    "##gff-version 3\n"
                    "ctg\tBakta\tCDS\t1\t3\t.\t+\t0\tID=gene1;gene=x;product=Protein\n"
                    "##FASTA\n>ctg\nATG\n",
                    encoding="utf-8",
                )
            dna = os.path.join(td, "ordered_dna.fasta")
            read_combined_files(td, ".gff3", "CDS", dna, False, False)
            with open(dna, encoding="utf-8") as handle:
                first_header = next(line.strip() for line in handle if line.startswith(">"))
            self.assertEqual(first_header, ">GCA_000001.1|gene1")

    def test_separate_gff_fasta_writes_bakta_prokka_metadata(self):
        with tempfile.TemporaryDirectory() as td:
            Path(td, "GCA_000003.1.gff3").write_text(
                "##gff-version 3\n"
                "contig1\tProkka\tCDS\t1\t9\t.\t-\t0\t"
                "ID=cds1;gene=abcX;product=ABC transporter subunit X\n",
                encoding="utf-8",
            )
            Path(td, "GCA_000003.1.fna").write_text(
                ">contig1\nGGGTTTCAT\n", encoding="utf-8"
            )
            dna = os.path.join(td, "paired_dna.fasta")
            annotation = os.path.join(td, "gene_annotations.tsv.gz")
            read_separate_files(
                td, ".gff3", ".fna", "CDS", dna, False, False,
                annotation_out=annotation,
            )
            records = list(iter_gene_annotations(annotation))
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0]["sequence_id"], "GCA_000003.1|cds1")
            self.assertEqual(records[0]["gene_name"], "abcX")
            self.assertEqual(records[0]["annotation"], "ABC transporter subunit X")
            with open(dna, encoding="utf-8") as handle:
                self.assertIn("ATGAAACCC", "".join(line.strip() for line in handle if not line.startswith(">")))


class IncrementalNamingTests(unittest.TestCase):
    def test_incremental_outputs_use_stage_specific_names(self):
        with tempfile.TemporaryDirectory() as td:
            original = Path(td) / "original.fasta"
            original.write_text(
                ">GCA_000001.1|gene1\nATGAAACCC\n"
                ">GCA_000001.1|gene2\nATGCCCGGG\n",
                encoding="utf-8",
            )
            new = Path(td) / "new.fasta"
            new.write_text(">GCA_000002.1|gene1\nATGAAACCC\n", encoding="utf-8")
            clusters = Path(td) / "original.clstr"
            clusters.write_text(
                ">Cluster 0\n"
                "0 9nt, >GCA_000001.1|gene1... *\n"
                ">Cluster 1\n"
                "0 9nt, >GCA_000001.1|gene2... *\n",
                encoding="utf-8",
            )
            result = prepare_incremental_inputs(
                str(clusters), str(original), str(new), str(Path(td) / "out")
            )
            self.assertTrue(result["new_genes_fasta"].endswith("supplemental_new_genes_dna.fasta"))
            self.assertTrue(result["combined_fasta"].endswith("supplemented_full_sequences_dna.fasta"))
            self.assertTrue(result["representatives_dna_fasta"].endswith("reclustering_input_dna.fasta"))
            self.assertEqual(result["original_sequence_count"], 2)
            self.assertEqual(result["new_gene_count"], 1)


class GPAMetadataTests(unittest.TestCase):
    def test_family_metadata_consensus_is_deterministic(self):
        self.assertEqual(_family_metadata_consensus(["abc", "xyz", "abc", ""]), "abc")
        self.assertEqual(_family_metadata_consensus(["xyz", "abc"]), "abc")
        self.assertEqual(_family_metadata_consensus(["", ""]), "")

    def _metadata(self, path):
        with open(path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sequence_id", "gene_name", "annotation", "source_gff"])
            writer.writerow(["GCA_000001.1|gene1", "abcA", "Alpha enzyme", "a.gff"])
            writer.writerow(["GCA_000002.1|gene1", "abcA", "Alpha enzyme", "b.gff"])

    def test_gpa_includes_metadata_for_species_and_genus(self):
        class Family:
            family_id = "group_test"
            prevalence = 2
            members = ["GCA_000001.1|gene1", "GCA_000002.1|gene1"]

        class SequenceStore:
            def sequence_length(self, identifier):
                return 9

        class Resolver:
            sequence_delimiter = "|"
            def genome(self, identifier):
                return identifier.split("|", 1)[0]
            def taxon(self, identifier):
                return "Streptococcus" if identifier.startswith("GCA_000001") else "Enterococcus"

        class Ledger:
            genomes = ("GCA_000001.1", "GCA_000002.1")
            taxa = ("Enterococcus", "Streptococcus")
            resolver = Resolver()
            sequences = SequenceStore()
            def sorted_families(self, _family_ids):
                return [Family()]

        with tempfile.TemporaryDirectory() as td:
            metadata = os.path.join(td, "gene_annotations.tsv")
            self._metadata(metadata)
            for mode in ("Species", "Genus"):
                output = os.path.join(td, f"{mode}.csv")
                _write_gpa(output, Ledger(), ["group_test"], mode, [metadata])
                with open(output, encoding="utf-8", newline="") as handle:
                    rows = list(csv.reader(handle))
                self.assertEqual(rows[0][:len(GPA_METADATA_COLUMNS)], GPA_METADATA_COLUMNS)
                self.assertEqual(rows[1][0], "group_test")
                self.assertEqual(rows[1][1], "abcA")
                self.assertEqual(rows[1][2], "Alpha enzyme")

    def test_shared_workflow_species_and_genus_smoke(self):
        with tempfile.TemporaryDirectory() as td:
            fasta = Path(td) / "genes.fasta"
            fasta.write_text(
                ">GCA_000001.1|gene1\nATGAAACCC\n"
                ">GCA_000002.1|gene1\nATGAAACCC\n",
                encoding="utf-8",
            )
            clstr = Path(td) / "clusters.clstr"
            clstr.write_text(
                ">Cluster 0\n"
                "0 9nt, >GCA_000001.1|gene1... *\n"
                "1 9nt, >GCA_000002.1|gene1... at 100.00%\n",
                encoding="utf-8",
            )
            annotations = Path(td) / "gene_annotations.tsv"
            self._metadata(annotations)
            taxon_map = Path(td) / "taxon_map.tsv"
            taxon_map.write_text(
                "identifier\tgenome\ttaxon\n"
                "GCA_000001.1\tGCA_000001.1\tStreptococcus\n"
                "GCA_000002.1\tGCA_000002.1\tEnterococcus\n",
                encoding="utf-8",
            )

            for mode in ("Species", "Genus"):
                outdir = Path(td) / mode
                options = SimpleNamespace(
                    fasta=str(fasta), clusters=str(clstr), cluster_format="CD-HIT",
                    reclustered=None, sequence_tag="StORF", taxon_map=str(taxon_map) if mode == "Genus" else None,
                    sequence_delimiter="|", mem=100, sequence_buffer_cap_mb=4, tmp_dir=None,
                    sequence_db=None, sequence_db_out=None, output_dir=str(outdir),
                    species_groups="99,95,15,0", core_groups="1,2,3,4,5,10",
                    gene_presence_absence_out=True, single_copy_only=False,
                    single_copy_tolerance=0.0, gene_annotation_files=[str(annotations)],
                    gene_group_table=False, core_group_table=False, write_groups=None,
                    write_individual_groups=False, align_core=False, align_aa=False,
                    visualise=False, verbose=False, added_genomes=None,
                )
                run_pangenome_workflow(options, mode)
                gpa = outdir / "gene_presence_absence.csv"
                self.assertTrue(gpa.is_file())
                with open(gpa, encoding="utf-8", newline="") as handle:
                    rows = list(csv.reader(handle))
                self.assertEqual(rows[1][1], "abcA")
                self.assertEqual(rows[1][2], "Alpha enzyme")


if __name__ == "__main__":
    unittest.main()
