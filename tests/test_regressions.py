import csv
import gzip
import sys
from types import SimpleNamespace

from PyamilySeq.PyamilySeq import main as pyamilyseq_main
from PyamilySeq.PyamilySeq_Species import calc_First_only_core, get_cores
from PyamilySeq.clusterings import cluster_BLAST
from PyamilySeq.location_tracking import parse_clusters
from PyamilySeq.utils import read_combined_files, read_fasta


def test_blast_partial_backfills_fasta_singletons(tmp_path, monkeypatch):
    cluster_file = tmp_path / "blast.tsv"
    fasta_file = tmp_path / "seqs.fa"
    output_dir = tmp_path / "out"

    cluster_file.write_text("G1|geneA\tG2|geneB\n")
    fasta_file.write_text(
        ">G1|geneA\nATGAAATAG\n"
        ">G2|geneB\nATGAAATAG\n"
        ">G3|geneC\nATGCCCTAG\n"
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "PyamilySeq",
            "Partial",
            "-clustering_format",
            "BLAST",
            "-cluster_file",
            str(cluster_file),
            "-original_fasta",
            str(fasta_file),
            "-output_dir",
            str(output_dir),
            "-gene_group_table",
            "-visualise",
            "-visualise_heatmap_limit",
            "10",
        ],
    )
    pyamilyseq_main()

    table = output_dir / "gene_group_table.tsv"
    assignments = dict(line.rstrip().split("\t") for line in table.read_text().splitlines())
    assert assignments["G1|geneA"] == assignments["G2|geneB"]
    assert assignments["G3|geneC"].startswith("group_")
    assert set(assignments) == {"G1|geneA", "G2|geneB", "G3|geneC"}

    with (output_dir / "gene_presence_absence.csv").open(newline="") as handle:
        rows = list(csv.reader(handle))
    genes_in_matrix = set()
    for row in rows[1:]:
        for cell in row[14:]:
            genes_in_matrix.update(cell.split())
    assert genes_in_matrix == {"geneA", "geneB", "geneC"}

    visual_dir = output_dir / "Visualisations"
    expected_visuals = [
        "pangenome_visualisations.html",
        "core_categories.svg",
        "family_prevalence_histogram.svg",
        "presence_absence_heatmap.svg",
        "category_counts.tsv",
        "family_prevalence.tsv",
    ]
    for visual_name in expected_visuals:
        assert (visual_dir / visual_name).exists()
    assert (visual_dir / "core_categories.svg").read_text().startswith("<svg")
    assert "Presence/Absence" in (visual_dir / "pangenome_visualisations.html").read_text()


def test_blast_parser_counts_representative_genomes_without_self_hits(tmp_path):
    cluster_file = tmp_path / "blast.tsv"
    fasta_file = tmp_path / "seqs.fa"
    cluster_file.write_text("G1|geneA\tG2|geneB\n")
    fasta_file.write_text(
        ">G1|geneA\nATGAAATAG\n"
        ">G2|geneB\nATGAAATAG\n"
        ">G3|geneC\nATGCCCTAG\n"
    )

    options = SimpleNamespace(
        clusters=str(cluster_file),
        fasta=str(fasta_file),
        reclustered=None,
        sequence_tag="StORF",
    )
    genome_dict, clusters, _, sequences, _ = cluster_BLAST(options, "|")

    assert set(genome_dict) == {"G1", "G2", "G3"}
    assert sorted(len(genomes) for genomes in clusters.values()) == [1, 2]
    assert {seq for seqs in sequences.values() for seq in seqs} == {
        "G1|geneA",
        "G2|geneB",
        "G3|geneC",
    }


def test_species_core_bins_use_threshold_boundaries():
    genome_dict = {f"G{i}": 1 for i in range(100)}
    options = SimpleNamespace(species_groups="99,95,15,0", reclustered=None)

    for size, expected_key in [
        (100, "First_core_99"),
        (99, "First_core_99"),
        (98, "First_core_95"),
        (95, "First_core_95"),
        (94, "First_core_15"),
        (15, "First_core_15"),
        (14, "First_core_0"),
    ]:
        cores, groups = get_cores(options, genome_dict)
        calc_First_only_core("cluster", size, groups, cores)
        assert cores[expected_key] == ["cluster"]


def test_location_tracking_two_column_cluster_with_comma_members(tmp_path):
    cluster_file = tmp_path / "clusters.txt"
    cluster_file.write_text("fam1\tgeneA,geneX\nfam2\tgeneB,geneY\n")

    assert parse_clusters(str(cluster_file)) == {
        "fam1": ["geneA", "geneX"],
        "fam2": ["geneB", "geneY"],
    }


def test_genus_mode_partial_does_not_require_align(tmp_path, monkeypatch):
    cluster_file = tmp_path / "test.clstr"
    fasta_file = tmp_path / "seqs.fa"
    output_dir = tmp_path / "out"
    cluster_file.write_text(">Cluster 0\n0\t9nt, >Genus_species|geneA... *\n")
    fasta_file.write_text(">Genus_species|geneA\nATGAAATAG\n")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "PyamilySeq",
            "Partial",
            "-clustering_format",
            "CD-HIT",
            "-cluster_file",
            str(cluster_file),
            "-original_fasta",
            str(fasta_file),
            "-output_dir",
            str(output_dir),
            "-group_mode",
            "Genus",
            "--visualize",
            "--visualize-heatmap-limit",
            "5",
        ],
    )
    pyamilyseq_main()

    assert (output_dir / "gene_presence_absence.csv").exists()
    assert (output_dir / "Visualisations" / "pangenome_visualisations.html").exists()


def test_read_fasta_supports_gzip(tmp_path):
    fasta_file = tmp_path / "seqs.fa.gz"
    with gzip.open(fasta_file, "wt") as handle:
        handle.write(">G1|geneA\nATG\n>G2|geneB\nCCC\n")

    assert read_fasta(str(fasta_file)) == {
        "G1|geneA": "ATG",
        "G2|geneB": "CCC",
    }


def test_combined_gff_extraction_handles_gene_ident_duplicates_and_negative_bounds(tmp_path):
    gff_file = tmp_path / "Genome1.gff3"
    output_fasta = tmp_path / "combined.fasta"
    gff_file.write_text(
        "##gff-version 3\n"
        "ctg\tRef\tCDS\t1\t6\t.\t-\t0\tID=dup\n"
        "ctg\tRef\ttRNA\t1\t3\t.\t+\t0\tID=dup\n"
        "ctg\tRef\tgene\t4\t6\t.\t+\t0\tID=skip\n"
        "##FASTA\n"
        ">ctg\n"
        "AAACCC\n"
    )

    read_combined_files(
        str(tmp_path),
        ".gff3",
        "CDS,tRNA",
        str(output_fasta),
        translate=False,
        run_as_combiner=True,
    )

    assert read_fasta(str(output_fasta)) == {
        "Genome1|dup": "GGGTTT",
        "Genome1|dup_1": "AAA",
    }
