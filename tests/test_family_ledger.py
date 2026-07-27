import csv
import json
from types import SimpleNamespace

import pytest

from PyamilySeq.clusterings import cluster_CDHIT
from PyamilySeq.families import (
    FamilyInputError,
    build_family_ledger,
    parse_blast,
    parse_mmseqs,
    read_fasta_records,
)
from PyamilySeq.workflow import _write_alignment, classify_families, run_pangenome_workflow


def options(tmp_path, fasta, clusters, reclustered=None, cluster_format="MMseqs"):
    return SimpleNamespace(
        fasta=str(fasta),
        clusters=str(clusters),
        reclustered=str(reclustered) if reclustered else None,
        cluster_format=cluster_format,
        sequence_tag="StORF",
        taxon_map=None,
        sequence_delimiter="|",
        species_groups="95,99,15,0",
        core_groups="95,99,15,0",
        single_copy_only=False,
        single_copy_tolerance=0.0,
        output_dir=str(tmp_path / "out"),
        gene_presence_absence_out=True,
        gene_group_table=True,
        core_group_table=True,
        write_groups=None,
        write_individual_groups=False,
        align_core=False,
        align_aa=False,
        sequence_type=None,
        threads=1,
        visualise=True,
        visualise_heatmap_limit=20,
        visualise_permutations=10,
        visualise_seed=7,
        blast_graph_mode="reciprocal-best",
        blast_min_identity=0.0,
        blast_max_evalue=10.0,
        blast_min_query_coverage=0.0,
        blast_min_subject_coverage=0.0,
        blast_columns="qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore",
    )


def test_reclustered_ledger_assigns_every_sequence_once_and_drives_all_outputs(tmp_path):
    fasta = tmp_path / "all.fa"
    initial = tmp_path / "initial.tsv"
    reclustered = tmp_path / "reclustered.tsv"
    fasta.write_text(
        ">G1|a description\nATGAAA\n"
        ">G2|b\nATGCCC\n"
        ">G3|c\nATGGGG\n"
        ">G1|StORF_x\nATGTTT\n"
        ">G4|StORF_y\nATGATA\n"
    )
    initial.write_text("G1|a\tG2|b\nG3|c\tG3|c\n")
    # No self-hit is required for G1|a; repeated representative rows need not be adjacent.
    reclustered.write_text(
        "G1|a\tG1|StORF_x\n"
        "G4|StORF_y\tG4|StORF_y\n"
        "G3|c\tG3|c\n"
    )
    opts = options(tmp_path, fasta, initial, reclustered)

    ledger = run_pangenome_workflow(opts, "Species")

    assert len(ledger.assignment) == 5
    assert set(ledger.assignment) == set(read_fasta_records(str(fasta)))
    assert sorted(family.stage for family in ledger.families.values()) == [
        "First",
        "extended",
        "only_Second",
    ]
    assert sorted(len(family.members) for family in ledger.families.values()) == [1, 1, 3]

    output = tmp_path / "out"
    assignments = list(csv.DictReader((output / "gene_family_assignments.tsv").open(), delimiter="\t"))
    assert len(assignments) == 5
    assert {row["assignment_status"] for row in assignments} == {"assigned"}
    assert sum(
        int(row["sequence_count"])
        for row in csv.DictReader((output / "family_summary.tsv").open(), delimiter="\t")
    ) == 5
    with (output / "gene_presence_absence.csv").open(newline="") as handle:
        gpa = list(csv.reader(handle))
    assert sum(int(row[4]) for row in gpa[1:]) == 5
    assert sum(1 for line in (output / "pan_genome_reference.fa").read_text().splitlines() if line.startswith(">")) == 3
    prevalence = list(csv.DictReader(
        (output / "Visualisations" / "family_prevalence.tsv").open(), delimiter="\t"
    ))
    assert sum(int(row["sequence_count"]) for row in prevalence) == 5
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["exactly_once_assignment"] is True


def test_recluster_provenance_uses_first_round_membership_not_name_tags(tmp_path):
    fasta = tmp_path / "all.fa"
    initial = tmp_path / "initial.tsv"
    reclustered = tmp_path / "reclustered.tsv"
    fasta.write_text(">G1|old\nATG\n>G2|novel_without_tag\nCCC\n>G3|another_novel\nGGG\n")
    initial.write_text("G1|old\tG1|old\n")
    reclustered.write_text(
        "G1|old\tG2|novel_without_tag\n"
        "G3|another_novel\tG3|another_novel\n"
    )

    ledger = build_family_ledger(
        str(fasta), str(initial), "MMseqs", recluster_path=str(reclustered), sequence_tag="StORF"
    )

    assert sorted(family.stage for family in ledger.families.values()) == ["extended", "only_Second"]


def test_duplicate_fasta_identifiers_fail_instead_of_overwriting(tmp_path):
    fasta = tmp_path / "duplicate.fa"
    fasta.write_text(">G1|a first\nAAA\n>G1|a second\nCCC\n")
    with pytest.raises(FamilyInputError, match="Duplicate FASTA identifier"):
        read_fasta_records(str(fasta))


def test_mmseqs_rows_are_aggregated_independent_of_row_order(tmp_path):
    path = tmp_path / "clusters.tsv"
    path.write_text("rep\tchild1\nother\tmember\nrep\tchild2\n")
    parsed = {cluster.label: set(cluster.members) for cluster in parse_mmseqs(str(path))}
    assert parsed["rep"] == {"rep", "child1", "child2"}


def test_scored_blast_uses_reciprocal_best_hits_to_avoid_weak_bridges(tmp_path):
    path = tmp_path / "blast.tsv"

    def row(query, subject, bitscore):
        return f"{query}\t{subject}\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t{bitscore}\n"

    path.write_text(
        row("A", "B", 100)
        + row("B", "A", 100)
        + row("C", "D", 100)
        + row("D", "C", 100)
        + row("B", "C", 1)
        + row("C", "B", 1)
    )
    opts = SimpleNamespace(
        blast_graph_mode="reciprocal-best",
        blast_min_identity=0,
        blast_max_evalue=10,
        blast_min_query_coverage=0,
        blast_min_subject_coverage=0,
        blast_columns="qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore",
    )
    families = {frozenset(cluster.members) for cluster in parse_blast(str(path), opts)}
    assert families == {frozenset({"A", "B"}), frozenset({"C", "D"})}


def test_family_ids_are_deterministic_across_edge_row_order(tmp_path):
    fasta = tmp_path / "seqs.fa"
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    fasta.write_text(">G1|a\nAAA\n>G2|b\nCCC\n>G3|c\nGGG\n")
    first.write_text("G1|a\tG2|b\nG3|c\tG3|c\n")
    second.write_text("G3|c\tG3|c\nG1|a\tG2|b\n")
    opts1 = options(tmp_path, fasta, first)
    opts2 = options(tmp_path, fasta, second)
    ids1 = set(build_family_ledger(str(fasta), str(first), "MMseqs", options=opts1).families)
    ids2 = set(build_family_ledger(str(fasta), str(second), "MMseqs", options=opts2).families)
    assert ids1 == ids2


def test_custom_species_groups_are_validated_and_sorted(tmp_path):
    fasta = tmp_path / "seqs.fa"
    clusters = tmp_path / "clusters.tsv"
    fasta.write_text(">G1|a\nAAA\n")
    clusters.write_text("G1|a\tG1|a\n")
    ledger = build_family_ledger(str(fasta), str(clusters), "MMseqs")
    categories, groups = classify_families(ledger, "Species", "15,99,95")
    assert groups == [99, 95, 15, 0]
    assert categories["First_core_99"]
    with pytest.raises(FamilyInputError):
        classify_families(ledger, "Species", "101,95")


def test_genus_group_ranges_are_expanded(tmp_path):
    fasta = tmp_path / "seqs.fa"
    clusters = tmp_path / "clusters.tsv"
    fasta.write_text(">Genus_one|a\nAAA\n")
    clusters.write_text("Genus_one|a\tGenus_one|a\n")
    ledger = build_family_ledger(str(fasta), str(clusters), "MMseqs", mode="Genus")

    categories, groups = classify_families(ledger, "Genus", "1-3,5")

    assert groups == [1, 2, 3, 5]
    assert categories["First_genera_1"]


def test_protein_group_exports_are_not_retranslated_or_mislabeled(tmp_path):
    fasta = tmp_path / "proteins.fa"
    clusters = tmp_path / "clusters.tsv"
    fasta.write_text(">G1|protein\nMKTWYV\n")
    clusters.write_text("G1|protein\tG1|protein\n")
    opts = options(tmp_path, fasta, clusters)
    opts.sequence_type = "AA"
    opts.write_groups = "99"

    run_pangenome_workflow(opts, "Species")

    output = tmp_path / "out" / "Gene_Groups_Output"
    protein_output = output / "combined_group_sequences_99_aa.fasta"
    assert protein_output.exists()
    assert "MKTWYV" in protein_output.read_text()
    assert not (output / "combined_group_sequences_99_dna.fasta").exists()


def test_core_alignment_rejects_multicopy_families_before_mafft(tmp_path):
    fasta = tmp_path / "seqs.fa"
    clusters = tmp_path / "clusters.tsv"
    fasta.write_text(">G1|a\nATG\n>G1|b\nATG\n>G2|c\nATG\n")
    clusters.write_text("G1|a\tG1|b\nG1|a\tG2|c\n")
    opts = options(tmp_path, fasta, clusters)
    ledger = build_family_ledger(str(fasta), str(clusters), "MMseqs")

    with pytest.raises(FamilyInputError, match="exactly one sequence from every genome"):
        _write_alignment(opts, ledger, list(ledger.families))


def test_legacy_cdhit_adapter_finalizes_last_representative(tmp_path):
    cluster_file = tmp_path / "test.clstr"
    fasta = tmp_path / "seqs.fa"
    cluster_file.write_text(
        ">Cluster 0\n0\t3nt, >G1|a... *\n"
        ">Cluster 1\n0\t3nt, >G2|b... *\n1\t3nt, >G3|c... at 100.00%\n"
    )
    fasta.write_text(">G1|a\nAAA\n>G2|b\nCCC\n>G3|c\nGGG\n")
    opts = SimpleNamespace(clusters=str(cluster_file), fasta=str(fasta), reclustered=None)
    _, _, _, _, reps = cluster_CDHIT(opts, "|")
    assert reps["G2|b"] == [2, 2]
