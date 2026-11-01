import csv

import pytest

from PyamilySeq.order_analysis import Genome, build_reference_reports, main, parse_gff


def cohort(tmp_path, names=("A", "B")):
    gd, fd = tmp_path / "gff", tmp_path / "fasta"
    gd.mkdir(); fd.mkdir()
    genomes = []
    metadata = []
    for name in names:
        (gd / f"{name}.gff3").write_text(
            "ctg\tRef\tregion\t1\t100\t.\t+\t.\tID=chromosome;gbkey=Src;Is_circular=true\n"
            "ctg\tRef\tCDS\t1\t90\t.\t+\t0\tID=g1\n")
        (fd / f"{name}.fa").write_text(">ctg chromosome\n" + "A" * 100 + "\n")
        genomes.append(parse_gff(gd / f"{name}.gff3"))
        metadata.append({"genome": name, "contig": "ctg", "replicon_id": "chromosome_1",
                         "replicon_type": "chromosome", "replicon_topology": "circular",
                         "conflict_status": "none"})
    return genomes, metadata, fd


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_two_subcohorts_explicit_references_and_separate_replicons(tmp_path):
    genomes, metadata, fd = cohort(tmp_path, ("A", "B", "C", "D"))
    (tmp_path / "groups.tsv").write_text("genome\tsubcohort\nA\tg1\nB\tg1\nC\tg2\nD\tg2\n")
    (tmp_path / "refs.tsv").write_text("subcohort\treference_genome\ng1\tA\ng2\tC\n")
    # A carries a separate plasmid; its coordinate identity remains a separate manifest entry.
    metadata.append({"genome": "A", "contig": "p", "replicon_id": "plasmid_p", "replicon_type": "plasmid",
                     "replicon_topology": "circular", "conflict_status": "none"})
    result = build_reference_reports(genomes, metadata, tmp_path / "out", fd,
        tmp_path / "groups.tsv", tmp_path / "refs.tsv")
    assert result["selected_references"] == {"g1": "A", "g2": "C"}
    candidates = read(tmp_path / "out" / "reference_candidates.tsv")
    manifest = read(tmp_path / "out" / "reference_manifest.tsv")
    assert len(candidates) == 4
    assert {r["selected_reference_genome"] for r in manifest} == {"A", "C"}
    assert '"plasmid_p"' in next(r for r in manifest if r["subcohort"] == "g1")["reference_replicons"]


def test_default_cohort_has_candidates_but_no_implicit_reference(tmp_path):
    genomes, metadata, fd = cohort(tmp_path)
    build_reference_reports(genomes, metadata, tmp_path / "out", fd)
    rows = read(tmp_path / "out" / "reference_manifest.tsv")
    assert len(rows) == 1
    assert rows[0]["selected_reference_genome"] == ""
    assert rows[0]["selection_method"] == "not_selected"
    assert "No reference was specified" in rows[0]["selection_warnings"]


def test_missing_and_duplicate_subcohort_assignments_are_rejected(tmp_path):
    genomes, metadata, fd = cohort(tmp_path)
    mapping = tmp_path / "groups.tsv"
    mapping.write_text("genome\tsubcohort\nA\tx\n")
    with pytest.raises(ValueError, match="missing genome assignments"):
        build_reference_reports(genomes, metadata, tmp_path / "out", fd, mapping)
    mapping.write_text("genome\tsubcohort\nA\tx\nA\ty\nB\tx\n")
    with pytest.raises(ValueError, match="Duplicate sub-cohort"):
        build_reference_reports(genomes, metadata, tmp_path / "out", fd, mapping)


def test_unknown_genome_and_empty_subcohort_rejected(tmp_path):
    genomes, metadata, fd = cohort(tmp_path)
    mapping = tmp_path / "groups.tsv"
    mapping.write_text("genome\tsubcohort\nA\tx\nB\t\n")
    with pytest.raises(ValueError, match="empty genome or subcohort"):
        build_reference_reports(genomes, metadata, tmp_path / "out", fd, mapping)
    mapping.write_text("genome\tsubcohort\nA\tx\nZ\ty\n")
    with pytest.raises(ValueError, match="unknown genome"):
        build_reference_reports(genomes, metadata, tmp_path / "out", fd, mapping)


def test_reference_absent_or_wrong_subcohort_rejected(tmp_path):
    genomes, metadata, fd = cohort(tmp_path)
    groups, refs = tmp_path / "groups.tsv", tmp_path / "refs.tsv"
    groups.write_text("genome\tsubcohort\nA\tx\nB\ty\n")
    refs.write_text("subcohort\treference_genome\nx\tZ\ny\tB\n")
    with pytest.raises(ValueError, match="absent from"):
        build_reference_reports(genomes, metadata, tmp_path / "out", fd, groups, refs)
    refs.write_text("subcohort\treference_genome\nx\tB\ny\tA\n")
    with pytest.raises(ValueError, match="assigned to"):
        build_reference_reports(genomes, metadata, tmp_path / "out", fd, groups, refs)


def test_reference_missing_fasta_or_gff_annotation_rejected(tmp_path):
    genomes, metadata, fd = cohort(tmp_path)
    refs = tmp_path / "refs.tsv"
    refs.write_text("subcohort\treference_genome\ndefault\tA\n")
    (fd / "A.fa").unlink()
    with pytest.raises(ValueError, match="matching FASTA"):
        build_reference_reports(genomes, metadata, tmp_path / "out", fd, reference_map=refs)
    (fd / "A.fa").write_text(">ctg\nAAAA\n")
    no_annotation = Genome(name="A", path=tmp_path / "A.gff3")
    with pytest.raises(ValueError, match="no usable GFF annotation"):
        build_reference_reports([no_annotation, genomes[1]], metadata, tmp_path / "out", fd, reference_map=refs)


def test_missing_assembly_metadata_is_reported_without_completeness_inference(tmp_path):
    genomes, metadata, fd = cohort(tmp_path)
    build_reference_reports(genomes, metadata, tmp_path / "out", fd)
    rows = read(tmp_path / "out" / "reference_candidates.tsv")
    assert all(r["assembly_level"] == "" and r["completeness"] == "" for r in rows)
    assert all("assembly level not supplied" in r["reference_suitability_warnings"] for r in rows)
    assert all("completeness metadata not supplied" in r["reference_suitability_warnings"] for r in rows)


def test_candidate_ranking_is_deterministic_and_not_an_auto_selection(tmp_path):
    genomes, metadata, fd = cohort(tmp_path)
    for _ in range(2):
        out = tmp_path / f"out{_}"
        build_reference_reports(genomes, metadata, out, fd)
    first = read(tmp_path / "out0" / "reference_candidates.tsv")
    second = read(tmp_path / "out1" / "reference_candidates.tsv")
    assert [(r["genome"], r["candidate_review_rank"]) for r in first] == [
        (r["genome"], r["candidate_review_rank"]) for r in second]
    assert all(r["recommendation"].endswith("not selected") for r in first)


def test_reference_map_requires_every_subcohort_exactly_once(tmp_path):
    genomes, metadata, fd = cohort(tmp_path)
    groups, refs = tmp_path / "groups.tsv", tmp_path / "refs.tsv"
    groups.write_text("genome\tsubcohort\nA\tx\nB\ty\n")
    refs.write_text("subcohort\treference_genome\nx\tA\n")
    with pytest.raises(ValueError, match="missing sub-cohort assignments"):
        build_reference_reports(genomes, metadata, tmp_path / "out", fd, groups, refs)


def test_reference_conflict_rejected_and_candidate_warning_retained(tmp_path):
    genomes, metadata, fd = cohort(tmp_path)
    metadata[0]["conflict_status"] = "identity"
    refs = tmp_path / "refs.tsv"
    refs.write_text("subcohort\treference_genome\ndefault\tA\n")
    with pytest.raises(ValueError, match="conflicting replicon metadata"):
        build_reference_reports(genomes, metadata, tmp_path / "out", fd, reference_map=refs)


def test_sourmash_matrix_reports_only_within_subcohort_statistics(tmp_path):
    genomes, metadata, fd = cohort(tmp_path)
    matrix = tmp_path / "similarity.csv"
    matrix.write_text(",A,B\nA,1,0.9\nB,0.9,1\n")
    build_reference_reports(genomes, metadata, tmp_path / "out", fd, similarity_matrix=matrix)
    assert all(r["sourmash_within_subcohort_n"] == "1" for r in read(tmp_path / "out" / "reference_candidates.tsv"))


def test_cli_writes_explicit_reference_outputs(tmp_path):
    genomes, _metadata, fd = cohort(tmp_path)
    (tmp_path / "groups.tsv").write_text("genome\tsubcohort\nA\tnear\nB\tnear\n")
    (tmp_path / "refs.tsv").write_text("subcohort\treference_genome\nnear\tA\n")
    result = main(["--gff-dir", str(tmp_path / "gff"), "--fasta-dir", str(fd),
                   "--subcohort-map", str(tmp_path / "groups.tsv"), "--reference-map", str(tmp_path / "refs.tsv"),
                   "--outdir", str(tmp_path / "cli-out"), "--no-plots"])
    assert result == 0
    candidate_rows = read(tmp_path / "cli-out" / "reference_candidates.tsv")
    manifest_rows = read(tmp_path / "cli-out" / "reference_manifest.tsv")
    assert len(candidate_rows) == 2
    assert manifest_rows[0]["selected_reference_genome"] == "A"
    assert manifest_rows[0]["selection_method"] == "user-specified"

