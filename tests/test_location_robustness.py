import gzip

from PyamilySeq.location_tracking import (
    compute_family_metrics,
    compute_neighbors_for_contigs,
    map_family_members_to_features,
    parse_clusters,
    parse_gff,
)


def test_parent_gene_and_child_cds_collapse_to_one_locus(tmp_path):
    gff = tmp_path / "Genome1.gff3"
    gff.write_text(
        "##sequence-region ctg 1 1000\n"
        "ctg\tRef\tgene\t100\t300\t.\t+\t.\tID=gene1;product=misleading_alias\n"
        "ctg\tRef\tCDS\t120\t280\t.\t+\t0\tID=cds1;Parent=gene1;product=misleading_alias\n"
    )
    contigs, aliases = parse_gff(str(gff))

    assert len(contigs["ctg"]) == 1
    locus = contigs["ctg"][0]
    assert locus["type"] == "CDS"
    assert locus["collapsed_feature_types"] == ["CDS", "gene"]
    assert locus["position_basis"] == "bp"
    assert aliases["gene1"] is locus
    assert aliases["cds1"] is locus
    assert "misleading_alias" not in aliases


def test_gene_transcript_cds_hierarchy_collapses_through_filtered_parent(tmp_path):
    gff = tmp_path / "Genome1.gff3"
    gff.write_text(
        "ctg\tRef\tgene\t100\t400\t.\t+\t.\tID=gene1\n"
        "ctg\tRef\tmRNA\t100\t400\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "ctg\tRef\tCDS\t120\t200\t.\t+\t0\tID=cds1;Parent=tx1\n"
        "ctg\tRef\tCDS\t250\t380\t.\t+\t0\tID=cds1;Parent=tx1\n"
    )

    contigs, aliases = parse_gff(str(gff))

    assert len(contigs["ctg"]) == 1
    locus = contigs["ctg"][0]
    assert (locus["start"], locus["end"], locus["feature_parts"]) == (120, 380, 2)
    assert aliases["gene1"] is locus
    assert aliases["tx1"] is locus
    assert aliases["cds1"] is locus


def test_negative_strand_neighbors_are_transcription_relative(tmp_path):
    gff = tmp_path / "Genome1.gff"
    gff.write_text(
        "ctg\tRef\tCDS\t10\t20\t.\t+\t0\tID=left\n"
        "ctg\tRef\tCDS\t30\t40\t.\t-\t0\tID=target\n"
        "ctg\tRef\tCDS\t50\t60\t.\t-\t0\tID=right\n"
    )
    contigs, aliases = parse_gff(str(gff))
    families = {
        "fam_left": ["left"],
        "fam_target": ["target"],
        "fam_right": ["right"],
    }
    occurrences, _ = map_family_members_to_features(families, aliases)
    compute_neighbors_for_contigs(contigs, k=1)
    target = occurrences["fam_target"][0]
    assert target["neighbors_up"] == ["right"]
    assert target["neighbors_down"] == ["left"]
    assert target["neighbor_families_up"] == ["fam_right"]
    assert target["signed_neighbor_family_context"] == ["-1:fam_right", "+1:fam_left"]


def test_circular_contigs_wrap_without_boundary_censoring(tmp_path):
    gff = tmp_path / "Genome1.gff3"
    gff.write_text(
        "##sequence-region ctg 1 100\n"
        "ctg\tRef\tCDS\t1\t10\t.\t+\t0\tID=first\n"
        "ctg\tRef\tCDS\t40\t50\t.\t+\t0\tID=middle\n"
        "ctg\tRef\tCDS\t80\t90\t.\t+\t0\tID=last\n"
    )
    contigs, _ = parse_gff(str(gff), circular_contigs={"ctg"})
    compute_neighbors_for_contigs(contigs, k=1)
    first = contigs["ctg"][0]
    assert first["neighbors_up"] == ["last"]
    assert first["boundary_censored_up"] is False
    assert first["intergenic_bp_up"] == [10]
    assert first["contig_topology"] == "circular"
    assert first["topology_source"] == "user"


def test_gff_topology_declarations_remain_distinct_from_unspecified(tmp_path):
    gff = tmp_path / "Genome1.gff3"
    gff.write_text(
        "circular\tRef\tregion\t1\t100\t.\t.\t.\tID=circular;Is_circular=true\n"
        "circular\tRef\tCDS\t10\t20\t.\t+\t0\tID=a\n"
        "linear\tRef\tregion\t1\t100\t.\t.\t.\tID=linear;Is_circular=false\n"
        "linear\tRef\tCDS\t10\t20\t.\t+\t0\tID=b\n"
        "unknown\tRef\tCDS\t10\t20\t.\t+\t0\tID=c\n"
    )

    contigs, _ = parse_gff(str(gff))

    assert (contigs["circular"][0]["contig_topology"], contigs["circular"][0]["topology_source"]) == (
        "circular", "gff"
    )
    assert (contigs["linear"][0]["contig_topology"], contigs["linear"][0]["topology_source"]) == (
        "linear", "gff"
    )
    assert (contigs["unknown"][0]["contig_topology"], contigs["unknown"][0]["topology_source"]) == (
        "unspecified", "none"
    )


def test_gzip_gff_and_mapping_confidence_threshold(tmp_path):
    gff = tmp_path / "Genome1.gff3.gz"
    with gzip.open(gff, "wt") as handle:
        handle.write("ctg\tRef\tCDS\t1\t9\t.\t+\t0\tID=abc_def\n")
    _, aliases = parse_gff(str(gff))
    occurrences, unmapped, diagnostics = map_family_members_to_features(
        {"fam": ["abc-def"]},
        aliases,
        return_diagnostics=True,
        min_confidence=0.8,
    )
    assert occurrences["fam"] == []
    assert unmapped["fam"] == ["abc-def"]
    assert diagnostics[0]["status"] == "low_confidence"
    assert diagnostics[0]["confidence"] == 0.6


def test_location_metrics_report_signed_architecture_and_boundary_flags(tmp_path):
    gff = tmp_path / "Genome1.gff"
    gff.write_text(
        "ctg\tRef\tCDS\t10\t20\t.\t+\t0\tID=a\n"
        "ctg\tRef\tCDS\t30\t40\t.\t+\t0\tID=b\n"
    )
    contigs, aliases = parse_gff(str(gff))
    occurrences, _ = map_family_members_to_features({"A": ["a"], "B": ["b"]}, aliases)
    compute_neighbors_for_contigs(contigs, k=1)
    rows, summaries = compute_family_metrics(occurrences, k=1)
    row_a = next(row for row in rows if row["family_id"] == "A")
    summary_a = next(row for row in summaries if row["family_id"] == "A")
    assert row_a["signed_neighbor_family_context"] == "+1:B"
    assert row_a["boundary_censored_up"] is True
    assert summary_a["dominant_signed_family_architecture"] == "+1:B"
    assert summary_a["boundary_censored_fraction"] == 1.0


def test_cluster_csv_is_not_mistaken_for_gpa_and_mmseqs_rows_accumulate(tmp_path):
    csv_path = tmp_path / "families.csv"
    csv_path.write_text("family_id,member_id\nfam1,geneA\nfam1,geneB\n")
    assert parse_clusters(str(csv_path)) == {"fam1": ["geneA", "geneB"]}

    mmseqs = tmp_path / "mmseqs.tsv"
    mmseqs.write_text("rep\tchild1\nother\tmember\nrep\tchild2\n")
    assert parse_clusters(str(mmseqs), cluster_format="MMseqs")["rep"] == [
        "rep",
        "child1",
        "child2",
    ]
