from PyamilySeq.location_tracking import (
    compute_family_metrics,
    compute_neighbors_for_contigs,
    map_family_members_to_features,
    parse_clusters,
    parse_gff,
    parse_gff_dir,
    write_mapping_diagnostics,
    write_location_visualisations,
    write_run_report,
)


def test_location_tracking_small(tmp_path):
    # create a tiny GFF
    gff_text = """##gff-version 3
contig1	Ref	CDS	100	200	.	+	0	ID=geneA
contig1	Ref	CDS	300	400	.	+	0	ID=geneB
contig1	Ref	CDS	500	600	.	+	0	ID=geneC
contig2	Ref	CDS	50	150	.	-	0	ID=geneX
contig2	Ref	CDS	200	300	.	-	0	ID=geneY
"""
    gff_file = tmp_path / "test.gff"
    gff_file.write_text(gff_text)

    # simple cluster file: family -> members (members use gene IDs)
    cluster_text = """fam1\tgeneA,geneX\nfam2\tgeneB,geneY\n"""
    cluster_file = tmp_path / "clusters.txt"
    cluster_file.write_text(cluster_text)

    contigs, id_map = parse_gff(str(gff_file))
    families = parse_clusters(str(cluster_file))
    fam_occ, unmapped = map_family_members_to_features(families, id_map)
    compute_neighbors_for_contigs(contigs, k=1)
    occ_rows, sum_rows = compute_family_metrics(fam_occ, k=1, rel_delta=0.5)

    # expectations: fam1 has 2 occurrences, fam2 has 2
    fams = {r['family_id']: r for r in sum_rows}
    assert 'fam1' in fams and fams['fam1']['occurrences'] == 2
    assert 'fam2' in fams and fams['fam2']['occurrences'] == 2

    # occurrences rows should include geneA and geneX
    found_ids = set([r['gene_id'] for r in occ_rows])
    assert 'geneA' in found_ids and 'geneX' in found_ids


def test_location_tracking_refuses_prefix_false_matches(tmp_path):
    gff_file = tmp_path / "test.gff"
    gff_file.write_text("ctg\tRef\tCDS\t1\t9\t.\t+\t0\tID=gene10\n")
    families = {"fam1": ["gene1"]}

    _, id_map = parse_gff(str(gff_file))
    fam_occ, unmapped, diagnostics = map_family_members_to_features(
        families,
        id_map,
        return_diagnostics=True,
    )

    assert fam_occ["fam1"] == []
    assert unmapped["fam1"] == ["gene1"]
    assert diagnostics[0]["status"] == "unmapped"


def test_location_tracking_genome_qualified_ids_disambiguate_duplicates(tmp_path):
    gff_dir = tmp_path / "gffs"
    gff_dir.mkdir()
    (gff_dir / "Genome1.gff3").write_text("ctg\tRef\tCDS\t10\t20\t.\t+\t0\tID=geneA\n")
    (gff_dir / "Genome2.gff3").write_text("ctg\tRef\tCDS\t30\t40\t.\t-\t0\tID=geneA\n")

    contigs, id_map, _ = parse_gff_dir(str(gff_dir))
    qualified = {"famA": ["Genome1|geneA", "Genome2|geneA"]}
    fam_occ, unmapped, diagnostics = map_family_members_to_features(
        qualified,
        id_map,
        return_diagnostics=True,
    )
    assert len(fam_occ["famA"]) == 2
    assert not unmapped
    assert {row["genome"] for row in diagnostics} == {"Genome1", "Genome2"}

    ambiguous = {"famA": ["geneA"]}
    fam_occ, unmapped, diagnostics = map_family_members_to_features(
        ambiguous,
        id_map,
        return_diagnostics=True,
    )
    assert fam_occ["famA"] == []
    assert unmapped["famA"] == ["geneA"]
    assert diagnostics[0]["status"] == "ambiguous"


def test_location_tracking_decodes_gff3_ids_and_reports_mapping(tmp_path):
    gff_file = tmp_path / "encoded.gff3"
    outdir = tmp_path / "out"
    gff_file.write_text("ctg\tRef\tCDS\t5\t15\t.\t+\t0\tID=gene%3Aalpha;Name=alpha\n")

    contigs, id_map = parse_gff(str(gff_file))
    fam_occ, unmapped, diagnostics = map_family_members_to_features(
        {"fam_encoded": ["gene:alpha"]},
        id_map,
        return_diagnostics=True,
    )

    assert len(fam_occ["fam_encoded"]) == 1
    assert not unmapped
    assert diagnostics[0]["status"] == "mapped"
    mapping_path, problem_path = write_mapping_diagnostics(str(outdir), diagnostics)
    report_path = write_run_report(str(outdir), {"fam_encoded": ["gene:alpha"]}, fam_occ, diagnostics, contigs)
    assert "fam_encoded" in (outdir / "member_mapping.tsv").read_text()
    assert (outdir / "unmapped_members.tsv").read_text().splitlines() == [
        "\t".join([
            "family_id",
            "member_id",
            "status",
            "match_type",
            "matched_key",
            "gene_id",
            "genome",
            "contig",
            "start",
            "end",
            "strand",
            "note",
        ])
    ]
    assert mapping_path.endswith("member_mapping.tsv")
    assert problem_path.endswith("unmapped_members.tsv")
    assert "member_mapping_rate\t1.000000" in (outdir / "location_tracking_report.txt").read_text()
    assert report_path.endswith("location_tracking_report.txt")


def test_location_tracking_reports_neighbor_family_context(tmp_path):
    gff_dir = tmp_path / "gffs"
    gff_dir.mkdir()
    (gff_dir / "Genome1.gff3").write_text(
        "ctg\tRef\tCDS\t10\t20\t.\t+\t0\tID=left1\n"
        "ctg\tRef\tCDS\t30\t40\t.\t+\t0\tID=target1\n"
        "ctg\tRef\tCDS\t50\t60\t.\t+\t0\tID=right1\n"
    )
    (gff_dir / "Genome2.gff3").write_text(
        "ctg\tRef\tCDS\t10\t20\t.\t+\t0\tID=left2\n"
        "ctg\tRef\tCDS\t30\t40\t.\t+\t0\tID=target2\n"
        "ctg\tRef\tCDS\t50\t60\t.\t+\t0\tID=right2\n"
    )
    families = {
        "fam_left": ["Genome1|left1", "Genome2|left2"],
        "fam_target": ["Genome1|target1", "Genome2|target2"],
        "fam_right": ["Genome1|right1", "Genome2|right2"],
    }

    contigs, id_map, _ = parse_gff_dir(str(gff_dir))
    fam_occ, unmapped = map_family_members_to_features(families, id_map)
    assert not unmapped
    compute_neighbors_for_contigs(contigs, k=1)
    occ_rows, sum_rows = compute_family_metrics(fam_occ, k=1, rel_delta=0.1)

    target_rows = [row for row in occ_rows if row["family_id"] == "fam_target"]
    assert {row["neighbor_families_up"] for row in target_rows} == {"fam_left"}
    assert {row["neighbor_families_down"] for row in target_rows} == {"fam_right"}
    target_summary = {row["family_id"]: row for row in sum_rows}["fam_target"]
    assert target_summary["neighbor_family_conserved_fraction"] == 1.0
    assert target_summary["dominant_neighbor_family_context"] == "fam_left;fam_right"

    html_path, svg_path, neighbor_summary_path = write_location_visualisations(
        str(tmp_path / "visuals"),
        occ_rows,
        sum_rows,
        family_limit=3,
        rows_per_family=3,
    )
    assert html_path.endswith("location_tracking.html")
    assert svg_path.endswith("family_neighborhoods.svg")
    assert neighbor_summary_path.endswith("neighbor_family_summary.tsv")
    assert "fam_target" in (tmp_path / "visuals" / "family_neighborhoods.svg").read_text()
    assert "neighbor_family_summary.tsv" in (tmp_path / "visuals" / "location_tracking.html").read_text()
