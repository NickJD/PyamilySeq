from PyamilySeq.validation import validate_inputs


def test_validate_inputs_reports_fasta_and_gff_mapping_issues(tmp_path):
    cluster_file = tmp_path / "clusters.tsv"
    fasta_file = tmp_path / "seqs.fa"
    gff_dir = tmp_path / "gffs"
    outdir = tmp_path / "validation"

    gff_dir.mkdir()
    cluster_file.write_text(
        "fam1\tGenome1|geneA,Genome1|missing\n"
        "fam2\tGenome1|geneA\n"
    )
    fasta_file.write_text(">Genome1|geneA\nATG\n>Genome1|extra\nATG\n")
    (gff_dir / "Genome1.gff3").write_text("ctg\tRef\tCDS\t1\t9\t.\t+\t0\tID=geneA\n")

    summary, issues, paths = validate_inputs(
        str(cluster_file),
        str(outdir),
        fasta=str(fasta_file),
        gff_dir=str(gff_dir),
    )

    codes = {issue["code"] for issue in issues}
    assert "MEMBER_IN_MULTIPLE_FAMILIES" in codes
    assert "CLUSTER_MEMBERS_MISSING_FROM_FASTA" in codes
    assert "GFF_MAPPING_UNMAPPED" in codes
    assert summary["families"] == 2
    assert summary["member_mapping_mapped"] == 2
    assert summary["member_mapping_unmapped"] == 1

    for path in paths:
        assert path
    assert (outdir / "validation_report.json").exists()
    assert "GFF_MAPPING_UNMAPPED" in (outdir / "validation_issues.tsv").read_text()
