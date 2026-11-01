import gzip
import pytest

from PyamilySeq.order_analysis import collect_contig_replicon_metadata, main, parse_gff
import csv
import json

HDR = "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"


def report_row(name, role="assembled-molecule", molecule="1", location="Chromosome", gb="CM000001.1", ref="NC_000001.1"):
    return f"{name}\t{role}\t{molecule}\t{location}\t{gb}\t=\t{ref}\tPrimary Assembly\t4\t-\n"


def make_cohort(tmp_path, fasta_text, gff_text, report_text=None, compressed=False):
    fd, gd = tmp_path / "fasta", tmp_path / "gff"
    fd.mkdir(); gd.mkdir()
    if compressed:
        with gzip.open(fd / "Ecoli_A.fna.gz", "wt") as f: f.write(fasta_text)
        with gzip.open(gd / "Ecoli_A.gff3.gz", "wt") as f: f.write(gff_text)
    else:
        (fd / "Ecoli_A.fna").write_text(fasta_text)
        (gd / "Ecoli_A.gff3").write_text(gff_text)
    genome = parse_gff(next(gd.iterdir()))
    mapping = None
    if report_text is not None:
        report = tmp_path / "meta" / "assembly_report.txt"
        report.parent.mkdir(); report.write_text(report_text)
        mapping = tmp_path / "mapping.tsv"
        mapping.write_text(f"genome\tassembly_report\nEcoli_A\t{report.relative_to(tmp_path)}\n")
    return collect_contig_replicon_metadata([genome], fd, mapping)


def test_collects_exact_accession_matches_and_shared_fragment_identity(tmp_path):
    fasta = ">NC_000001.1 chromosome\nAAAA\n>CP000001.2 chromosome fragment\nCCCC\n"
    gff = """NC_000001.1\tRef\tregion\t1\t4\t.\t+\t.\tID=r1;gbkey=Src;Is_circular=true
NC_000001.1\tRef\tCDS\t1\t3\t.\t+\t0\tID=g1
CP000001.2\tRef\tregion\t1\t4\t.\t+\t.\tID=r2;gbkey=Src;Is_circular=true
CP000001.2\tRef\tCDS\t1\t3\t.\t+\t0\tID=g2
"""
    report = HDR + report_row("NC_000001.1") + report_row("CP000001.2", gb="CP000001.2", ref="NC_000002.1")
    results = make_cohort(tmp_path, fasta, gff, report)
    assert len(results) == 2
    assert results[0]["replicon_id"] == results[1]["replicon_id"]
    assert all(x["replicon_type"] == "chromosome" and x["replicon_topology"] == "circular" for x in results)
    assert results[0]["matched_refseq_accessions"] == ["NC_000001.1"]
    assert results[1]["matched_genbank_accessions"] == ["CP000001.2"]


def test_fasta_header_gff_and_unknown_records_without_report(tmp_path):
    results = make_cohort(tmp_path, ">chr chromosome linear\nAAAA\n>p plasmid pA\nCCCC\n>bare\nGGGG\n",
                    "chr\tRef\tregion\t1\t4\t.\t+\t.\tID=r;gbkey=Src;Is_circular=false\n"
                    "orphan\tRef\tCDS\t1\t3\t.\t+\t0\tID=o\n")
    by_id = {x["contig"]: x for x in results}
    assert set(by_id) == {"chr", "p", "bare"}
    assert by_id["chr"]["replicon_topology"] == "linear"
    assert by_id["p"]["replicon_type"] == "plasmid"
    assert by_id["bare"]["replicon_type"] == by_id["bare"]["replicon_topology"] == "unknown"
    assert by_id["bare"]["gff_annotation_status"] == "unannotated"


def test_multiple_gff_evidence_and_conflicts_are_retained(tmp_path):
    gff = """ctg\tRef\tregion\t1\t4\t.\t+\t.\tID=r1;gbkey=Src;Is_circular=true
ctg\tRef\tregion\t1\t4\t.\t+\t.\tID=r2;gbkey=plasmid;Is_circular=false
ctg\tRef\tCDS\t1\t3\t.\t+\t0\tID=g
"""
    result = make_cohort(tmp_path, ">ctg chromosome\nAAAA\n", gff,
                   HDR + report_row("ctg"))[0]
    assert result["conflict_status"] != "none"
    assert len(result["metadata_evidence"]["gff_region"]) == 2
    assert result["matched_refseq_accessions"] == []


def test_multiple_exact_accession_tokens_preserve_ambiguity_as_conflict(tmp_path):
    fasta = ">NC_000001.1 CP000009.1\nAAAA\n"
    report = HDR + report_row("chr", "assembled-molecule", "1", "Chromosome", "CM000001.1", "NC_000001.1") + report_row(
        "plasmid", "assembled-molecule", "plasmid p", "Plasmid", "CP000009.1", "NC_100009.1")
    result = make_cohort(tmp_path, fasta, "ctg\tRef\tCDS\t1\t3\t.\t+\t0\tID=g\n", report)[0]
    assert result["matched_refseq_accessions"] == ["NC_000001.1"]
    assert result["matched_genbank_accessions"] == ["CP000009.1"]
    assert result["replicon_type"] == "unknown"
    assert "type" in result["conflict_status"].split(",")


def test_compressed_fasta_gff_and_explicit_report_mapping(tmp_path):
    results = make_cohort(tmp_path, ">CP000001.1 plasmid\nAAAA\n",
                    "CP000001.1\tRef\tregion\t1\t4\t.\t+\t.\tID=r;gbkey=plasmid;Is_circular=true\n"
                    "CP000001.1\tRef\tCDS\t1\t3\t.\t+\t0\tID=g\n",
                    HDR + report_row("ctg", molecule="plasmid pA", location="Plasmid", gb="CP000001.1", ref="NC_100001.1"), True)
    assert len(results) == 1 and results[0]["replicon_type"] == "plasmid"
    assert results[0]["assembly_report_source"].endswith("assembly_report.txt")


def test_missing_report_mapping_file_and_duplicate_genome_assignment_errors(tmp_path):
    fd, gd = tmp_path / "fasta", tmp_path / "gff"
    fd.mkdir(); gd.mkdir()
    (fd / "Ecoli_A.fa").write_text(">c\nAAAA\n")
    (gd / "Ecoli_A.gff3").write_text("c\tRef\tCDS\t1\t3\t.\t+\t0\tID=g\n")
    genome = parse_gff(gd / "Ecoli_A.gff3")
    mapping = tmp_path / "map.tsv"
    mapping.write_text("genome\tassembly_report\nEcoli_A\tmissing.txt\n")
    with pytest.raises(FileNotFoundError): collect_contig_replicon_metadata([genome], fd, mapping)
    report = tmp_path / "ok.txt"
    report.write_text(HDR + report_row("c"))
    mapping.write_text("genome\tassembly_report\nEcoli_A\tok.txt\nEcoli_A\tok.txt\n")
    with pytest.raises(ValueError, match="Duplicate assembly-report mapping"):
        collect_contig_replicon_metadata([genome], fd, mapping)


def test_track_genome_order_writes_resolved_metadata_to_both_output_tables(tmp_path):
    fasta_dir, gff_dir = tmp_path / "fasta", tmp_path / "gff"
    fasta_dir.mkdir(); gff_dir.mkdir()
    records = [
        ("NC_000001.1", "chr1"), ("NC_000002.1", "chr2"), ("NC_000003.1", "chr3"),
        ("CP000001.1", "circ plasmid"), ("CP000002.1", "linear plasmid"),
        ("NW_000001.1", "unlocalized"), ("NZ_000001.1", "unplaced"),
        ("FASTA_ONLY", "fasta only"), ("UNKNOWN", "unknown"), ("CM000099.1", "conflict"),
    ]
    (fasta_dir / "Ecoli_A.fa.gz").write_bytes(gzip.compress("".join(f">{i} {desc}\nACGT\n" for i, desc in records).encode()))
    (fasta_dir / "Ecoli_B.fa.gz").write_bytes(gzip.compress(b">B_unknown\nAAAAA\n"))
    rows = []
    report_records = [
        report_row("NC_000001.1", molecule="1", location="Chromosome", gb="CM000001.1", ref="NC_000001.1"),
        report_row("NC_000002.1", molecule="1", location="Chromosome", gb="CM000002.1", ref="NC_000002.1"),
        report_row("NC_000003.1", molecule="1", location="Chromosome", gb="CM000003.1", ref="NC_000003.1"),
        report_row("CP000001.1", molecule="plasmid pA", location="Plasmid", gb="CP000001.1", ref="NC_100001.1"),
        report_row("CP000002.1", molecule="plasmid pB", location="Plasmid", gb="CP000002.1", ref="NC_100002.1"),
        report_row("NW_000001.1", "unlocalized-scaffold", "1", "Chromosome", "NW_000001.1", "-"),
        report_row("NZ_000001.1", "unplaced-scaffold", "-", "-", "NZ_000001.1", "-"),
        report_row("CM000099.1", molecule="2", location="Chromosome", gb="CM000099.1", ref="NC_999999.1"),
    ]
    report = tmp_path / "assembly_report.txt"
    report.write_text(HDR + "".join(report_records))
    mapping = tmp_path / "report_map.tsv"
    mapping.write_text(f"genome\tassembly_report\nEcoli_A\t{report.name}\n")
    for contig, desc in records:
        if contig == "FASTA_ONLY":
            continue
        attrs = "ID=region;replicon_id=Primary Assembly:chromosome:1;gbkey=Src;Is_circular=true" if contig.startswith("NC_00000") else ""
        if contig == "CP000001.1": attrs = "ID=r;gbkey=plasmid;Is_circular=true;replicon_id=Primary Assembly:plasmid:plasmid pA"
        if contig == "CP000002.1": attrs = "ID=r;gbkey=plasmid;Is_circular=false;replicon_id=Primary Assembly:plasmid:plasmid pB"
        if contig == "CM000099.1": attrs = "ID=r;gbkey=plasmid;Is_circular=true;replicon_id=conflicting_plasmid"
        if attrs:
            rows.append(f"{contig}\tRef\tregion\t1\t4\t.\t+\t.\t{attrs}\n")
        rows.append(f"{contig}\tRef\tCDS\t1\t3\t.\t+\t0\tID=gene_{contig}\n")
    rows.append("GFF_ONLY\tRef\tCDS\t1\t3\t.\t+\t0\tID=gff_only\n")
    (gff_dir / "Ecoli_A.gff3.gz").write_bytes(gzip.compress(("".join(rows)).encode()))
    (gff_dir / "Ecoli_B.gff3.gz").write_bytes(gzip.compress(b"B_unknown\tRef\tCDS\t1\t3\t.\t+\t0\tID=b1\n"))
    outdir = tmp_path / "out"
    rc = main(["--gff-dir", str(gff_dir), "--fasta-dir", str(fasta_dir), "--assembly-report-map", str(mapping),
               "--outdir", str(outdir), "--no-plots"])
    assert rc == 0
    def read(name):
        with (outdir / name).open() as f:
            return list(csv.DictReader(f, delimiter="\t"))
    manifest, endpoints = read("replicon_manifest.tsv"), read("physical_endpoints.tsv")
    assert len(manifest) == len(records) + 1
    assert len(endpoints) == 2 * len(manifest)
    by_key = {(r["genome"], r["contig"]): r for r in manifest}
    assert ("Ecoli_A", "GFF_ONLY") not in by_key
    assert by_key[("Ecoli_A", "FASTA_ONLY")]["has_gff_annotation"] == "False"
    assert by_key[("Ecoli_A", "UNKNOWN")]["replicon_type"] == "unknown"
    assert by_key[("Ecoli_A", "NZ_000001.1")]["replicon_id"] == "unknown"
    assert by_key[("Ecoli_A", "CM000099.1")]["conflict_status"] != "none"
    assert len({by_key[("Ecoli_A", c)]["replicon_id"] for c in ("NC_000001.1", "NC_000002.1", "NC_000003.1")}) == 1
    assert by_key[("Ecoli_A", "CP000001.1")]["replicon_id"] != by_key[("Ecoli_A", "CP000002.1")]["replicon_id"]
    for key, row in by_key.items():
        pair = [e for e in endpoints if (e["genome"], e["contig"]) == key]
        assert len(pair) == 2
        assert {e["endpoint_side"]: int(e["physical_coordinate"]) for e in pair} == {"left": 1, "right": int(row["contig_length"])}
        assert all(e["replicon_id"] == row["replicon_id"] and e["metadata_evidence"] == row["metadata_evidence"] for e in pair)
