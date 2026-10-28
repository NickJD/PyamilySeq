import csv
import gzip
from pathlib import Path

import pytest

from PyamilySeq.cohort_analysis import (
    _stratify_breakpoints,
    _write_similarity_groups,
    prepare_cohort,
    run_pangenome,
    run_sourmash,
)


def _write_inputs(root: Path):
    fasta_dir = root / "fasta"
    gff_dir = root / "gff"
    fasta_dir.mkdir()
    gff_dir.mkdir()
    for genome, seq in (("sampleA", "ACGTACGT"), ("sampleB", "ACGTAAAA")):
        (fasta_dir / f"{genome}.fna").write_text(f">contig1\n{seq}\n")
        (gff_dir / f"{genome}.gff3").write_text(
            "ctg\tRef\tCDS\t1\t4\t.\t+\t0\tID=g1\n"
        )
    return fasta_dir, gff_dir


def test_prepare_builds_manifest_and_staged_pangenome_inputs(tmp_path):
    fasta_dir, gff_dir = _write_inputs(tmp_path)
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("genome\tphylogroup\nsampleA\tA\nsampleB\tB\n")
    outdir = tmp_path / "cohort"
    manifest = prepare_cohort(fasta_dir, gff_dir, outdir, metadata, checksums=True)
    with manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["genome"] for row in rows] == ["sampleA", "sampleB"]
    assert rows[0]["assembly_bp"] == "8"
    assert len(rows[0]["sha256"]) == 64
    assert rows[0]["metadata"] == '{"phylogroup":"A"}'
    stage = outdir / "pangenome_inputs"
    assert (stage / ".pyamilyseq-cohort-stage").exists()
    assert (stage / "sampleA.gff3").is_symlink()
    assert (stage / "sampleA.fna").is_symlink()


def test_prepare_decompresses_gff_and_rejects_unpaired_genomes(tmp_path):
    fasta_dir, gff_dir = _write_inputs(tmp_path)
    text = (gff_dir / "sampleA.gff3").read_text()
    (gff_dir / "sampleA.gff3").unlink()
    with gzip.open(gff_dir / "sampleA.gff3.gz", "wt") as handle:
        handle.write(text)
    manifest = prepare_cohort(fasta_dir, gff_dir, tmp_path / "cohort")
    assert "CDS" in (tmp_path / "cohort" / "pangenome_inputs" / "sampleA.gff3").read_text()

    (fasta_dir / "unpaired.fna").write_text(">x\nACGT\n")
    with pytest.raises(ValueError, match="must match exactly"):
        prepare_cohort(fasta_dir, gff_dir, tmp_path / "other")


def test_similarity_groups_are_connected_components_and_keep_threshold_edges(tmp_path):
    matrix = tmp_path / "matrix.csv"
    matrix.write_text(",a,b,c,d\na,1,0.8,0.1,0.1\nb,0.8,1,0.7,0.1\nc,0.1,0.7,1,0.2\nd,0.1,0.1,0.2,1\n")
    groups, edges = tmp_path / "groups.tsv", tmp_path / "edges.tsv"
    assignments = _write_similarity_groups(matrix, tmp_path / "labels.csv", groups, edges, 0.65)
    assert assignments["a"] == assignments["b"] == assignments["c"]
    assert assignments["d"] != assignments["a"]
    assert "connected_components_single_linkage" in groups.read_text()
    assert len(edges.read_text().splitlines()) == 3


def test_similarity_matrix_reads_sourmash_header_plus_numeric_rows(tmp_path):
    from PyamilySeq.cohort_analysis import _read_similarity_matrix

    matrix = tmp_path / "matrix.csv"
    matrix.write_text("a,b,c\n1,0.8,0.1\n0.8,1,0.7\n0.1,0.7,1\n")
    labels = tmp_path / "labels.csv"
    labels.write_text("sort_order,md5,label,name,filename,signature_file\n"
                      "1,x,a,a,a,a\n2,y,b,b,b,b\n3,z,c,c,c,c\n")
    parsed_labels, values = _read_similarity_matrix(matrix, labels)
    assert parsed_labels == ["a", "b", "c"]
    assert values[1][2] == 0.7


def test_sourmash_stage_uses_manifest_and_writes_matrix_groups(tmp_path):
    fasta_dir, gff_dir = _write_inputs(tmp_path)
    manifest = prepare_cohort(fasta_dir, gff_dir, tmp_path / "cohort")
    fake = tmp_path / "sourmash"
    fake.write_text("""#!/usr/bin/env python3
import csv, pathlib, sys
args = sys.argv[1:]
if args[:2] == ['sketch', 'fromfile']:
    data = list(csv.DictReader(open(args[2])))
    out = pathlib.Path(args[args.index('-o') + 1])
    out.write_text('fake signatures')
elif args[0] == 'compare':
    out = pathlib.Path(args[args.index('--csv') + 1])
    labels = ['sampleA', 'sampleB']
    with out.open('w') as f: f.write('1,0.9\\n0.9,1\\n')
    pathlib.Path(args[args.index('--labels-to') + 1]).write_text('index,name\\n0,sampleA\\n1,sampleB\\n')
else: raise SystemExit(3)
""")
    fake.chmod(0o755)
    matrix, signatures = run_sourmash(manifest, tmp_path / "cohort", str(fake), min_similarity=0.8)
    assert matrix.is_file() and signatures.is_file()
    assert (tmp_path / "cohort" / "similarity" / "similarity_groups.tsv").is_file()


def test_global_pangenome_invokes_supported_full_mode_and_checks_output(tmp_path):
    fasta_dir, gff_dir = _write_inputs(tmp_path)
    manifest = prepare_cohort(fasta_dir, gff_dir, tmp_path / "cohort")
    fake = tmp_path / "pyamilyseq"
    fake.write_text("""#!/usr/bin/env python3
import pathlib, sys
args = sys.argv[1:]
assert args[:4] == ['Full', '-clustering_format', 'CD-HIT', '-input_type']
outdir = pathlib.Path(args[args.index('-output_dir') + 1])
outdir.mkdir(parents=True, exist_ok=True)
(outdir / 'gene_presence_absence.csv').write_text('Gene,Non-unique Gene name,Annotation,sampleA,sampleB\\n')
""")
    fake.chmod(0o755)
    gpa = run_pangenome(manifest, tmp_path / "cohort", str(fake))
    assert gpa.is_file()
    assert '"global_family_catalogue": true' in (tmp_path / "cohort" / "pangenome_run.json").read_text()


def test_stratify_breakpoints_reports_group_denominators(tmp_path):
    order_dir = tmp_path / "order"
    order_dir.mkdir()
    (order_dir / "contig_breakpoints_vs_reference.tsv").write_text(
        "genome\tcontig\treference_contig\tleft_boundary_reference_gene_order\tleft_boundary_support_eligible\t"
        "right_boundary_reference_gene_order\tright_boundary_support_eligible\n"
        "a\tca\tr1\t10\tTrue\t90\tFalse\n"
        "b\tcb\tr1\t12\tTrue\t91\tFalse\n"
        "c\tcc\tr1\t80\tFalse\t92\tFalse\n"
    )
    similarity_dir = tmp_path / "similarity"
    similarity_dir.mkdir()
    groups = similarity_dir / "similarity_groups.tsv"
    groups.write_text("genome\tsimilarity_group\na\tSM0001\nb\tSM0001\nc\tSM0001\n")
    (tmp_path / "genome_manifest.tsv").write_text(
        "genome\tmetadata\na\t{\"phylogroup\":\"A\"}\nb\t{\"phylogroup\":\"A\"}\nc\t{\"phylogroup\":\"B\"}\n"
    )
    output = tmp_path / "stratified.tsv"
    _stratify_breakpoints(order_dir, groups, output, window_genes=5)
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    similarity = next(row for row in rows if row["stratification"] == "similarity_group")
    assert similarity["n_distinct_genomes"] == "2"
    assert similarity["n_genomes_in_group"] == "3"
    assert similarity["fraction_group_with_breakpoint"] == "0.666667"
    phylogroup_a = next(row for row in rows if row["stratification"] == "metadata:phylogroup" and row["stratum"] == "A")
    assert phylogroup_a["n_distinct_genomes"] == "2"
    assert phylogroup_a["n_genomes_in_group"] == "2"
