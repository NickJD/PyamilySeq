"""Prepare and run cohort-scale pangenome, Sourmash, and gene-order analyses.

This module deliberately keeps one global PyamilySeq pangenome. Sourmash groups
are metadata for stratified summaries; they never partition family discovery.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
from collections import defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from .order_analysis import FASTA_SUFFIXES, GFF_SUFFIXES, _normalise_genome_name, _stem_without_suffix

VERSION = "0.1.0"
MANIFEST_FIELDS = [
    "genome", "fasta_path", "gff_path", "fasta_size_bytes", "n_contigs", "assembly_bp",
    "ambiguous_bp", "sha256", "metadata",
]


def _fasta_stats(path: Path, checksum: bool = False) -> tuple[int, int, int, str]:
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    n_contigs = assembly_bp = ambiguous_bp = 0
    seen: set[str] = set()
    active = False
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(">"):
                header = line[1:].strip().split()
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                contig_id = header[0]
                if contig_id in seen:
                    raise ValueError(f"Duplicate FASTA contig ID {contig_id!r} in {path}")
                seen.add(contig_id)
                n_contigs += 1
                active = True
            else:
                sequence = "".join(line.split())
                if sequence and not active:
                    raise ValueError(f"Sequence data before first FASTA header in {path}:{line_number}")
                assembly_bp += len(sequence)
                ambiguous_bp += sum(base not in "ACGTacgt" for base in sequence)
    if not n_contigs or not assembly_bp:
        raise ValueError(f"FASTA has no sequence records or bases: {path}")
    digest_hex = ""
    if checksum:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        digest_hex = digest.hexdigest()
    return n_contigs, assembly_bp, ambiguous_bp, digest_hex


def _scan_files(directory: Path, suffixes: Iterable[str]) -> list[Path]:
    if not directory.exists():
        raise FileNotFoundError(f"Input directory does not exist: {directory}")
    if not directory.is_dir():
        raise NotADirectoryError(f"Input must be a directory: {directory}")
    accepted = tuple(s.lower() for s in suffixes)
    return sorted((p for p in directory.rglob("*") if p.is_file() and p.name.lower().endswith(accepted)),
                  key=lambda p: (p.name.casefold(), str(p).casefold()))


def _read_metadata(path: Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "genome" not in reader.fieldnames:
            raise ValueError("Metadata TSV must have a 'genome' column")
        records = {}
        for row in reader:
            genome = (row.get("genome") or "").strip()
            if not genome:
                raise ValueError("Metadata TSV contains a blank genome ID")
            key = _normalise_genome_name(genome)
            if key in records:
                raise ValueError(f"Duplicate genome {genome!r} in metadata TSV")
            records[key] = {k: (v or "").strip() for k, v in row.items() if k != "genome"}
        return records


def prepare_cohort(fasta_dir: Path, gff_dir: Path, outdir: Path, metadata_path: Path | None = None,
                   checksums: bool = False) -> Path:
    fasta_paths = _scan_files(fasta_dir, FASTA_SUFFIXES)
    gff_paths = _scan_files(gff_dir, GFF_SUFFIXES)
    if not fasta_paths:
        raise FileNotFoundError(f"No recognised FASTA files found under {fasta_dir}")
    if not gff_paths:
        raise FileNotFoundError(f"No recognised GFF/GFF3 files found under {gff_dir}")

    def unique_by_stem(paths: list[Path], kind: str) -> dict[str, Path]:
        result = {}
        displays = {}
        for path in paths:
            key = _normalise_genome_name(_stem_without_suffix(path, FASTA_SUFFIXES if kind == "FASTA" else GFF_SUFFIXES))
            if key in result:
                raise ValueError(f"Duplicate {kind} genome stem {path.stem!r}: {result[key]} and {path}")
            result[key] = path
            displays[key] = _stem_without_suffix(path, FASTA_SUFFIXES if kind == "FASTA" else GFF_SUFFIXES)
        return result

    fasta_by_key = unique_by_stem(fasta_paths, "FASTA")
    gff_by_key = unique_by_stem(gff_paths, "GFF")
    unmatched_fasta = sorted(set(fasta_by_key) - set(gff_by_key))
    unmatched_gff = sorted(set(gff_by_key) - set(fasta_by_key))
    if unmatched_fasta or unmatched_gff:
        raise ValueError(
            f"FASTA/GFF genome stems must match exactly. Unmatched FASTA: {unmatched_fasta[:10]}; "
            f"unmatched GFF: {unmatched_gff[:10]}"
        )
    metadata = _read_metadata(metadata_path)
    unknown_metadata = sorted(set(metadata) - set(fasta_by_key))
    if unknown_metadata:
        raise ValueError(f"Metadata contains genome IDs absent from FASTA/GFF inputs: {unknown_metadata[:10]}")

    outdir.mkdir(parents=True, exist_ok=True)
    staging = outdir / "pangenome_inputs"
    marker = staging / ".pyamilyseq-cohort-stage"
    if staging.exists() and not marker.exists():
        raise FileExistsError(f"Refusing to modify unowned staging directory: {staging}")
    staging.mkdir(exist_ok=True)
    marker.write_text("Managed by PyamilySeq cohort workflow.\n", encoding="utf-8")

    rows = []
    for key in sorted(fasta_by_key):
        fasta = fasta_by_key[key]
        gff = gff_by_key[key]
        genome = _stem_without_suffix(fasta, FASTA_SUFFIXES)
        n_contigs, assembly_bp, ambiguous_bp, sha256 = _fasta_stats(fasta, checksums)
        staged_gff = staging / f"{genome}.gff3"
        if staged_gff.exists() or staged_gff.is_symlink():
            staged_gff.unlink()
        if gff.name.lower().endswith(".gz"):
            with gzip.open(gff, "rb") as source, staged_gff.open("wb") as target:
                shutil.copyfileobj(source, target)
        else:
            staged_gff.symlink_to(gff.resolve())
        fasta_suffix = ".fna.gz" if fasta.name.lower().endswith(".gz") else ".fna"
        staged_fasta = staging / f"{genome}{fasta_suffix}"
        if staged_fasta.exists() or staged_fasta.is_symlink():
            staged_fasta.unlink()
        staged_fasta.symlink_to(fasta.resolve())
        rows.append({
            "genome": genome, "fasta_path": str(fasta.resolve()), "gff_path": str(gff.resolve()),
            "fasta_size_bytes": fasta.stat().st_size, "n_contigs": n_contigs,
            "assembly_bp": assembly_bp, "ambiguous_bp": ambiguous_bp, "sha256": sha256,
            "metadata": json.dumps(metadata.get(key, {}), sort_keys=True, separators=(",", ":")),
        })

    manifest = outdir / "genome_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    config = {
        "tool": "PyamilySeq cohort analysis", "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "fasta_dir": str(fasta_dir.resolve()), "gff_dir": str(gff_dir.resolve()),
        "metadata_tsv": str(metadata_path.resolve()) if metadata_path else "",
        "genomes": len(rows), "manifest": manifest.name, "pangenome_inputs": staging.name,
    }
    (outdir / "cohort_config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def _load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(MANIFEST_FIELDS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Cohort manifest is missing columns: {', '.join(sorted(missing))}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Cohort manifest is empty: {path}")
    for row in rows:
        for field in ("fasta_path", "gff_path"):
            if not Path(row[field]).is_file():
                raise FileNotFoundError(f"Manifest {field} does not exist for {row['genome']}: {row[field]}")
    return rows


def _run(command: list[str], label: str) -> None:
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(f"Could not run {label}; executable not found: {command[0]}") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"{label} failed with exit code {exc.returncode}: {' '.join(command)}") from exc


def run_sourmash(manifest_path: Path, outdir: Path, sourmash: str = "sourmash", ksize: int = 31,
                 scaled: int = 1000, min_similarity: float | None = None,
                 force: bool = False) -> tuple[Path, Path]:
    if ksize < 1 or scaled < 1:
        raise ValueError("ksize and scaled must be positive integers")
    if min_similarity is not None and not 0 <= min_similarity <= 1:
        raise ValueError("min_similarity must be between 0 and 1")
    rows = _load_manifest(manifest_path)
    similarity_dir = outdir / "similarity"
    similarity_dir.mkdir(parents=True, exist_ok=True)
    signature_path = similarity_dir / f"cohort-k{ksize}-scaled{scaled}.sig.zip"
    matrix_path = similarity_dir / f"similarity-k{ksize}-scaled{scaled}.csv"
    dataset_path = similarity_dir / f"sourmash-inputs-k{ksize}-scaled{scaled}.csv"
    labels_path = similarity_dir / f"similarity-labels-k{ksize}-scaled{scaled}.csv"
    groups_path = similarity_dir / "similarity_groups.tsv"
    edges_path = similarity_dir / "similarity_edges.tsv"
    if not force and any(p.exists() for p in (signature_path, matrix_path, groups_path, edges_path)):
        raise FileExistsError("Sourmash outputs already exist; pass --force to replace them")
    if force:
        for path in (signature_path, matrix_path, labels_path, groups_path, edges_path):
            path.unlink(missing_ok=True)
    with dataset_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "genome_filename", "protein_filename"], lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({"name": row["genome"], "genome_filename": row["fasta_path"], "protein_filename": ""})
    _run([sourmash, "sketch", "fromfile", str(dataset_path), "-p", f"dna,k={ksize},scaled={scaled}", "-o", str(signature_path)], "Sourmash sketching")
    compare = [sourmash, "compare", str(signature_path), "--ksize", str(ksize), "--csv", str(matrix_path),
               "--labels-to", str(labels_path)]
    _run(compare, "Sourmash comparison")
    if min_similarity is not None:
        groups = _write_similarity_groups(matrix_path, labels_path, groups_path, edges_path, min_similarity)
        print(f"Sourmash similarity groups: {len(groups)} groups at threshold {min_similarity:.4f}")
    return matrix_path, signature_path


def _read_matrix_labels(path: Path, expected: int) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        raise ValueError(f"Sourmash label sidecar is empty: {path}")
    header = [value.strip().casefold() for value in rows[0]]
    label_candidates = ("name", "label", "labeltext", "signature", "filename")
    label_index = next((header.index(name) for name in label_candidates if name in header), None)
    if label_index is not None:
        data = rows[1:]
        index_index = next((header.index(name) for name in ("index", "idx") if name in header), None)
        if index_index is not None:
            data.sort(key=lambda row: int(row[index_index]))
        labels = [row[label_index].strip() for row in data]
    elif all(len(row) == 1 for row in rows):
        labels = [row[0].strip() for row in rows]
    elif len(rows[0]) == 2 and rows[0][0].strip().casefold() in {"index", "idx"}:
        labels = [row[1].strip() for row in rows[1:]]
    else:
        raise ValueError(f"Could not identify genome labels in Sourmash sidecar: {path}")
    if len(labels) != expected or any(not label for label in labels):
        raise ValueError(f"Sourmash label sidecar has {len(labels)} labels for a {expected}x{expected} matrix")
    if len(labels) != len(set(labels)):
        raise ValueError("Sourmash label sidecar contains duplicate genome labels")
    return labels


def _read_similarity_matrix(path: Path, labels_path: Path | None = None) -> tuple[list[str], list[list[float]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        raise ValueError(f"Sourmash produced an empty similarity matrix: {path}")
    # Some Sourmash versions emit a labeled matrix; others emit a numeric
    # matrix and put signature names in the --labels-to sidecar.
    header = rows[0]
    labeled = len(header) > 1 and header[0].strip() == "" and all(
        row and row[0].strip() and len(row) == len(header) for row in rows[1:]
    )
    if labeled:
        labels = [x.strip() for x in header[1:]]
        data_rows = rows[1:]
        row_labels = [row[0].strip() for row in data_rows]
        if labels != row_labels:
            raise ValueError("Sourmash similarity matrix must have matching row/column labels")
        values = [row[1:] for row in data_rows]
    elif len(rows) > 1 and len(header) == len(rows) - 1 and all(len(row) == len(header) for row in rows[1:]):
        # Sourmash 4.9 compare --csv writes genome names as the header row,
        # followed by one numeric row per signature (without row-name cells).
        labels = [value.strip() for value in header]
        values = rows[1:]
        if labels_path and labels != _read_matrix_labels(labels_path, len(labels)):
            raise ValueError("Sourmash matrix header labels do not match the --labels-to sidecar order")
    else:
        values = rows
        if any(len(row) != len(values) for row in values):
            raise ValueError("Numeric Sourmash similarity matrix must be square")
        labels = _read_matrix_labels(labels_path, len(values)) if labels_path else []
        if not labels:
            raise ValueError(f"Numeric Sourmash matrix requires its --labels-to sidecar: {path}")
    try:
        matrix = [[float(value) for value in row] for row in values]
    except ValueError as exc:
        raise ValueError(f"Non-numeric similarity in Sourmash matrix {path}") from exc
    if any(len(row) != len(labels) for row in matrix):
        raise ValueError("Sourmash similarity matrix must be square")
    if len(labels) != len(set(labels)):
        raise ValueError("Sourmash similarity matrix contains duplicate genome labels")
    return labels, matrix


def _write_similarity_groups(matrix_path: Path, labels_path: Path, groups_path: Path, edges_path: Path,
                             threshold: float) -> dict[str, str]:
    labels, matrix = _read_similarity_matrix(matrix_path, labels_path)
    adjacency = {label: set() for label in labels}
    edges = []
    for i, left in enumerate(labels):
        for j in range(i + 1, len(labels)):
            similarity = matrix[i][j]
            if similarity >= threshold:
                right = labels[j]
                adjacency[left].add(right)
                adjacency[right].add(left)
                edges.append((left, right, similarity))
    components = []
    unseen = set(labels)
    while unseen:
        seed = min(unseen, key=str.casefold)
        component = set()
        queue = deque([seed])
        unseen.remove(seed)
        while queue:
            current = queue.popleft()
            component.add(current)
            for neighbor in sorted(adjacency[current], key=str.casefold):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    queue.append(neighbor)
        components.append(sorted(component, key=str.casefold))
    components.sort(key=lambda values: values[0].casefold())
    assignments = {genome: f"SM{index:04d}" for index, group in enumerate(components, 1) for genome in group}
    with groups_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome", "similarity_group", "group_size", "threshold", "grouping_method"])
        sizes = {group_id: sum(value == group_id for value in assignments.values()) for group_id in set(assignments.values())}
        for genome in sorted(assignments, key=str.casefold):
            group_id = assignments[genome]
            writer.writerow([genome, group_id, sizes[group_id], threshold, "connected_components_single_linkage"])
    with edges_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome_a", "genome_b", "jaccard_similarity", "threshold"])
        writer.writerows([*edge, threshold] for edge in edges)
    return assignments


def run_pangenome(manifest_path: Path, outdir: Path, executable: str = "PyamilySeq", threads: int = 8,
                  memory_mb: int = 4000, identity: float = 0.9, length_diff: float = 0.8,
                  force: bool = False) -> Path:
    if threads < 1 or memory_mb < 1 or not 0 < identity <= 1 or not 0 < length_diff <= 1:
        raise ValueError("threads and memory must be positive; identity and length difference must be in (0, 1]")
    rows = _load_manifest(manifest_path)
    config_path = manifest_path.parent / "cohort_config.json"
    if not config_path.exists():
        raise FileNotFoundError(f"Cohort config missing; run cohort prepare first: {config_path}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    input_dir = manifest_path.parent / config["pangenome_inputs"]
    output_dir = outdir / "pangenome"
    gpa_path = output_dir / "gene_presence_absence.csv"
    if gpa_path.exists() and not force:
        raise FileExistsError(f"Pangenome output already exists: {gpa_path}; pass --force to replace")
    if not input_dir.is_dir() or not (input_dir / ".pyamilyseq-cohort-stage").exists():
        raise ValueError(f"Managed pangenome inputs missing or unmarked: {input_dir}; rerun cohort prepare")
    command = [executable, "Full", "-clustering_format", "CD-HIT", "-input_type", "separate", "-input_dir", str(input_dir),
               "-name_split_gff", ".gff3", "-output_dir", str(output_dir), "-sequence_type", "AA",
               "-gene_ident", "CDS", "-c", str(identity), "-s", str(length_diff), "-T", str(threads), "-M", str(memory_mb)]
    _run(command, "PyamilySeq global pangenome")
    if not gpa_path.is_file():
        raise RuntimeError(f"PyamilySeq completed without producing expected presence/absence table: {gpa_path}")
    (outdir / "pangenome_run.json").write_text(json.dumps({
        "command": command, "genomes": len(rows), "identity": identity, "length_diff": length_diff,
        "global_family_catalogue": True, "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }, indent=2) + "\n", encoding="utf-8")
    return gpa_path


def _stratify_breakpoints(order_dir: Path, groups_path: Path, output_path: Path, window_genes: int) -> None:
    with groups_path.open("r", encoding="utf-8-sig", newline="") as handle:
        assignments = {row["genome"]: row["similarity_group"] for row in csv.DictReader(handle, delimiter="\t")}
    strata_by_genome: dict[str, list[tuple[str, str]]] = {
        genome: [("similarity_group", group)] for genome, group in assignments.items()
    }
    manifest_path = groups_path.parent.parent / "genome_manifest.tsv"
    if manifest_path.is_file():
        with manifest_path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                try:
                    metadata = json.loads(row.get("metadata", "{}") or "{}")
                except json.JSONDecodeError:
                    metadata = {}
                for field, value in metadata.items():
                    if str(value).strip():
                        strata_by_genome.setdefault(row["genome"], []).append((f"metadata:{field}", str(value)))
    stratum_members = defaultdict(set)
    for genome, strata in strata_by_genome.items():
        for stratum in strata:
            stratum_members[stratum].add(genome)
    input_path = order_dir / "contig_breakpoints_vs_reference.tsv"
    if not input_path.exists():
        raise FileNotFoundError(f"Order-analysis breakpoint table not found: {input_path}")
    endpoints = defaultdict(list)
    with input_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            genome = row["genome"]
            strata = strata_by_genome.get(genome, [])
            if not strata:
                continue
            for side in ("left", "right"):
                if row[f"{side}_boundary_support_eligible"].casefold() != "true":
                    continue
                for stratum in strata:
                    endpoints[(*stratum, row["reference_contig"])].append({
                        "genome": genome, "contig": row["contig"], "side": side,
                        "rank": int(row[f"{side}_boundary_reference_gene_order"]),
                    })
    report_rows = []
    for (stratification, stratum, ref_contig), values in sorted(endpoints.items()):
        values.sort(key=lambda row: (row["rank"], row["genome"].casefold(), row["contig"].casefold()))
        clusters: list[list[dict]] = []
        for value in values:
            if not clusters or value["rank"] - clusters[-1][0]["rank"] > window_genes:
                clusters.append([value])
            else:
                clusters[-1].append(value)
        for cluster in clusters:
            genomes = sorted({row["genome"] for row in cluster}, key=str.casefold)
            ranks = [row["rank"] for row in cluster]
            n_group = len(stratum_members[(stratification, stratum)])
            report_rows.append({
                "stratification": stratification, "stratum": stratum, "reference_contig": ref_contig,
                "reference_gene_order_min": min(ranks), "reference_gene_order_max": max(ranks),
                "n_endpoint_observations": len(cluster), "n_distinct_genomes": len(genomes),
                "n_genomes_in_group": n_group,
                "fraction_group_with_breakpoint": f"{len(genomes) / n_group:.6f}" if n_group else "",
                "genomes": ",".join(genomes), "boundary_sides": ",".join(sorted({r["side"] for r in cluster})),
                "interpretation": "Candidate recurrent endpoint region within this cohort stratum; not causal evidence.",
            })
    fields = ["stratification", "stratum", "reference_contig", "reference_gene_order_min", "reference_gene_order_max",
              "n_endpoint_observations", "n_distinct_genomes", "n_genomes_in_group", "fraction_group_with_breakpoint",
              "genomes", "boundary_sides", "interpretation"]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(report_rows)


def run_order_analysis(manifest_path: Path, outdir: Path, family_map: Path | None = None,
                       order_executable: str = "Track-Genome-Order", reference_genome: str | None = None,
                       no_plots: bool = False, window_genes: int = 10,
                       force: bool = False) -> Path:
    rows = _load_manifest(manifest_path)
    config = json.loads((manifest_path.parent / "cohort_config.json").read_text(encoding="utf-8"))
    family_map = family_map or (outdir / "pangenome" / "gene_presence_absence.csv")
    if not family_map.is_file():
        raise FileNotFoundError(f"Family map not found; run the pangenome stage first: {family_map}")
    order_dir = outdir / "order_analysis"
    expected = order_dir / "run_manifest.json"
    if expected.exists() and not force:
        raise FileExistsError(f"Order-analysis output already exists: {expected}; pass --force to replace")
    # Use the managed staging directory so order analysis sees the exact same
    # genome IDs as the global pangenome, even if source GFF/FASTA stems differ
    # only by case or compressed GFFs needed normalization.
    input_dir = manifest_path.parent / config["pangenome_inputs"]
    command = [order_executable, "--gff-dir", str(input_dir), "--fasta-dir", str(input_dir),
               "--outdir", str(order_dir), "--family-map", str(family_map), "--max-genomes-plotted", "50",
               "--breakpoint-window-genes", str(window_genes)]
    if reference_genome:
        command.extend(["--reference-genome", reference_genome])
    if no_plots:
        command.append("--no-plots")
    _run(command, "Genome order and breakpoint analysis")
    groups_path = outdir / "similarity" / "similarity_groups.tsv"
    if groups_path.exists():
        _stratify_breakpoints(order_dir, groups_path, outdir / "breakpoint_support_by_similarity_group.tsv", window_genes)
    (outdir / "cohort_analysis_summary.json").write_text(json.dumps({
        "genomes": len(rows), "pangenome_family_map": str(family_map), "order_analysis_dir": str(order_dir),
        "similarity_stratification": str(groups_path) if groups_path.exists() else None,
        "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }, indent=2) + "\n", encoding="utf-8")
    return order_dir


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="PyamilySeq-Cohort",
        description="Prepare and run global pangenome, Sourmash similarity, and breakpoint/order analyses.")
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare", help="Validate paired assemblies/annotations and write a cohort manifest.")
    prepare.add_argument("--fasta-dir", type=Path, required=True)
    prepare.add_argument("--gff-dir", type=Path, required=True)
    prepare.add_argument("--metadata", type=Path, help="Optional TSV with one row per genome and a genome column.")
    prepare.add_argument("--outdir", type=Path, required=True)
    prepare.add_argument("--checksums", action="store_true", help="Compute SHA-256 over each FASTA while scanning it.")

    similarity = sub.add_parser("similarity", help="Build Sourmash signatures, matrix, and optional connected-component groups.")
    similarity.add_argument("--manifest", type=Path, required=True)
    similarity.add_argument("--outdir", type=Path, required=True)
    similarity.add_argument("--sourmash", default="sourmash")
    similarity.add_argument("--ksize", type=int, default=31)
    similarity.add_argument("--scaled", type=int, default=1000)
    similarity.add_argument("--min-similarity", type=float, required=True,
                            help="Jaccard threshold for exploratory single-linkage similarity batches (0-1).")
    similarity.add_argument("--force", action="store_true")

    pan = sub.add_parser("pangenome", help="Run one global PyamilySeq full-mode pangenome from the prepared cohort.")
    pan.add_argument("--manifest", type=Path, required=True)
    pan.add_argument("--outdir", type=Path, required=True)
    pan.add_argument("--executable", default="PyamilySeq")
    pan.add_argument("--threads", type=int, default=8)
    pan.add_argument("--memory-mb", type=int, default=4000)
    pan.add_argument("--identity", type=float, default=0.9)
    pan.add_argument("--length-diff", type=float, default=0.8)
    pan.add_argument("--force", action="store_true")

    order = sub.add_parser("order", help="Run whole-cohort order analysis and stratify breakpoints by Sourmash group.")
    order.add_argument("--manifest", type=Path, required=True)
    order.add_argument("--outdir", type=Path, required=True)
    order.add_argument("--family-map", type=Path)
    order.add_argument("--order-executable", default="Track-Genome-Order")
    order.add_argument("--reference-genome")
    order.add_argument("--window-genes", type=int, default=10)
    order.add_argument("--no-plots", action="store_true")
    order.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "prepare":
            manifest = prepare_cohort(args.fasta_dir, args.gff_dir, args.outdir, args.metadata, args.checksums)
            print(f"Prepared {sum(1 for _ in manifest.open()) - 1} genomes: {manifest}")
        elif args.command == "similarity":
            matrix, signatures = run_sourmash(args.manifest, args.outdir, args.sourmash, args.ksize,
                                               args.scaled, args.min_similarity, args.force)
            print(f"Similarity matrix: {matrix}\nSignatures: {signatures}")
        elif args.command == "pangenome":
            gpa = run_pangenome(args.manifest, args.outdir, args.executable, args.threads, args.memory_mb,
                                args.identity, args.length_diff, args.force)
            print(f"Global pangenome complete: {gpa}")
        else:
            order_dir = run_order_analysis(args.manifest, args.outdir, args.family_map, args.order_executable,
                                           args.reference_genome, args.no_plots,
                                           args.window_genes, args.force)
            print(f"Cohort order analysis complete: {order_dir}")
        return 0
    except (FileNotFoundError, NotADirectoryError, FileExistsError, ValueError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
