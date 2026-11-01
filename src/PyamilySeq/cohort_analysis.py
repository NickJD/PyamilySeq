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
import re
import shutil
import shlex
import statistics
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
    "ambiguous_bp", "n50_bp", "l50_contigs", "longest_contig_bp", "median_contig_bp", "sha256", "metadata",
]


def _replicon_header_metadata(header: str) -> tuple[str, str, bool]:
    """Return replicon class, topology and completeness from explicit FASTA wording.

    Unknown is retained unless the submitter header explicitly says chromosome/
    plasmid and complete. These labels describe the record, not read-level proof.
    """
    text = header.casefold()
    kind = "plasmid" if re.search(r"\bplasmid\b", text) else (
        "chromosome" if re.search(r"\bchromosome\b", text) else "unknown")
    complete = bool(re.search(r"\bcomplete (?:genome|sequence)\b", text))
    topology = "circular" if re.search(r"\bcircular\b", text) else (
        "complete_unknown_topology" if complete else "unknown")
    return kind, topology, complete


def _fasta_record_metadata(path: Path) -> list[dict[str, str]]:
    rows = []
    with _open_fasta(path) as handle:
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip()
                seqid = header.split(maxsplit=1)[0]
                kind, topology, complete = _replicon_header_metadata(header)
                rows.append({"contig": seqid, "description": header,
                             "replicon_type": kind, "topology": topology,
                             "complete_record": str(complete).lower()})
    return rows


def _reference_score(row: dict[str, str]) -> tuple[int, int, int, str]:
    """Prefer one explicit complete chromosome, then assembly continuity."""
    records = _fasta_record_metadata(Path(row["fasta_path"]))
    complete_chromosomes = sum(
        record["replicon_type"] == "chromosome" and record["complete_record"] == "true"
        for record in records
    )
    # A single chromosome record is a useful coordinate reference; plasmid records
    # do not disqualify it. Topology remains unknown unless explicitly annotated.
    preferred = int(complete_chromosomes == 1)
    return (preferred, int(row.get("n50_bp") or 0), -int(row.get("n_contigs") or 0), row["genome"].casefold())


def _n50_l50(lengths: list[int]) -> tuple[int, int]:
    ordered = sorted((length for length in lengths if length > 0), reverse=True)
    target = sum(ordered) / 2
    cumulative = 0
    for index, length in enumerate(ordered, 1):
        cumulative += length
        if cumulative >= target:
            return length, index
    return 0, 0


def _fasta_stats(path: Path, checksum: bool = False) -> tuple[int, int, int, str, int, int, int, float]:
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    n_contigs = assembly_bp = ambiguous_bp = 0
    seen: set[str] = set()
    contig_lengths: list[int] = []
    active = False
    length = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(">"):
                header = line[1:].strip().split()
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                contig_id = header[0]
                if contig_id in seen:
                    raise ValueError(f"Duplicate FASTA contig ID {contig_id!r} in {path}")
                if active:
                    contig_lengths.append(length)
                seen.add(contig_id)
                n_contigs += 1
                active = True
                length = 0
            else:
                sequence = "".join(line.split())
                if sequence and not active:
                    raise ValueError(f"Sequence data before first FASTA header in {path}:{line_number}")
                assembly_bp += len(sequence)
                ambiguous_bp += sum(base not in "ACGTacgt" for base in sequence)
                length += len(sequence)
    if not n_contigs or not assembly_bp:
        raise ValueError(f"FASTA has no sequence records or bases: {path}")
    if active:
        contig_lengths.append(length)
    n50, l50 = _n50_l50(contig_lengths)
    digest_hex = ""
    if checksum:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        digest_hex = digest.hexdigest()
    return (n_contigs, assembly_bp, ambiguous_bp, digest_hex, n50, l50,
            max(contig_lengths, default=0), statistics.median(contig_lengths) if contig_lengths else 0)


def _scan_files(directory: Path, suffixes: Iterable[str]) -> list[Path]:
    if not directory.exists():
        raise FileNotFoundError(f"Input directory does not exist: {directory}")
    if not directory.is_dir():
        raise NotADirectoryError(f"Input must be a directory: {directory}")
    accepted = tuple(s.lower() for s in suffixes)
    return sorted((p for p in directory.rglob("*") if p.is_file() and p.name.lower().endswith(accepted)),
                  key=lambda p: (p.name.casefold(), str(p).casefold()))


def split_combined_fasta(input_fasta: Path, output_dir: Path,
                         genome_header_regex: str | None = None) -> tuple[Path, Path]:
    """Split a concatenated multi-genome FASTA into per-genome FASTAs.

    With no custom regex, headers containing a BV-BRC/RefSeq-style bracketed
    label such as ``[Escherichia coli strain X | 562.160571]`` are grouped by
    the text before the pipe. Custom regexes must contain a named ``genome``
    capture group. Sequence record IDs are preserved verbatim.
    """
    if not input_fasta.is_file():
        raise FileNotFoundError(f"Combined assembly FASTA does not exist: {input_fasta}")
    custom = re.compile(genome_header_regex) if genome_header_regex else None
    default = re.compile(r"\[([^\]|]+?)\s*\|")
    output_dir.mkdir(parents=True, exist_ok=True)
    handles: dict[str, object] = {}
    names: dict[str, str] = {}
    paths_by_key: dict[str, Path] = {}
    ids_by_genome: dict[str, set[str]] = defaultdict(set)
    manifest = output_dir / "combined_fasta_genomes.tsv"
    previously_managed: set[Path] = set()
    if manifest.is_file():
        with manifest.open("r", encoding="utf-8-sig", newline="") as old_handle:
            for old_row in csv.DictReader(old_handle, delimiter="\t"):
                old_name = Path(old_row.get("fasta", "")).name
                if old_name:
                    previously_managed.add(output_dir / old_name)
    active_handle = None
    active_genome = ""
    active_id = ""
    records = 0
    try:
        with _open_fasta(input_fasta) as source:
            for line_number, raw in enumerate(source, 1):
                if raw.startswith(">"):
                    header = raw[1:].strip()
                    record_id = header.split(maxsplit=1)[0] if header else ""
                    if not record_id:
                        raise ValueError(f"Empty FASTA header at {input_fasta}:{line_number}")
                    if custom:
                        match = custom.search(header)
                        genome = match.group("genome").strip() if match else ""
                    else:
                        match = default.search(header)
                        genome = match.group(1).strip() if match else ""
                    if not genome:
                        raise ValueError(
                            f"Could not identify a genome from FASTA header at {input_fasta}:{line_number}: {header!r}. "
                            "Use --genome-header-regex with a named (?P<genome>...) group."
                        )
                    key = _normalise_genome_name(genome)
                    if record_id in ids_by_genome[key]:
                        raise ValueError(f"Duplicate sequence ID {record_id!r} within genome {genome!r}")
                    ids_by_genome[key].add(record_id)
                    if key not in names:
                        names[key] = genome
                        slug = re.sub(r"[^A-Za-z0-9._-]+", "_", genome).strip("._-") or f"genome_{len(names)}"
                        if any(path.stem.casefold() == slug.casefold() for path in paths_by_key.values()):
                            slug = f"{slug}_{len(names)}"
                        fasta_path = output_dir / f"{slug}.fna"
                        paths_by_key[key] = fasta_path
                        handle = fasta_path.open("w", encoding="utf-8")
                        handles[key] = handle
                    active_genome, active_id = key, record_id
                    active_handle = handles[key]
                    active_handle.write(raw if raw.endswith("\n") else raw + "\n")
                    records += 1
                else:
                    if active_handle is None:
                        if raw.strip():
                            raise ValueError(f"Sequence data before first FASTA header at {input_fasta}:{line_number}")
                        continue
                    active_handle.write(raw if raw.endswith("\n") else raw + "\n")
        if not records:
            raise ValueError(f"No FASTA records found in {input_fasta}")
    finally:
        for handle in handles.values():
            handle.close()
    with manifest.open("w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome", "genome_key", "fasta", "n_sequences"])
        for key in sorted(names):
            fasta_path = paths_by_key[key]
            writer.writerow([names[key], key, fasta_path.name, len(ids_by_genome[key])])
    current_paths = set(paths_by_key.values())
    for stale_path in previously_managed - current_paths:
        if stale_path.parent == output_dir:
            stale_path.unlink(missing_ok=True)
    return output_dir, manifest


def _open_fasta(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.name.lower().endswith(".gz") else path.open("r", encoding="utf-8", errors="replace")


def _iter_fasta_records(path: Path):
    header = None
    sequence = []
    with _open_fasta(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:].strip().split(maxsplit=1)[0] if line[1:].strip() else ""
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                sequence = []
            else:
                if header is None:
                    raise ValueError(f"Sequence data before first FASTA header in {path}:{line_number}")
                sequence.append("".join(line.split()))
    if header is not None:
        yield header, "".join(sequence)


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


def prepare_cohort(fasta_dir: Path, gff_dir: Path | None, outdir: Path, metadata_path: Path | None = None,
                   checksums: bool = False) -> Path:
    fasta_paths = _scan_files(fasta_dir, FASTA_SUFFIXES)
    gff_paths = _scan_files(gff_dir, GFF_SUFFIXES) if gff_dir is not None else []
    if not fasta_paths:
        raise FileNotFoundError(f"No recognised FASTA files found under {fasta_dir}")
    if gff_dir is not None and not gff_paths:
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
    gff_by_key = unique_by_stem(gff_paths, "GFF") if gff_paths else {}
    unmatched_fasta = sorted(set(fasta_by_key) - set(gff_by_key))
    unmatched_gff = sorted(set(gff_by_key) - set(fasta_by_key))
    if gff_paths and (unmatched_fasta or unmatched_gff):
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

    expected_staged_names = set()
    for key in sorted(fasta_by_key):
        genome = _stem_without_suffix(fasta_by_key[key], FASTA_SUFFIXES)
        fasta_suffix = ".fna.gz" if fasta_by_key[key].name.lower().endswith(".gz") else ".fna"
        expected_staged_names.add(f"{genome}{fasta_suffix}")
        if key in gff_by_key:
            expected_staged_names.add(f"{genome}.gff3")
    for staged_path in staging.iterdir():
        if staged_path.name != marker.name and staged_path.name not in expected_staged_names:
            if staged_path.is_file() or staged_path.is_symlink():
                staged_path.unlink()
            else:
                raise ValueError(f"Unexpected directory in managed pangenome staging area: {staged_path}")

    rows = []
    assembly_stats = []
    for key in sorted(fasta_by_key):
        fasta = fasta_by_key[key]
        gff = gff_by_key.get(key)
        genome = _stem_without_suffix(fasta, FASTA_SUFFIXES)
        (n_contigs, assembly_bp, ambiguous_bp, sha256, n50, l50,
         longest, median) = _fasta_stats(fasta, checksums)
        staged_gff_path = ""
        if gff is not None:
            staged_gff = staging / f"{genome}.gff3"
            if staged_gff.exists() or staged_gff.is_symlink():
                staged_gff.unlink()
            if gff.name.lower().endswith(".gz"):
                with gzip.open(gff, "rb") as source, staged_gff.open("wb") as target:
                    shutil.copyfileobj(source, target)
            else:
                staged_gff.symlink_to(gff.resolve())
            staged_gff_path = str(gff.resolve())
        fasta_suffix = ".fna.gz" if fasta.name.lower().endswith(".gz") else ".fna"
        staged_fasta = staging / f"{genome}{fasta_suffix}"
        if staged_fasta.exists() or staged_fasta.is_symlink():
            staged_fasta.unlink()
        staged_fasta.symlink_to(fasta.resolve())
        rows.append({
            "genome": genome, "fasta_path": str(fasta.resolve()), "gff_path": staged_gff_path,
            "fasta_size_bytes": fasta.stat().st_size, "n_contigs": n_contigs,
            "assembly_bp": assembly_bp, "ambiguous_bp": ambiguous_bp, "n50_bp": n50,
            "l50_contigs": l50, "longest_contig_bp": longest, "median_contig_bp": median, "sha256": sha256,
            "metadata": json.dumps(metadata.get(key, {}), sort_keys=True, separators=(",", ":")),
        })
        assembly_stats.append({"genome": genome, "n_contigs": n_contigs, "assembly_bp": assembly_bp,
                               "ambiguous_bp": ambiguous_bp, "n50_bp": n50, "l50_contigs": l50,
                               "longest_contig_bp": longest, "median_contig_bp": median})

    manifest = outdir / "genome_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with (outdir / "assembly_metrics.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(assembly_stats[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(assembly_stats)
    config = {
        "tool": "PyamilySeq cohort analysis", "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "fasta_dir": str(fasta_dir.resolve()), "gff_dir": str(gff_dir.resolve()) if gff_dir else "",
        "metadata_tsv": str(metadata_path.resolve()) if metadata_path else "",
        "genomes": len(rows), "manifest": manifest.name, "pangenome_inputs": staging.name,
    }
    (outdir / "cohort_config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def _flatten_managed_gffs(gff_dir: Path, prefixes: Iterable[str]) -> Path:
    """Refresh the workflow-owned flattened GFF directory without stale links."""
    flattened = gff_dir / "gff3"
    flattened.mkdir(exist_ok=True)
    marker = flattened / ".pyamilyseq-gff-links"
    marker.write_text("Managed GFF links created by PyamilySeq cohort workflow.\n", encoding="utf-8")
    expected = {f"{prefix}.gff3" for prefix in prefixes}
    root = gff_dir.resolve()
    for link in flattened.glob("*.gff3"):
        if link.is_symlink():
            target = link.resolve(strict=False)
            try:
                target.relative_to(root)
            except ValueError:
                continue
            if link.name not in expected:
                link.unlink()
    for prefix in sorted(prefixes, key=str.casefold):
        source = gff_dir / prefix / f"{prefix}.gff3"
        if not source.is_file():
            raise RuntimeError(f"Expected managed GFF3 output is missing: {source}")
        target = flattened / source.name
        if target.exists() or target.is_symlink():
            target.unlink()
        target.symlink_to(source.resolve())
    return flattened


def annotate_with_bakta(fasta_dir: Path, gff_dir: Path, executable: str = "bakta",
                        database: Path | None = None, threads: int = 1,
                        force: bool = False) -> Path:
    """Create standardized per-genome GFF3 annotations with a local Bakta install."""
    if threads < 1:
        raise ValueError("Bakta threads must be positive")
    fasta_paths = _scan_files(fasta_dir, FASTA_SUFFIXES)
    if not fasta_paths:
        raise FileNotFoundError(f"No recognised FASTA files found under {fasta_dir}")
    gff_dir.mkdir(parents=True, exist_ok=True)
    for fasta in fasta_paths:
        prefix = _stem_without_suffix(fasta, FASTA_SUFFIXES)
        output = gff_dir / prefix
        expected_gff = output / f"{prefix}.gff3"
        if expected_gff.is_file() and not force:
            continue
        command = [executable, "--output", str(output), "--prefix", prefix, "--threads", str(threads)]
        if database:
            command.extend(["--db", str(database)])
        if force:
            command.append("--force")
        command.append(str(fasta))
        _run(command, f"Bakta annotation for {prefix}")
        if not expected_gff.is_file():
            raise RuntimeError(f"Bakta completed without the expected GFF3 output: {expected_gff}")
    return _flatten_managed_gffs(gff_dir, [
        _stem_without_suffix(fasta, FASTA_SUFFIXES) for fasta in _scan_files(fasta_dir, FASTA_SUFFIXES)
    ])


def annotate_with_pyrodigal(fasta_dir: Path, gff_dir: Path, force: bool = False) -> Path:
    """Predict CDS features per genome using a genome-trained Pyrodigal model."""
    try:
        import pyrodigal
    except ImportError as exc:
        raise RuntimeError("Pyrodigal annotation requires `pip install PyamilySeq[annotation]` or pyrodigal in this environment") from exc
    fasta_paths = _scan_files(fasta_dir, FASTA_SUFFIXES)
    if not fasta_paths:
        raise FileNotFoundError(f"No recognised FASTA files found under {fasta_dir}")
    gff_dir.mkdir(parents=True, exist_ok=True)
    for fasta in fasta_paths:
        prefix = _stem_without_suffix(fasta, FASTA_SUFFIXES)
        output = gff_dir / prefix
        output.mkdir(parents=True, exist_ok=True)
        expected_gff = output / f"{prefix}.gff3"
        if expected_gff.is_file() and not force:
            continue
        records = list(_iter_fasta_records(fasta))
        if not records or any(not sequence for _, sequence in records):
            raise ValueError(f"FASTA contains no records or an empty sequence: {fasta}")
        training_sequence = max(records, key=lambda item: len(item[1]))[1]
        finder = pyrodigal.GeneFinder()
        training_info = finder.train(training_sequence)
        finder = pyrodigal.GeneFinder(training_info=training_info, closed=False)
        with expected_gff.open("w", encoding="utf-8") as handle:
            handle.write("##gff-version 3\n")
            for sequence_id, sequence in records:
                genes = finder.find_genes(sequence)
                genes.write_gff(handle, sequence_id=sequence_id, header=False)
    return _flatten_managed_gffs(gff_dir, [
        _stem_without_suffix(fasta, FASTA_SUFFIXES) for fasta in fasta_paths
    ])


def prepare_cohort_from_combined_fasta(combined_fasta: Path, gff_dir: Path | None, outdir: Path,
                                       metadata_path: Path | None = None, checksums: bool = False,
                                       genome_header_regex: str | None = None,
                                       annotate: bool = False, annotate_pyrodigal_enabled: bool = False,
                                       bakta: str = "bakta",
                                       bakta_database: Path | None = None, bakta_threads: int = 1,
                                       force: bool = False) -> Path:
    split_dir = outdir / "split_genomes"
    split_combined_fasta(combined_fasta, split_dir, genome_header_regex)
    if annotate:
        gff_dir = annotate_with_bakta(split_dir, outdir / "bakta_annotations", bakta,
                                      bakta_database, bakta_threads, force)
    elif annotate_pyrodigal_enabled:
        gff_dir = annotate_with_pyrodigal(split_dir, outdir / "pyrodigal_annotations", force)
    return prepare_cohort(split_dir, gff_dir, outdir, metadata_path, checksums)


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
            if field == "gff_path" and not row[field]:
                continue
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
                 force: bool = False, grouping_method: str = "greedy_complete_linkage") -> tuple[Path, Path]:
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
        groups = _write_similarity_groups(matrix_path, labels_path, groups_path, edges_path,
                                           min_similarity, grouping_method)
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
                             threshold: float, grouping_method: str = "greedy_complete_linkage") -> dict[str, str]:
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
    if grouping_method == "single_linkage":
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
    elif grouping_method == "greedy_complete_linkage":
        label_index = {label: index for index, label in enumerate(labels)}
        for label in sorted(labels, key=str.casefold):
            eligible = []
            for index, cluster in enumerate(components):
                minimum = min(matrix[label_index[label]][label_index[member]] for member in cluster)
                if minimum >= threshold:
                    eligible.append((minimum, index))
            if eligible:
                _minimum, index = max(eligible, key=lambda item: (item[0], -item[1]))
                components[index].append(label)
            else:
                components.append([label])
    else:
        raise ValueError(f"Unknown similarity grouping method: {grouping_method}")
    components.sort(key=lambda values: values[0].casefold())
    assignments = {genome: f"SM{index:04d}" for index, group in enumerate(components, 1) for genome in group}
    with groups_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome", "similarity_group", "group_size", "threshold", "grouping_method"])
        sizes = {group_id: sum(value == group_id for value in assignments.values()) for group_id in set(assignments.values())}
        for genome in sorted(assignments, key=str.casefold):
            group_id = assignments[genome]
            writer.writerow([genome, group_id, sizes[group_id], threshold, grouping_method])
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
    if any(not row.get("gff_path") for row in rows):
        raise ValueError("Pangenome analysis requires a GFF annotation for every genome; prepare with --gff-dir or --annotate-with-bakta")
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
    mapped_spans = defaultdict(list)
    with input_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            genome = row["genome"]
            strata = strata_by_genome.get(genome, [])
            if not strata:
                continue
            if row.get("reference_order_min_anchor") and row.get("reference_order_max_anchor"):
                for stratum in strata:
                    mapped_spans[(*stratum, row["reference_contig"])].append(
                        (genome, int(row["reference_order_min_anchor"]), int(row["reference_order_max_anchor"]))
                    )
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
            region_start = max(1, min(ranks) - window_genes)
            region_end = max(ranks) + window_genes
            callable_genomes = set(genomes)
            for genome, span_start, span_end in mapped_spans[(stratification, stratum, ref_contig)]:
                if span_start <= region_start and span_end >= region_end:
                    callable_genomes.add(genome)
            n_group = len(stratum_members[(stratification, stratum)])
            n_callable = len(callable_genomes)
            report_rows.append({
                "stratification": stratification, "stratum": stratum, "reference_contig": ref_contig,
                "reference_gene_order_min": min(ranks), "reference_gene_order_max": max(ranks),
                "n_endpoint_observations": len(cluster), "n_distinct_genomes": len(genomes),
                "n_genomes_in_group": n_group,
                "n_genomes_callable_for_region": n_callable,
                "fraction_callable_with_contig_end": f"{len(genomes) / n_callable:.6f}" if n_callable else "",
                "genomes": ",".join(genomes), "boundary_sides": ",".join(sorted({r["side"] for r in cluster})),
                "interpretation": "Candidate recurrent endpoint region within this cohort stratum; not causal evidence.",
            })
    fields = ["stratification", "stratum", "reference_contig", "reference_gene_order_min", "reference_gene_order_max",
              "n_endpoint_observations", "n_distinct_genomes", "n_genomes_in_group",
              "n_genomes_callable_for_region", "fraction_callable_with_contig_end",
              "genomes", "boundary_sides", "interpretation"]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(report_rows)


def run_order_analysis(manifest_path: Path, outdir: Path, family_map: Path | None = None,
                       order_executable: str = "Track-Genome-Order", reference_genome: str | None = None,
                       no_plots: bool = False, window_genes: int = 10,
                       force: bool = False, per_similarity_group: bool = False) -> Path:
    rows = _load_manifest(manifest_path)
    if any(not row.get("gff_path") for row in rows):
        raise ValueError("Order analysis requires a GFF annotation for every genome; prepare with --gff-dir or --annotate-with-bakta")
    config = json.loads((manifest_path.parent / "cohort_config.json").read_text(encoding="utf-8"))
    family_map = family_map or (outdir / "pangenome" / "gene_presence_absence.csv")
    if not family_map.is_file():
        raise FileNotFoundError(f"Family map not found; run the pangenome stage first: {family_map}")
    order_dir = outdir / "order_analysis"
    expected = order_dir / "run_manifest.json"
    groups_path = outdir / "similarity" / "similarity_groups.tsv"
    if per_similarity_group:
        if not groups_path.is_file():
            raise FileNotFoundError("Per-group order analysis requires similarity groups; run cohort similarity with --min-similarity first")
        group_members: dict[str, list[str]] = defaultdict(list)
        with groups_path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                group_members[row["similarity_group"]].append(row["genome"])
        group_manifest = order_dir / "group_runs.tsv"
        if group_manifest.exists() and not force:
            raise FileExistsError(f"Per-group order outputs already exist: {order_dir}; pass --force to replace")
        if force and order_dir.exists():
            # Only replace the specific directories this mode owns. The output
            # root may also contain unrelated order-analysis results.
            if not group_manifest.is_file():
                raise FileExistsError(
                    f"Cannot safely replace per-group outputs without their manifest: {order_dir}"
                )
            for managed in (order_dir / "group_inputs", order_dir / "groups"):
                if managed.is_dir():
                    shutil.rmtree(managed)
            group_manifest.unlink()
        order_dir.mkdir(parents=True, exist_ok=True)
        input_dir = manifest_path.parent / config["pangenome_inputs"]
        rows_by_genome = {row["genome"]: row for row in rows}
        group_run_rows = []
        for group_id, members in sorted(group_members.items()):
            missing = sorted(set(members) - set(rows_by_genome))
            if missing:
                raise ValueError(f"Similarity group {group_id} names genomes absent from the manifest: {missing[:5]}")
            if reference_genome:
                if reference_genome not in members:
                    raise ValueError(f"Requested reference genome {reference_genome!r} is not in similarity group {group_id}")
                selected_reference = reference_genome
            else:
                selected_reference = max(members, key=lambda genome: _reference_score(rows_by_genome[genome]))
            group_input = order_dir / "group_inputs" / group_id
            group_gff = group_input / "gff"
            group_fasta = group_input / "fasta"
            group_gff.mkdir(parents=True, exist_ok=True)
            group_fasta.mkdir(parents=True, exist_ok=True)
            for genome in members:
                row = rows_by_genome[genome]
                source_gff = input_dir / f"{genome}.gff3"
                source_fasta = Path(row["fasta_path"])
                (group_gff / source_gff.name).symlink_to(source_gff.resolve())
                (group_fasta / source_fasta.name).symlink_to(source_fasta.resolve())
            group_out = order_dir / "groups" / group_id
            command = [order_executable, "--gff-dir", str(group_gff), "--fasta-dir", str(group_fasta),
                       "--outdir", str(group_out), "--family-map", str(family_map),
                       "--reference-genome", selected_reference, "--breakpoint-window-genes", str(window_genes)]
            if no_plots:
                command.append("--no-plots")
            _run(command, f"Genome order analysis for {group_id}")
            group_run_rows.append({"similarity_group": group_id, "n_genomes": len(members),
                                   "reference_genome": selected_reference, "order_analysis_dir": str(group_out)})
        with (order_dir / "group_runs.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["similarity_group", "n_genomes", "reference_genome", "order_analysis_dir"],
                                    delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(group_run_rows)
        (outdir / "cohort_analysis_summary.json").write_text(json.dumps({
            "genomes": len(rows), "pangenome_family_map": str(family_map),
            "order_analysis_dir": str(order_dir), "per_similarity_group": True,
            "group_runs": len(group_run_rows),
            "similarity_stratification": str(groups_path),
            "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }, indent=2) + "\n", encoding="utf-8")
        return order_dir
    if expected.exists() and not force:
        raise FileExistsError(f"Order-analysis output already exists: {expected}; pass --force to replace")
    # Use the managed staging directory so order analysis sees the exact same
    # genome IDs as the global pangenome, even if source GFF/FASTA stems differ
    # only by case or compressed GFFs needed normalization.
    input_dir = manifest_path.parent / config["pangenome_inputs"]
    selected_reference = reference_genome or max(rows, key=_reference_score)["genome"]
    command = [order_executable, "--gff-dir", str(input_dir), "--fasta-dir", str(input_dir),
               "--outdir", str(order_dir), "--family-map", str(family_map), "--max-genomes-plotted", "50",
               "--reference-genome", selected_reference, "--breakpoint-window-genes", str(window_genes)]
    if no_plots:
        command.append("--no-plots")
    _run(command, "Genome order and breakpoint analysis")
    if groups_path.exists():
        _stratify_breakpoints(order_dir, groups_path, outdir / "breakpoint_support_by_similarity_group.tsv", window_genes)
    (outdir / "cohort_analysis_summary.json").write_text(json.dumps({
        "genomes": len(rows), "pangenome_family_map": str(family_map), "order_analysis_dir": str(order_dir),
        "similarity_stratification": str(groups_path) if groups_path.exists() else None,
        "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }, indent=2) + "\n", encoding="utf-8")
    return order_dir


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def create_cohort_subset(manifest_path: Path, genome_list_path: Path, outdir: Path) -> Path:
    """Create an independently runnable cohort directory from selected genome IDs."""
    rows = _load_manifest(manifest_path)
    if not genome_list_path.is_file():
        raise FileNotFoundError(f"Selected-genome list not found: {genome_list_path}")
    with genome_list_path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(2048)
        handle.seek(0)
        first_line = next((line.strip() for line in sample.splitlines() if line.strip()), "")
        if "\t" in sample or first_line.casefold() == "genome":
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames or "genome" not in reader.fieldnames:
                raise ValueError("Selected-genome TSV must contain a 'genome' column")
            selected = [(row.get("genome") or "").strip() for row in reader]
        else:
            selected = [line.strip() for line in handle if line.strip() and not line.lstrip().startswith("#")]
    selected_keys = {_normalise_genome_name(name) for name in selected if name}
    if not selected_keys:
        raise ValueError("Selected-genome list is empty")
    if len(selected_keys) != len([name for name in selected if name]):
        raise ValueError("Selected-genome list contains duplicate IDs")
    rows_by_key = {_normalise_genome_name(row["genome"]): row for row in rows}
    unknown = selected_keys - set(rows_by_key)
    if unknown:
        raise ValueError(f"Selected genomes are absent from the manifest: {sorted(unknown)[:10]}")
    subset_rows = [rows_by_key[key] for key in sorted(selected_keys)]
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"Subset output directory is not empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    staging = outdir / "pangenome_inputs"
    staging.mkdir()
    (staging / ".pyamilyseq-cohort-stage").write_text("Managed by PyamilySeq cohort workflow.\n", encoding="utf-8")
    for row in subset_rows:
        genome = row["genome"]
        fasta = Path(row["fasta_path"])
        gff_text = row.get("gff_path", "")
        if not fasta.is_file():
            raise FileNotFoundError(f"Genome FASTA is no longer available: {fasta}")
        fasta_suffix = ".fna.gz" if fasta.name.lower().endswith(".gz") else ".fna"
        (staging / f"{genome}{fasta_suffix}").symlink_to(fasta.resolve())
        if gff_text:
            gff = Path(gff_text)
            if not gff.is_file():
                raise FileNotFoundError(f"Genome GFF is no longer available: {gff}")
            staged_gff = staging / f"{genome}.gff3"
            if gff.name.lower().endswith(".gz"):
                with gzip.open(gff, "rb") as source, staged_gff.open("wb") as target:
                    shutil.copyfileobj(source, target)
            else:
                staged_gff.symlink_to(gff.resolve())
    manifest = outdir / "genome_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(subset_rows)
    source_metrics = manifest_path.parent / "assembly_metrics.tsv"
    metric_rows = []
    if source_metrics.is_file():
        with source_metrics.open("r", encoding="utf-8-sig", newline="") as source:
            metric_reader = csv.DictReader(source, delimiter="\t")
            metric_rows = [row for row in metric_reader if _normalise_genome_name(row.get("genome", "")) in selected_keys]
    if metric_rows:
        with (outdir / "assembly_metrics.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(metric_rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(metric_rows)
    parent_config_path = manifest_path.parent / "cohort_config.json"
    parent_config = json.loads(parent_config_path.read_text(encoding="utf-8")) if parent_config_path.is_file() else {}
    config = {
        "tool": "PyamilySeq cohort analysis", "version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "parent_manifest": str(manifest_path.resolve()), "selected_genomes": len(subset_rows),
        "fasta_dir": parent_config.get("fasta_dir", ""), "gff_dir": parent_config.get("gff_dir", ""),
        "metadata_tsv": "", "genomes": len(subset_rows), "manifest": manifest.name,
        "pangenome_inputs": staging.name,
    }
    (outdir / "cohort_config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def run_location_analysis(manifest_path: Path, outdir: Path, executable: str = "Track-Family-Locations",
                          family_map: Path | None = None, neighbors: int = 3,
                          minimum_mapping_rate: float = 0.0) -> Path:
    """Run the local family-neighborhood report using cohort inputs and global families."""
    if neighbors < 0 or not 0 <= minimum_mapping_rate <= 1:
        raise ValueError("Neighbor count must be non-negative and minimum mapping rate must be between 0 and 1")
    rows = _load_manifest(manifest_path)
    if any(not row.get("gff_path") for row in rows):
        raise ValueError("Location analysis requires a GFF annotation for every genome")
    config_path = manifest_path.parent / "cohort_config.json"
    if not config_path.is_file():
        raise FileNotFoundError(f"Cohort config missing; run cohort prepare first: {config_path}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    input_dir = manifest_path.parent / config["pangenome_inputs"]
    family_map = family_map or (outdir / "pangenome" / "gene_presence_absence.csv")
    if not family_map.is_file():
        raise FileNotFoundError(f"Family map not found; run pangenome first: {family_map}")
    location_dir = outdir / "location_analysis"
    command = [executable, "--clusters", str(family_map), "--cluster-format", "gpa",
               "--gff-dir", str(input_dir), "--outdir", str(location_dir), "--k", str(neighbors),
               "--minimum-mapping-rate", str(minimum_mapping_rate)]
    _run(command, "Gene-family neighborhood analysis")
    html_path = location_dir / "location_tracking.html"
    if not html_path.is_file():
        raise RuntimeError(f"Location analysis completed without its HTML report: {html_path}")
    return html_path


def build_cohort_html_report(manifest_path: Path, outdir: Path, report_path: Path | None = None) -> Path:
    """Build a self-contained interactive dashboard and genome selection tool."""
    manifest_rows = _load_manifest(manifest_path)
    cohort_root = manifest_path.parent
    order_dir = outdir / "order_analysis"
    order_genomes = {row["genome"]: row for row in _read_tsv_rows(order_dir / "genome_summary.tsv")}
    group_rows = _read_tsv_rows(outdir / "similarity" / "similarity_groups.tsv")
    groups = {row["genome"]: row.get("similarity_group", "") for row in group_rows}
    breakpoints = _read_tsv_rows(order_dir / "contig_breakpoints_vs_reference.tsv")
    endpoint_counts: dict[str, int] = defaultdict(int)
    endpoint_points = []
    contig_spans = []
    for row in breakpoints:
        left_rank = int(row.get("left_boundary_reference_gene_order") or 0)
        right_rank = int(row.get("right_boundary_reference_gene_order") or 0)
        if left_rank and right_rank:
            contig_spans.append({"genome": row["genome"], "contig": row.get("contig", ""),
                                 "reference_contig": row.get("reference_contig", ""),
                                 "start": min(left_rank, right_rank), "end": max(left_rank, right_rank)})
        for side in ("left", "right"):
            rank = int(row.get(f"{side}_boundary_reference_gene_order") or 0)
            if rank:
                endpoint_counts[row["genome"]] += 1
                endpoint_points.append({
                    "genome": row["genome"], "contig": row.get("contig", ""),
                    "reference_contig": row.get("reference_contig", ""),
                    "rank": rank, "side": side,
                    "family": row.get(f"{side}_boundary_anchor_family", ""),
                    "near_end_anchor": row.get(f"{side}_boundary_support_eligible", "").casefold() == "true",
                })
    candidates = _read_tsv_rows(order_dir / "contig_end_gene_candidates.tsv")
    end_gene_counts: dict[str, int] = defaultdict(int)
    family_end_counts: dict[str, set[str]] = defaultdict(set)
    for row in candidates:
        end_gene_counts[row["genome"]] += 1
        family = row.get("family", "")
        if family and (row.get("near_left_end", "").casefold() == "true" or
                       row.get("near_right_end", "").casefold() == "true"):
            family_end_counts[family].add(row["genome"])
    genome_data = []
    for manifest_row in manifest_rows:
        name = manifest_row["genome"]
        meta = {}
        try:
            meta = json.loads(manifest_row.get("metadata", "{}") or "{}")
        except json.JSONDecodeError:
            pass
        summary = order_genomes.get(name, {})
        record_metadata = _fasta_record_metadata(Path(manifest_row["fasta_path"]))
        genome_data.append({
            "genome": name, "n_contigs": int(manifest_row.get("n_contigs") or 0),
            "assembly_bp": int(manifest_row.get("assembly_bp") or 0),
            "n50_bp": int(manifest_row.get("n50_bp") or 0),
            "l50_contigs": int(manifest_row.get("l50_contigs") or 0),
            "n_genes": int(summary.get("n_genes") or 0),
            "n_endpoint_sides": endpoint_counts.get(name, 0),
            "n_end_candidate_genes": end_gene_counts.get(name, 0),
            "similarity_group": groups.get(name, ""), "metadata": meta,
            "length_basis": summary.get("length_basis", ""),
            "records": record_metadata,
            "physical_end_count": 2 * len(record_metadata),
        })
    support_rows = _read_tsv_rows(order_dir / "assembly_breakpoint_support.tsv")
    similarity_maps: dict[str, dict[str, float]] = {}
    similarity_matrix_path = outdir / "similarity" / "similarity-k31-scaled1000.csv"
    similarity_labels_path = outdir / "similarity" / "similarity-labels-k31-scaled1000.csv"
    if similarity_matrix_path.is_file():
        labels, matrix = _read_similarity_matrix(similarity_matrix_path, similarity_labels_path)
        similarity_maps = {seed: {other: value for other, value in zip(labels, values)}
                           for seed, values in zip(labels, matrix)}
    family_rows = []
    family_map_path = outdir / "pangenome" / "gene_presence_absence.csv"
    if family_map_path.is_file():
        outdir.mkdir(parents=True, exist_ok=True)
        with family_map_path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                family = row.get("Gene", "")
                try:
                    prevalence = int(row.get("No. isolates", 0) or 0)
                except ValueError:
                    prevalence = 0
                n_edge = len(family_end_counts.get(family, set()))
                if family:
                    family_rows.append({"family": family, "prevalence": prevalence, "n_edge_genomes": n_edge,
                                       "edge_fraction": n_edge / prevalence if prevalence else 0,
                                       "annotation": row.get("Annotation", "")})
        family_rows.sort(key=lambda row: (-row["n_edge_genomes"], -row["edge_fraction"], row["family"]))
        edge_summary_path = outdir / "family_end_association.tsv"
        with edge_summary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["family", "n_genomes_with_family",
                "n_genomes_with_near_end_family_member", "fraction_family_genomes_near_end", "annotation"],
                delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in family_rows:
                writer.writerow({"family": row["family"], "n_genomes_with_family": row["prevalence"],
                                 "n_genomes_with_near_end_family_member": row["n_edge_genomes"],
                                 "fraction_family_genomes_near_end": f"{row['edge_fraction']:.6f}",
                                 "annotation": row["annotation"]})
        family_rows = family_rows[:250]
    dataset = {
        "genomes": genome_data, "endpoints": endpoint_points, "contig_spans": contig_spans,
        "regions": [{"contig": row.get("reference_contig", ""),
                     "start": int(row.get("reference_gene_order_min") or 0),
                     "end": int(row.get("reference_gene_order_max") or 0),
                     "n_genomes": int(row.get("n_distinct_genomes") or 0),
                     "genomes": row.get("genomes", ""),
                     "description": row.get("interpretation", "")}
                    for row in support_rows],
        "families": family_rows,
        "groups": sorted({row["similarity_group"] for row in genome_data if row["similarity_group"]}),
        "similarities": similarity_maps,
        "provenance": {"manifest": str(manifest_path.resolve()), "cohort_root": str(cohort_root.resolve()),
                       "order_reports_available": bool(order_genomes), "pangenome_available": family_map_path.is_file(),
                       "location_report_available": (outdir / "location_analysis" / "location_tracking.html").is_file(),
                       "subset_command": ("PyamilySeq Analyze cohort subset --manifest " +
                           shlex.quote(str(manifest_path.resolve())) +
                           " --genomes-file selected_genomes.tsv --outdir selected_cohort")},
    }
    payload = json.dumps(dataset, ensure_ascii=True, separators=(",", ":")).replace("</", "<\\/")
    report_path = report_path or outdir / "cohort_report.html"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    html_doc = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>PyamilySeq cohort report</title><style>
:root{font:15px/1.45 system-ui,sans-serif;color:#17212b;background:#f4f7f9}body{margin:0}header{background:#102d3a;color:white;padding:24px 4vw}main{max-width:1500px;margin:20px auto;padding:0 18px}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px}.card,section{background:white;border:1px solid #d8e0e5;border-radius:10px;padding:16px;margin:14px 0}.card b{display:block;font-size:1.6rem}.controls{display:flex;gap:12px;flex-wrap:wrap;align-items:end}.controls label{display:flex;flex-direction:column;font-size:.85rem}input,select,button{font:inherit;padding:7px;border:1px solid #9cadb7;border-radius:5px}button{cursor:pointer;background:#e8f3f5}button.primary{background:#116276;color:white}table{width:100%;border-collapse:collapse;font-size:.88rem}th,td{text-align:left;border-bottom:1px solid #e2e8ec;padding:7px}th{position:sticky;top:0;background:#edf3f5}td.num{text-align:right;font-variant-numeric:tabular-nums}.scroll{max-height:480px;overflow:auto}canvas{width:100%;height:320px;border:1px solid #e3e9ec;border-radius:5px}.note{color:#52616b}.badge{background:#e4eff2;padding:2px 7px;border-radius:10px}.code{font-family:ui-monospace,monospace;background:#f1f3f4;padding:10px;overflow-wrap:anywhere}.two{display:grid;grid-template-columns:1fr 1fr;gap:12px}@media(max-width:850px){.two{grid-template-columns:1fr}}
</style></head><body><header><h1>PyamilySeq cohort report</h1><div>Assembly structure, reference-projected contig ends, family calls, and interactive cohort selection</div></header><main>
<div id="cards" class="cards"></div>
<section><h2>Genome structure</h2><p class="note">Every contig has two physical ends. Complete-record wording is retained as metadata, but topology is called circular only when explicitly stated. Plasmids are separate replicons.</p><canvas id="scatter" width="1200" height="260" aria-label="Contig count versus N50"></canvas></section>
<section><h2>Mapped contig ends across genomes</h2><p class="note">Each row is a genome; each tick is one physical contig end projected to reference-gene rank using a family anchor. Dashes are not silently converted into placements: unanchored ends remain in the physical-end counts, outside this map. Reference: <span id="refname"></span>.</p><canvas id="endpointMap" width="1200" height="500" aria-label="Genome by reference contig endpoint map"></canvas></section>
<section><h2>Select a sub-cohort</h2><div class="controls"><label>Search genomes or metadata<input id="search" type="search" placeholder="strain, phylogroup, metadata…"></label><label>Similarity group<select id="group"><option value="">All groups</option></select></label><label>Jaccard seed genome<select id="similaritySeed"><option value="">No seed filter</option></select></label><label>Minimum Jaccard <output id="thresholdText">0.90</output><input id="simThreshold" type="range" min="0.20" max="1.00" step="0.01" value="0.90"></label><label>Minimum contigs<input id="minContigs" type="number" min="0" value="0"></label><label>Maximum contigs<input id="maxContigs" type="number" min="0" placeholder="no limit"></label><label>Minimum N50 (bp)<input id="minN50" type="number" min="0" value="0"></label><button id="selectVisible">Select filtered</button><button id="clear">Clear selection</button><button id="download" class="primary">Download selected genome list</button></div><p><span id="counts"></span></p><div class="scroll"><table id="genomes"><thead><tr><th>Use</th><th>Genome</th><th>Similarity group</th><th>Sequence records</th><th>Physical ends</th><th>Mapped ends</th><th>Assembly bp</th><th>N50 bp</th><th>Genes</th><th>Metadata</th></tr></thead><tbody></tbody></table></div><p class="note">A seed filter retains genomes whose Sourmash Jaccard score to the selected seed meets the slider value. This is a screening control, not a taxonomic boundary. The download is a TSV with a <code>genome</code> column. Use the command below to create a separate cohort manifest and managed inputs.</p><div id="subsetCommand" class="code"></div></section>
<div class="two"><section><h2>Candidate recurrent assembly-end regions</h2><p class="note">Descriptive clusters only. Frequencies require callable denominators, lineage-aware null models, and sensitivity analysis before publication claims.</p><div class="scroll"><table><thead><tr><th>Reference contig</th><th>Gene range</th><th>Genomes</th><th>Members</th></tr></thead><tbody id="regions"></tbody></table></div></section><section><h2>Families with genes near contig ends</h2><p class="note">This summarizes called family members within the end-candidate report; it does not estimate false absence or prove that fragmentation changed a call.</p><div class="scroll"><table><thead><tr><th>Family</th><th>Genomes with family</th><th>Genomes with end-near member</th><th>End fraction</th><th>Annotation</th></tr></thead><tbody id="families"></tbody></table></div></section></div>
<section><h2>Further analysis</h2><div id="locationLink"></div><ol><li>Filter and select genomes above; download the selected genome list.</li><li>Create a self-contained sub-cohort manifest and managed inputs:</li></ol><div id="subsetCommand2" class="code"></div><ol start="3"><li>Run the subset through Sourmash, a new PyamilySeq pangenome, and order analysis. Family IDs will be recomputed within the subset, so compare methods using matched genomes and the same input annotations.</li></ol><p class="note">For a paper on fragmentation effects, benchmark draft assemblies and trusted complete references, and rerun each pangenome tool on the same genome subsets. Keep assembly method, annotation method, similarity stratum, and assembly quality as explicit covariates.</p></section>
</main><script id="report-data" type="application/json">__DATA__</script><script>
const D=JSON.parse(document.getElementById('report-data').textContent), selected=new Set(), $=id=>document.getElementById(id), esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt=n=>Number(n||0).toLocaleString(), groups=$('group'); D.groups.forEach(g=>groups.add(new Option(g,g)));
const seed=$('similaritySeed');Object.keys(D.similarities).sort().forEach(g=>seed.add(new Option(g,g)));
const n50s=D.genomes.map(g=>g.n50_bp).sort((a,b)=>a-b), median=n50s.length?n50s[Math.floor(n50s.length/2)]:0;
$('cards').innerHTML=[[D.genomes.length,'Genomes'],[fmt(D.genomes.reduce((a,g)=>a+g.n_contigs,0)),'Sequence records'],[fmt(D.genomes.reduce((a,g)=>a+g.assembly_bp,0)),'Assembly bp'],[fmt(median),'Median N50 bp'],[D.regions.length,'Candidate recurrent regions'],[D.families.length?fmt(D.families.length):'—','End-associated family rows']].map(x=>`<div class="card"><b>${x[0]}</b>${x[1]}</div>`).join('');
function filterRows(){let q=$('search').value.toLowerCase(), grp=groups.value,s=seed.value,threshold=+$('simThreshold').value,minC=+$('minContigs').value||0,maxC=+$('maxContigs').value||Infinity,minN=+$('minN50').value||0;return D.genomes.filter(g=>(!q||(g.genome+' '+JSON.stringify(g.metadata)).toLowerCase().includes(q))&&(!grp||g.similarity_group===grp)&&(!s||(D.similarities[s]?.[g.genome]??0)>=threshold)&&g.n_contigs>=minC&&g.n_contigs<=maxC&&g.n50_bp>=minN)}
function render(){let rows=filterRows();$('counts').textContent=`Showing ${rows.length} of ${D.genomes.length}; selected ${selected.size}.`; $('genomes').querySelector('tbody').innerHTML=rows.map(g=>`<tr><td><input type="checkbox" data-genome="${esc(g.genome)}" ${selected.has(g.genome)?'checked':''}></td><td>${esc(g.genome)}</td><td>${esc(g.similarity_group||'—')}</td><td class="num">${fmt(g.n_contigs)}</td><td class="num">${fmt(g.physical_end_count)}</td><td class="num">${fmt(g.n_endpoint_sides)}</td><td class="num">${fmt(g.assembly_bp)}</td><td class="num">${fmt(g.n50_bp)}</td><td class="num">${fmt(g.n_genes)}</td><td>${esc(JSON.stringify(g.metadata))}</td></tr>`).join('');$('genomes').querySelectorAll('input[data-genome]').forEach(c=>c.onchange=()=>{c.checked?selected.add(c.dataset.genome):selected.delete(c.dataset.genome);render()});}
['search','group','similaritySeed','minContigs','maxContigs','minN50','simThreshold'].forEach(id=>$(id).addEventListener(id==='search'||id==='simThreshold'?'input':'change',render));$('simThreshold').oninput=()=>{$('thresholdText').value=(+$('simThreshold').value).toFixed(2);render()};$('selectVisible').onclick=()=>{filterRows().forEach(g=>selected.add(g.genome));render()};$('clear').onclick=()=>{selected.clear();render()};$('download').onclick=()=>{if(!selected.size){alert('Select at least one genome.');return}let body='genome\n'+[...selected].sort().map(x=>x.replace(/[\t\r\n]/g,' ')).join('\n')+'\n',a=document.createElement('a');a.href=URL.createObjectURL(new Blob([body],{type:'text/tab-separated-values'}));a.download='selected_genomes.tsv';a.click();URL.revokeObjectURL(a.href)};
const cmd=D.provenance.subset_command;$('subsetCommand').textContent=cmd;$('subsetCommand2').textContent=cmd;
if(D.provenance.location_report_available)$('locationLink').innerHTML='<a href="location_analysis/location_tracking.html">Open family neighborhood and local synteny report</a>';else $('locationLink').textContent='Run the cohort locations stage to add a family-neighborhood report.';
$('regions').innerHTML=D.regions.sort((a,b)=>b.n_genomes-a.n_genomes).slice(0,100).map(r=>`<tr><td>${esc(r.contig)}</td><td>${fmt(r.start)}–${fmt(r.end)}</td><td>${fmt(r.n_genomes)}</td><td title="${esc(r.genomes)}">${esc(r.genomes.split(',').slice(0,4).join(', '))}${r.genomes.split(',').length>4?'…':''}</td></tr>`).join('')||'<tr><td colspan="4">Run order analysis to populate endpoint regions.</td></tr>';
$('families').innerHTML=D.families.slice(0,100).map(f=>`<tr><td>${esc(f.family)}</td><td>${fmt(f.prevalence)}</td><td>${fmt(f.n_edge_genomes)}</td><td>${(100*f.edge_fraction).toFixed(1)}%</td><td>${esc(f.annotation)}</td></tr>`).join('')||'<tr><td colspan="5">Run pangenome and order analyses to populate family summaries.</td></tr>';
function chart(){let c=$('scatter'),ctx=c.getContext('2d'),w=c.width,h=c.height,p=42;ctx.clearRect(0,0,w,h);let data=D.genomes.filter(g=>g.n_contigs>0),maxX=Math.max(1,...data.map(g=>g.n_contigs)),maxY=Math.max(1,...data.map(g=>g.n50_bp/1e6));ctx.strokeStyle='#aab7be';ctx.beginPath();ctx.moveTo(p,10);ctx.lineTo(p,h-p);ctx.lineTo(w-10,h-p);ctx.stroke();ctx.fillStyle='#42545e';ctx.fillText('N50 (Mb)',5,16);ctx.fillText('Sequence records',w/2,h-8);data.forEach(g=>{let x=p+(g.n_contigs/maxX)*(w-p-20),y=h-p-(g.n50_bp/1e6/maxY)*(h-p-24);ctx.beginPath();ctx.arc(x,y,4,0,Math.PI*2);ctx.fillStyle='#19758a';ctx.fill();});}
function endpointMap(){let c=$('endpointMap'),ctx=c.getContext('2d'),w=c.width,h=c.height,p={l:210,r:20,t:24,b:28},rows=[...D.genomes].sort((a,b)=>(a.similarity_group||'').localeCompare(b.similarity_group||'')||a.genome.localeCompare(b.genome)),pts=D.endpoints.filter(x=>x.rank>0),spans=D.contig_spans,max=Math.max(1,...pts.map(x=>x.rank)),dy=Math.min(18,(h-p.t-p.b)/Math.max(1,rows.length));ctx.clearRect(0,0,w,h);ctx.font='10px system-ui';ctx.strokeStyle='#d8e0e5';ctx.fillStyle='#263943';rows.forEach((g,i)=>{let y=p.t+(i+.5)*dy;ctx.fillText(g.genome.slice(0,30),4,y+3);ctx.beginPath();ctx.moveTo(p.l,y);ctx.lineTo(w-p.r,y);ctx.stroke();});ctx.strokeStyle='#6f8792';ctx.beginPath();ctx.moveTo(p.l,h-p.b);ctx.lineTo(w-p.r,h-p.b);ctx.stroke();ctx.fillStyle='#263943';ctx.fillText('Reference gene order',w/2,h-7);let idx=new Map(rows.map((g,i)=>[g.genome,i]));spans.forEach(s=>{let i=idx.get(s.genome);if(i===undefined)return;let y=p.t+(i+.5)*dy,x1=p.l+(s.start/max)*(w-p.l-p.r),x2=p.l+(s.end/max)*(w-p.l-p.r);ctx.strokeStyle='#76aeb5';ctx.lineWidth=2;ctx.beginPath();ctx.moveTo(x1,y);ctx.lineTo(x2,y);ctx.stroke();});ctx.lineWidth=1;pts.forEach(x=>{let i=idx.get(x.genome);if(i===undefined)return;let y=p.t+(i+.5)*dy,xx=p.l+(x.rank/max)*(w-p.l-p.r);ctx.strokeStyle='#c34f36';ctx.beginPath();ctx.moveTo(xx,y-4);ctx.lineTo(xx,y+4);ctx.stroke();});$('refname').textContent=pts.length?`${pts[0].reference_contig} (gene ranks 1–${max})`:'not available';}
render();chart();endpointMap();
</script></body></html>'''.replace("__DATA__", payload)
    report_path.write_text(html_doc, encoding="utf-8")
    return report_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="PyamilySeq-Cohort",
        description="Prepare cohorts, run pangenome and similarity analysis, and report reference-projected contig ends.")
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare", help="Validate paired assemblies/annotations and write a cohort manifest.")
    fasta_input = prepare.add_mutually_exclusive_group(required=True)
    fasta_input.add_argument("--fasta-dir", type=Path, help="Directory containing one multi-contig FASTA per genome.")
    fasta_input.add_argument("--combined-fasta", type=Path,
                             help="One FASTA containing records from multiple genomes; records are grouped by header label.")
    prepare.add_argument("--gff-dir", type=Path,
                         help="Optional directory of per-genome GFF3 annotations; required for pangenome and order analyses unless Bakta annotation is requested.")
    prepare.add_argument("--genome-header-regex",
                         help="Regex for --combined-fasta with a named (?P<genome>...) capture group. Defaults to BV-BRC/RefSeq bracket labels.")
    prepare.add_argument("--annotate-with-bakta", action="store_true",
                         help="Annotate split FASTAs with a local Bakta installation before preparing the cohort.")
    prepare.add_argument("--annotate-with-pyrodigal", action="store_true",
                         help="Predict CDS features per genome with Pyrodigal before preparing the cohort.")
    prepare.add_argument("--bakta", default="bakta", help="Bakta executable or path.")
    prepare.add_argument("--bakta-database", type=Path, help="Optional Bakta database directory.")
    prepare.add_argument("--bakta-threads", type=int, default=1)
    prepare.add_argument("--force", action="store_true", help="Replace existing Bakta annotations and managed cohort inputs.")
    prepare.add_argument("--metadata", type=Path, help="Optional TSV with one row per genome and a genome column.")
    prepare.add_argument("--outdir", type=Path, required=True)
    prepare.add_argument("--checksums", action="store_true", help="Compute SHA-256 over each FASTA while scanning it.")

    similarity = sub.add_parser("similarity", help="Build Sourmash signatures, matrix, and optional connected-component groups.")
    similarity.add_argument("--manifest", type=Path, required=True)
    similarity.add_argument("--outdir", type=Path, required=True)
    similarity.add_argument("--sourmash", default="sourmash")
    similarity.add_argument("--ksize", type=int, default=31)
    similarity.add_argument("--scaled", type=int, default=1000)
    similarity.add_argument("--min-similarity", type=float,
                            help="Optional Jaccard threshold (0-1); omit to inspect the matrix before choosing groups.")
    similarity.add_argument("--grouping-method", choices=("greedy_complete_linkage", "single_linkage"),
                            default="greedy_complete_linkage",
                            help="Complete-linkage groups require every pair to meet threshold; single-linkage can chain distant genomes.")
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

    locations = sub.add_parser("locations", help="Run gene-family neighborhood analysis and create HTML/SVG reports.")
    locations.add_argument("--manifest", type=Path, required=True)
    locations.add_argument("--outdir", type=Path, required=True)
    locations.add_argument("--family-map", type=Path)
    locations.add_argument("--executable", default="Track-Family-Locations")
    locations.add_argument("--neighbors", type=int, default=3)
    locations.add_argument("--minimum-mapping-rate", type=float, default=0.0)

    order = sub.add_parser("order", help="Run whole-cohort order analysis and stratify contig-end regions by Sourmash group.")
    order.add_argument("--manifest", type=Path, required=True)
    order.add_argument("--outdir", type=Path, required=True)
    order.add_argument("--family-map", type=Path)
    order.add_argument("--order-executable", default="Track-Genome-Order")
    order.add_argument("--reference-genome")
    order.add_argument("--window-genes", type=int, default=10)
    order.add_argument("--per-similarity-group", action="store_true",
                       help="Run independent local-reference order analyses for each Sourmash group.")
    order.add_argument("--no-plots", action="store_true")
    order.add_argument("--force", action="store_true")

    report = sub.add_parser("report", help="Build a self-contained interactive HTML cohort report.")
    report.add_argument("--manifest", type=Path, required=True)
    report.add_argument("--outdir", type=Path, required=True)
    report.add_argument("--report-path", type=Path, help="Optional HTML output path (defaults to OUTDIR/cohort_report.html).")

    subset = sub.add_parser("subset", help="Create a runnable cohort from selected genome IDs.")
    subset.add_argument("--manifest", type=Path, required=True)
    subset.add_argument("--genomes-file", type=Path, required=True,
                        help="TSV with a genome column or a plain one-genome-per-line text file.")
    subset.add_argument("--outdir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "prepare":
            if args.annotate_with_bakta and args.annotate_with_pyrodigal:
                raise ValueError("Choose only one annotation engine: --annotate-with-bakta or --annotate-with-pyrodigal")
            if args.combined_fasta:
                manifest = prepare_cohort_from_combined_fasta(
                    args.combined_fasta, args.gff_dir, args.outdir,
                    metadata_path=args.metadata, checksums=args.checksums,
                    genome_header_regex=args.genome_header_regex, annotate=args.annotate_with_bakta,
                    annotate_pyrodigal_enabled=args.annotate_with_pyrodigal, bakta=args.bakta,
                    bakta_database=args.bakta_database, bakta_threads=args.bakta_threads, force=args.force,
                )
            else:
                gff_dir = args.gff_dir
                if args.annotate_with_bakta:
                    gff_dir = annotate_with_bakta(args.fasta_dir, args.outdir / "bakta_annotations",
                                                   args.bakta, args.bakta_database, args.bakta_threads, args.force)
                elif args.annotate_with_pyrodigal:
                    gff_dir = annotate_with_pyrodigal(args.fasta_dir, args.outdir / "pyrodigal_annotations", args.force)
                manifest = prepare_cohort(args.fasta_dir, gff_dir, args.outdir, args.metadata, args.checksums)
            print(f"Prepared {sum(1 for _ in manifest.open()) - 1} genomes: {manifest}")
        elif args.command == "similarity":
            matrix, signatures = run_sourmash(args.manifest, args.outdir, args.sourmash, args.ksize,
                                               args.scaled, args.min_similarity, args.force, args.grouping_method)
            print(f"Similarity matrix: {matrix}\nSignatures: {signatures}")
        elif args.command == "pangenome":
            gpa = run_pangenome(args.manifest, args.outdir, args.executable, args.threads, args.memory_mb,
                                args.identity, args.length_diff, args.force)
            print(f"Global pangenome complete: {gpa}")
        elif args.command == "locations":
            html_path = run_location_analysis(args.manifest, args.outdir, args.executable,
                                              args.family_map, args.neighbors, args.minimum_mapping_rate)
            print(f"Gene-family neighborhood report: {html_path}")
        elif args.command == "order":
            order_dir = run_order_analysis(args.manifest, args.outdir, args.family_map, args.order_executable,
                                           args.reference_genome, args.no_plots,
                                           args.window_genes, args.force, args.per_similarity_group)
            print(f"Cohort order analysis complete: {order_dir}")
        elif args.command == "report":
            report_path = build_cohort_html_report(args.manifest, args.outdir, args.report_path)
            print(f"Interactive cohort report: {report_path}")
        else:
            subset_manifest = create_cohort_subset(args.manifest, args.genomes_file, args.outdir)
            print(f"Selected cohort prepared: {subset_manifest}")
        return 0
    except (FileNotFoundError, NotADirectoryError, FileExistsError, ValueError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
