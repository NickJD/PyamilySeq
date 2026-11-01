"""Nucleotide-based placement of assembly contigs against explicit references.

This module keeps minimap2's PAF coordinates and raw outputs intact. It reports
placement evidence only; it does not infer biological causes or project physical
endpoints.
"""
from __future__ import annotations

import csv
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

from .order_analysis import FASTA_SUFFIXES, _normalise_genome_name, _open_text, _stem_without_suffix

PRESETS = {"asm5", "asm10", "asm20"}
BLOCK_FIELDS = [
    "subcohort", "reference_genome", "query_genome", "query_contig", "query_length",
    "query_start_0based", "query_end_0based", "reference_contig", "reference_replicon_id",
    "reference_replicon_type", "reference_replicon_topology", "reference_start_0based",
    "reference_end_0based", "strand", "matching_bases", "alignment_block_length",
    "sequence_identity", "mapping_quality", "primary_secondary_status", "CIGAR", "cs_tag",
    "minimap2_preset",
]
SUMMARY_FIELDS = [
    "subcohort", "genome", "contig", "contig_length", "n_alignment_blocks",
    "n_reference_contigs", "n_reference_replicons", "n_primary_alignments",
    "n_secondary_alignments", "aligned_query_fraction", "best_mapping_quality",
    "placement_status", "placement_method", "mapping_warnings",
]
REFERENCE_FIELDS = ["internal_reference_id", "subcohort", "reference_genome", "original_fasta_header",
                    "original_contig_id", "replicon_id", "replicon_name", "replicon_type",
                    "replicon_topology", "sequence_length"]
QUERY_FIELDS = ["internal_query_id", "subcohort", "genome", "original_fasta_header", "original_contig_id", "sequence_length"]


def read_fasta_records(path: Path) -> list[dict]:
    """Read non-empty records while preserving full headers and duplicate IDs."""
    records = []
    header = None
    chunks = []

    def emit():
        if header is None:
            return
        sequence = "".join(chunks)
        if sequence:
            records.append({"header": header, "contig": header.split(maxsplit=1)[0], "sequence": sequence,
                            "length": len(sequence)})

    with _open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                emit()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header in {path}")
                chunks = []
            elif line:
                if header is None:
                    raise ValueError(f"FASTA sequence data appears before the first header in {path}")
                chunks.append(line)
    emit()
    return records


def _fasta_for(genome: str, fasta_dir: Path) -> Path | None:
    matches = [p for p in fasta_dir.rglob("*") if p.is_file() and p.name.lower().endswith(FASTA_SUFFIXES)
               and _normalise_genome_name(_stem_without_suffix(p, FASTA_SUFFIXES)) == _normalise_genome_name(genome)]
    if len(matches) > 1:
        raise ValueError(f"Multiple FASTA files match genome {genome!r}: {', '.join(map(str, matches))}")
    return matches[0] if matches else None


def parse_paf(path: Path) -> list[dict]:
    """Parse PAF mandatory columns and preserve minimap2 alignment tags."""
    records = []
    with path.open("r", encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF row {path}:{line_no}; expected at least 12 tab-delimited fields")
            try:
                qlen, qs, qe = map(int, fields[1:4])
                rlen, rs, re_ = map(int, fields[6:9])
                matches, block_len, mapq = map(int, fields[9:12])
            except ValueError as exc:
                raise ValueError(f"Malformed numeric PAF field at {path}:{line_no}") from exc
            if not (0 <= qs <= qe <= qlen and 0 <= rs <= re_ <= rlen) or fields[4] not in {"+", "-"}:
                raise ValueError(f"Invalid PAF coordinates or strand at {path}:{line_no}")
            tags = {}
            for tag in fields[12:]:
                parts = tag.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            records.append({"query_id": fields[0], "query_length": qlen, "query_start": qs, "query_end": qe,
                "reference_id": fields[5], "reference_length": rlen, "reference_start": rs, "reference_end": re_,
                "strand": fields[4], "matches": matches, "block_length": block_len, "mapq": mapq,
                "status": "primary" if tags.get("tp", "P") == "P" else "secondary",
                "cigar": tags.get("cg", ""), "cs": tags.get("cs", "")})
    return records


def _write_fasta(records: list[dict], path: Path) -> None:
    with path.open("w", encoding="utf-8") as out:
        for record in records:
            out.write(f">{record['internal_id']}\n")
            seq = record["sequence"]
            for offset in range(0, len(seq), 80):
                out.write(seq[offset:offset + 80] + "\n")


def _union_fraction(intervals: list[tuple[int, int]], length: int) -> float:
    if length <= 0:
        return 0.0
    intervals = sorted(intervals)
    total = 0
    left = right = None
    for start, end in intervals:
        if left is None:
            left, right = start, end
        elif start <= right:
            right = max(right, end)
        else:
            total += right - left
            left, right = start, end
    if left is not None:
        total += right - left
    return min(1.0, total / length)


def _placement_status(blocks: list[dict], query_length: int, min_fraction: float,
                      min_mapq: int, min_block_length: int) -> tuple[str, float, int]:
    if not blocks:
        return "unplaced", 0.0, 0
    fraction = _union_fraction([(b["query_start"], b["query_end"]) for b in blocks], query_length)
    usable = [b for b in blocks if b["block_length"] >= min_block_length]
    effective_mapqs = [b["mapq"] if b["mapq"] < 255 else 0 for b in usable]
    best_mapq = max(effective_mapqs, default=0)
    # Near-equivalent alignments that cover the same query interval indicate a
    # competing placement, not a second biological segment.
    for i, first in enumerate(usable):
        for second in usable[i + 1:]:
            overlap = max(0, min(first["query_end"], second["query_end"]) - max(first["query_start"], second["query_start"]))
            shorter = min(first["query_end"] - first["query_start"], second["query_end"] - second["query_start"])
            different_locus = (first["reference_id"] != second["reference_id"] or
                               abs(first["reference_start"] - second["reference_start"]) > max(100, shorter // 10))
            similar = min(first["matches"], second["matches"]) >= 0.9 * max(first["matches"], second["matches"])
            if shorter > 0 and overlap / shorter >= 0.8 and different_locus and similar:
                return "ambiguous", fraction, best_mapq
    # Construct a non-overlapping set greedily by alignment support. This retains
    # all blocks in the output while allowing a split status for disjoint segments.
    chosen = []
    for block in sorted(usable, key=lambda b: (-b["matches"], -(b["query_end"] - b["query_start"]))):
        if all(min(block["query_end"], other["query_end"]) <= max(block["query_start"], other["query_start"])
               for other in chosen):
            chosen.append(block)
    if len(chosen) >= 2 and _union_fraction([(b["query_start"], b["query_end"]) for b in chosen], query_length) >= min_fraction:
        loci = {(b["reference_id"], b["reference_start"] // 10000, b["strand"]) for b in chosen}
        if len(loci) >= 2:
            return "split", fraction, best_mapq
    if fraction >= min_fraction and best_mapq >= min_mapq and usable:
        return "placed", fraction, best_mapq
    return "low_confidence", fraction, best_mapq


def run_nucleotide_placement(genomes: list, fasta_dir: Path, outdir: Path,
                             assignments: dict[str, str], selected_references: dict[str, str],
                             presets: dict[str, str], replicon_metadata: list[dict],
                             minimap2: str = "minimap2", threads: int = 1,
                             min_aligned_fraction: float = 0.8, min_mapq: int = 20,
                             min_block_length: int = 1000) -> dict:
    """Run minimap2 independently per sub-cohort and write PAF, block and summary outputs."""
    if not selected_references:
        _write_tsv(outdir / "contig_placement_blocks.tsv", [], BLOCK_FIELDS)
        _write_tsv(outdir / "contig_placement_summary.tsv", [], SUMMARY_FIELDS)
        _write_tsv(outdir / "placement_reference_sequence_manifest.tsv", [], REFERENCE_FIELDS)
        _write_tsv(outdir / "placement_query_sequence_manifest.tsv", [], QUERY_FIELDS)
        return {"version": "", "n_blocks": 0, "n_contigs": 0, "warnings": ["No reference selected; nucleotide placement was not run."]}
    if fasta_dir is None:
        _write_tsv(outdir / "contig_placement_blocks.tsv", [], BLOCK_FIELDS)
        _write_tsv(outdir / "contig_placement_summary.tsv", [], SUMMARY_FIELDS)
        _write_tsv(outdir / "placement_reference_sequence_manifest.tsv", [], REFERENCE_FIELDS)
        _write_tsv(outdir / "placement_query_sequence_manifest.tsv", [], QUERY_FIELDS)
        return {"version": "", "n_blocks": 0, "n_contigs": 0,
                "warnings": ["Reference selected without --fasta-dir; nucleotide placement was not run."]}
    if threads < 1 or not 0 < min_aligned_fraction <= 1 or min_mapq < 0 or min_block_length < 1:
        raise ValueError("Placement threads and minimum block length must be positive; aligned fraction must be in (0,1]; MAPQ cannot be negative")
    executable = shutil.which(minimap2) if not Path(minimap2).is_file() else str(Path(minimap2).resolve())
    if not executable:
        raise FileNotFoundError(f"minimap2 executable not found: {minimap2}")
    version_result = subprocess.run([executable, "--version"], capture_output=True, text=True)
    if version_result.returncode:
        raise RuntimeError(f"minimap2 --version failed ({version_result.returncode}): {version_result.stderr.strip()}")
    version = version_result.stdout.strip() or version_result.stderr.strip()
    if not version:
        raise RuntimeError("minimap2 --version returned no version text")

    by_key = {_normalise_genome_name(g.name): g for g in genomes}
    meta_lists: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in replicon_metadata:
        meta_lists[(_normalise_genome_name(row["genome"]), row["contig"])].append(row)
    block_rows, summary_rows, ref_manifest, query_manifest = [], [], [], []
    warnings = []
    placement_root = outdir / "placement"
    for group_i, (subcohort, reference_name) in enumerate(sorted(selected_references.items(), key=lambda x: x[0].casefold()), 1):
        ref_key = _normalise_genome_name(reference_name)
        ref_path = _fasta_for(reference_name, fasta_dir)
        if not ref_path:
            raise FileNotFoundError(f"No FASTA assembly found for reference genome {reference_name!r}")
        ref_records = read_fasta_records(ref_path)
        if not ref_records:
            raise ValueError(f"Reference FASTA contains no non-empty records: {ref_path}")
        for index, record in enumerate(ref_records, 1):
            record["internal_id"] = f"R{index:06d}"
            record["original_id"] = record["contig"]
            record["metadata"] = meta_lists.get((ref_key, record["contig"]), []).pop(0) if meta_lists.get((ref_key, record["contig"])) else {}
            md = record["metadata"]
            ref_manifest.append({"internal_reference_id": record["internal_id"], "subcohort": subcohort,
                "reference_genome": reference_name, "original_fasta_header": record["header"],
                "original_contig_id": record["original_id"], "replicon_id": md.get("replicon_id", "unknown"),
                "replicon_name": md.get("replicon_name", ""), "replicon_type": md.get("replicon_type", "unknown"),
                "replicon_topology": md.get("replicon_topology", "unknown"), "sequence_length": record["length"]})
        safe_group = re.sub(r"[^A-Za-z0-9_.-]+", "_", subcohort).strip("._") or "cohort"
        group_dir = placement_root / f"{group_i:03d}_{safe_group}"
        group_dir.mkdir(parents=True, exist_ok=True)
        ref_fasta = group_dir / "reference_internal.fasta"
        _write_fasta(ref_records, ref_fasta)
        preset = presets.get(subcohort, "asm10")
        if preset not in PRESETS:
            raise ValueError(f"Unsupported minimap2 preset {preset!r} for subcohort {subcohort!r}; choose asm5, asm10, or asm20")
        cohort_genomes = [g for g in genomes if assignments.get(_normalise_genome_name(g.name)) == subcohort]
        for genome in sorted(cohort_genomes, key=lambda g: g.name.casefold()):
            query_key = _normalise_genome_name(genome.name)
            query_path = _fasta_for(genome.name, fasta_dir)
            query_records = read_fasta_records(query_path) if query_path else []
            if not query_records:
                continue
            for index, record in enumerate(query_records, 1):
                record["internal_id"] = f"Q{index:06d}"
                query_manifest.append({"internal_query_id": record["internal_id"], "subcohort": subcohort,
                    "genome": genome.name, "original_fasta_header": record["header"],
                    "original_contig_id": record["contig"], "sequence_length": record["length"]})
                # Metadata is consumed in FASTA record order to handle duplicate IDs safely.
                same = meta_lists.get((query_key, record["contig"]), [])
                record["metadata"] = same.pop(0) if same else {}
            if query_key == ref_key:
                paf_records = []
                for record in query_records:
                    target = next((r for r in ref_records if r["original_id"] == record["contig"] and r["sequence"] == record["sequence"]), None)
                    if target is None:
                        continue
                    paf_records.append({"query_id": record["internal_id"], "query_length": record["length"],
                        "query_start": 0, "query_end": record["length"], "reference_id": target["internal_id"],
                        "reference_length": target["length"], "reference_start": 0, "reference_end": target["length"],
                        "strand": "+", "matches": record["length"], "block_length": record["length"],
                        "mapq": 60, "status": "primary", "cigar": f"{record['length']}M", "cs": f":{record['length']}"})
                paf_path = group_dir / f"{re.sub(r'[^A-Za-z0-9_.-]+', '_', genome.name)}.identity.paf"
            else:
                query_fasta = group_dir / f"{re.sub(r'[^A-Za-z0-9_.-]+', '_', genome.name)}.internal.fasta"
                _write_fasta(query_records, query_fasta)
                paf_path = group_dir / f"{re.sub(r'[^A-Za-z0-9_.-]+', '_', genome.name)}.paf"
                command = [executable, "-x", preset, "-c", "--cs=long", "--secondary=yes", "-N", "1000",
                           "-t", str(threads), str(ref_fasta), str(query_fasta)]
                with paf_path.open("w", encoding="utf-8") as handle:
                    result = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True)
                if result.returncode:
                    raise RuntimeError(f"minimap2 failed for {genome.name!r} vs {reference_name!r} ({result.returncode}): {result.stderr.strip()}")
                paf_records = parse_paf(paf_path)
            query_lookup = {r["internal_id"]: r for r in query_records}
            ref_lookup = {r["internal_id"]: r for r in ref_records}
            by_query: dict[str, list[dict]] = defaultdict(list)
            for paf in paf_records:
                by_query[paf["query_id"]].append(paf)
                q = query_lookup.get(paf["query_id"])
                r = ref_lookup.get(paf["reference_id"])
                if not q or not r:
                    raise ValueError(f"PAF contains an unknown internal sequence identifier in {paf_path}")
                md = r["metadata"]
                block_rows.append({"subcohort": subcohort, "reference_genome": reference_name,
                    "query_genome": genome.name, "query_contig": q["contig"], "query_length": paf["query_length"],
                    "query_start_0based": paf["query_start"], "query_end_0based": paf["query_end"],
                    "reference_contig": r["original_id"], "reference_replicon_id": md.get("replicon_id", "unknown"),
                    "reference_replicon_type": md.get("replicon_type", "unknown"),
                    "reference_replicon_topology": md.get("replicon_topology", "unknown"),
                    "reference_start_0based": paf["reference_start"], "reference_end_0based": paf["reference_end"],
                    "strand": paf["strand"], "matching_bases": paf["matches"], "alignment_block_length": paf["block_length"],
                    "sequence_identity": paf["matches"] / paf["block_length"] if paf["block_length"] else 0,
                    "mapping_quality": paf["mapq"], "primary_secondary_status": paf["status"],
                    "CIGAR": paf["cigar"], "cs_tag": paf["cs"], "minimap2_preset": preset})
            for q in query_records:
                blocks = by_query.get(q["internal_id"], [])
                if query_key == ref_key:
                    status, fraction, best_mapq = "reference", 1.0, 60
                else:
                    status, fraction, best_mapq = _placement_status(blocks, q["length"], min_aligned_fraction,
                                                                    min_mapq, min_block_length)
                ref_contigs = {b["reference_id"] for b in blocks}
                ref_replicons = {ref_lookup[b["reference_id"]]["metadata"].get("replicon_id", "unknown")
                                 if ref_lookup[b["reference_id"]]["metadata"].get("replicon_id", "unknown") != "unknown"
                                 else f"unknown:{b['reference_id']}" for b in blocks}
                qtype = q.get("metadata", {}).get("replicon_type", "unknown")
                rtypes = {ref_lookup[b["reference_id"]]["metadata"].get("replicon_type", "unknown") for b in blocks}
                mapping_warnings = []
                if qtype in {"chromosome", "plasmid"} and any(t in {"chromosome", "plasmid"} and t != qtype for t in rtypes):
                    mapping_warnings.append("potential_replicon_type_mismatch")
                if any(t == "unknown" for t in rtypes):
                    mapping_warnings.append("reference_replicon_type_unknown")
                summary_rows.append({"subcohort": subcohort, "genome": genome.name, "contig": q["contig"],
                    "contig_length": q["length"], "n_alignment_blocks": len(blocks),
                    "n_reference_contigs": len(ref_contigs), "n_reference_replicons": len(ref_replicons),
                    "n_primary_alignments": sum(b["status"] == "primary" for b in blocks),
                    "n_secondary_alignments": sum(b["status"] == "secondary" for b in blocks),
                    "aligned_query_fraction": fraction, "best_mapping_quality": best_mapq,
                    "placement_status": status, "placement_method": "identity" if query_key == ref_key else "minimap2_nucleotide",
                    "mapping_warnings": ";".join(mapping_warnings)})
    outdir.mkdir(parents=True, exist_ok=True)
    _write_tsv(outdir / "contig_placement_blocks.tsv", block_rows, BLOCK_FIELDS)
    _write_tsv(outdir / "contig_placement_summary.tsv", summary_rows, SUMMARY_FIELDS)
    _write_tsv(outdir / "placement_reference_sequence_manifest.tsv", ref_manifest, REFERENCE_FIELDS)
    _write_tsv(outdir / "placement_query_sequence_manifest.tsv", query_manifest, QUERY_FIELDS)
    return {"version": version, "n_blocks": len(block_rows), "n_contigs": len(summary_rows), "warnings": warnings}


def _write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
