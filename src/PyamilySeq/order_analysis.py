"""Genome and gene-order analysis for PyamilySeq.

This module analyses GFF/GFF3 collections without changing the source annotations.
Comparative orientation normalisation is only attempted when a gene-family mapping
is supplied; GFF coordinates alone cannot establish homology or contig orientation.

Python >= 3.10. Core reports use the standard library. Figures require matplotlib.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import shutil
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, TextIO
from urllib.parse import unquote

VERSION = "0.1.0"
GFF_SUFFIXES = (".gff.gz", ".gff3.gz", ".gff", ".gff3", ".gtf.gz", ".gtf")
FASTA_SUFFIXES = (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa")
GENE_LIKE_TYPES = {"gene", "pseudogene", "ncRNA_gene", "rRNA", "tRNA", "tmRNA", "ncRNA", "RNA"}


@dataclass
class Feature:
    genome: str
    contig: str
    start: int
    end: int
    strand: str
    feature_type: str
    feature_id: str
    locus_tag: str = ""
    gene_name: str = ""
    product: str = ""
    family: str = ""
    label: str = ""
    attributes: dict[str, str] = field(default_factory=dict)
    gene_order_on_contig: int = 0
    gene_order_genome: int = 0
    normalized_gene_order: int = 0
    normalized_strand: str = "."
    bp_from_previous_gene: int | None = None
    cumulative_start: int = 0
    cumulative_end: int = 0
    orientation: str = "not_assessed"
    orientation_reference_contig: str = ""

    @property
    def length(self) -> int:
        return max(0, self.end - self.start + 1)

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2.0


@dataclass
class Contig:
    name: str
    length: int = 0
    length_source: str = "max_feature_end"
    input_order: int = 0
    genes: list[Feature] = field(default_factory=list)
    orientation: str = "not_assessed"
    reference_contig: str = ""
    orientation_anchors: int = 0
    kendall_tau_raw: float | None = None
    orientation_status: str = "not_assessed"
    replicon_id: str = ""
    replicon_name: str = ""
    replicon_type: str = "unknown"
    replicon_topology: str = "unknown"
    metadata_source: str = "unknown"
    topology_source: str = "unknown"
    classification_status: str = "unknown"
    has_fasta_record: bool = False
    has_gff_annotation: bool = False
    sequence_record_representation: str = "linear_fasta_record"
    sequence_role: str = ""
    assigned_molecule: str = ""
    assigned_molecule_location_type: str = ""
    assembly_unit: str = ""
    genbank_accession: str = ""
    refseq_accession: str = ""
    assembly_report_source: str = ""
    matched_assembly_accession: str = ""
    source_values: dict[str, str] = field(default_factory=dict)
    conflict_status: str = "none"
    resolution_explanation: str = "No explicit replicon metadata available."
    fasta_header: str = ""
    gff_region_evidence: list[dict] = field(default_factory=list)
    resolved_replicon_metadata: dict = field(default_factory=dict)


@dataclass
class Genome:
    name: str
    path: Path
    contigs: dict[str, Contig] = field(default_factory=dict)
    genes: list[Feature] = field(default_factory=list)
    selected_feature_type: str = ""
    warnings: list[str] = field(default_factory=list)


@dataclass
class FamilyMap:
    by_gene: dict[str, str] = field(default_factory=dict)
    by_genome_gene: dict[tuple[str, str], str] = field(default_factory=dict)
    ambiguous_gene_keys: set[str] = field(default_factory=set)
    source: str = ""
    gene_column: str = ""
    family_column: str = ""
    genome_column: str = ""
    is_presence_absence: bool = False

    def lookup(self, genome: str, keys: Iterable[str], delimiter: str = "|") -> str:
        norm_genome = _normalise_genome_name(genome)
        for key in keys:
            if not key:
                continue
            scoped = self.by_genome_gene.get((norm_genome, key))
            if scoped:
                return scoped
            # PyamilySeq sequence IDs commonly use genome|gene.
            if delimiter and delimiter in key:
                prefix, remainder = key.split(delimiter, 1)
                if _normalise_genome_name(prefix) == norm_genome:
                    scoped = self.by_genome_gene.get((norm_genome, remainder))
                    if scoped:
                        return scoped
            scoped = None if key in self.ambiguous_gene_keys else self.by_gene.get(key)
            if scoped:
                return scoped
            if delimiter and delimiter not in key:
                scoped = self.by_gene.get(f"{genome}{delimiter}{key}")
                if scoped:
                    return scoped
        return ""


def _normalise_genome_name(value: str) -> str:
    value = str(value).strip()
    for suffix in GFF_SUFFIXES + FASTA_SUFFIXES:
        if value.lower().endswith(suffix):
            value = value[: -len(suffix)]
            break
    return value.casefold()


def _open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _parse_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    if text == "." or not text:
        return attrs
    for part in text.split(";"):
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:  # tolerate simple GTF-like key "value" syntax
            key, value = part.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        key = unquote(key.strip())
        value = unquote(value.strip().strip('"'))
        if key:
            attrs[key] = value
    return attrs


def _stem_without_suffix(path: Path, suffixes: Iterable[str]) -> str:
    name = path.name
    for suffix in sorted(suffixes, key=len, reverse=True):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def discover_gffs(gff_dir: Path, pattern: str = "*") -> list[Path]:
    if not gff_dir.exists():
        raise FileNotFoundError(f"GFF directory does not exist: {gff_dir}")
    if not gff_dir.is_dir():
        raise NotADirectoryError(f"GFF input must be a directory: {gff_dir}")
    accepted = {s.lower() for s in GFF_SUFFIXES}
    paths = [p for p in gff_dir.rglob(pattern) if p.is_file() and any(p.name.lower().endswith(s) for s in accepted)]
    return sorted(paths, key=lambda p: (p.name.casefold(), str(p).casefold()))


def _choose_feature_type(rows: list[dict], requested: str) -> str:
    if requested.lower() != "auto":
        return requested
    types = {row["type"] for row in rows}
    if "CDS" in types:
        return "CDS"
    if "gene" in types:
        return "gene"
    for candidate in ("pseudogene", "ncRNA_gene", "rRNA", "tRNA", "tmRNA", "ncRNA", "RNA"):
        if candidate in types:
            return candidate
    available = ", ".join(sorted(types)) or "none"
    raise ValueError(
        "Could not identify genes in this GFF. No CDS/gene/RNA-like features were found "
        f"(available feature types: {available}). Set --feature-type to the desired column-3 type."
    )


def parse_gff(path: Path, requested_feature_type: str = "auto") -> Genome:
    genome_name = _stem_without_suffix(path, GFF_SUFFIXES)
    genome = Genome(name=genome_name, path=path)
    selected_rows: list[dict] = []
    available_types: set[str] = set()
    attrs_by_id: dict[tuple[str, str], dict[str, str]] = {}
    seq_region_lengths: dict[str, int] = {}
    region_lengths: dict[str, int] = {}
    seen_contigs: dict[str, int] = {}

    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        start, end = int(parts[2]), int(parts[3])
                        if start < 1 or end < start:
                            genome.warnings.append(f"Line {line_number}: invalid ##sequence-region coordinates {start}-{end}")
                            continue
                        seq_region_lengths[parts[1]] = end - start + 1
                        if parts[1] not in seen_contigs:
                            seen_contigs[parts[1]] = len(seen_contigs) + 1
                    except ValueError:
                        genome.warnings.append(f"Line {line_number}: invalid ##sequence-region directive")
                continue
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                genome.warnings.append(f"Line {line_number}: skipped GFF row with {len(parts)} columns (expected 9)")
                continue
            seqid, _source, ftype, start_s, end_s, _score, strand, _phase, attr_text = parts
            try:
                start, end = int(start_s), int(end_s)
            except ValueError:
                genome.warnings.append(f"Line {line_number}: skipped non-integer coordinates")
                continue
            if start < 1 or end < start:
                genome.warnings.append(f"Line {line_number}: skipped invalid 1-based coordinates {start}-{end}")
                continue
            attrs = _parse_attributes(attr_text)
            if seqid not in seen_contigs:
                seen_contigs[seqid] = len(seen_contigs) + 1
            if ftype.lower() == "region":
                region_lengths[seqid] = max(region_lengths.get(seqid, 0), end - start + 1)
                contig_meta = genome.contigs.setdefault(seqid, Contig(name=seqid, input_order=len(seen_contigs)))
                contig_meta.gff_region_evidence.append(dict(attrs))
            available_types.add(ftype)
            identifier = attrs.get("ID")
            if identifier:
                id_key = (seqid, identifier)
                existing_attrs = attrs_by_id.setdefault(id_key, {})
                for key, value in attrs.items():
                    if key in existing_attrs and existing_attrs[key] != value:
                        genome.warnings.append(
                            f"Line {line_number}: repeated ID {identifier!r} has conflicting {key!r} values; "
                            "the first value is retained for parent lookup"
                        )
                    else:
                        existing_attrs[key] = value
            # Retain only plausible selection candidates. This avoids holding every
            # annotation row when users select CDS or another explicit feature type.
            requested_lower = requested_feature_type.casefold()
            candidate = (
                requested_lower == "auto" and (ftype == "CDS" or ftype in GENE_LIKE_TYPES)
            ) or (requested_lower != "auto" and ftype == requested_feature_type)
            if candidate:
                selected_rows.append({
                    "seqid": seqid, "type": ftype, "start": start, "end": end,
                    "strand": strand if strand in {"+", "-", ".", "?"} else ".",
                    "attrs": attrs, "line_number": line_number,
                })

    # Contig lengths are selected later in this priority order: FASTA, sequence-region,
    # explicit region feature, then observed maximum feature end.
    for contig_name, first_seen in seen_contigs.items():
        length = seq_region_lengths.get(contig_name, 0)
        source = "sequence_region" if length else "max_feature_end"
        if not length and region_lengths.get(contig_name):
            length, source = region_lengths[contig_name], "region_feature"
        existing = genome.contigs.get(contig_name, Contig(name=contig_name))
        existing.length, existing.length_source, existing.input_order = length, source, first_seen
        existing.has_gff_annotation = True
        genome.contigs[contig_name] = existing

    feature_type = _choose_feature_type([{"type": t} for t in available_types], requested_feature_type)
    genome.selected_feature_type = feature_type
    selected = [r for r in selected_rows if r["type"] == feature_type]
    if not selected:
        raise ValueError(f"No features of type {feature_type!r} were found in {path}")

    # In NCBI-style GFF3, CDS rows often carry Parent=gene-ID while locus_tag/gene
    # are stored on the parent gene row. Inherit missing parent attributes to improve
    # stable identifiers and family-map matching without replacing CDS-specific fields.
    def parent_attributes(seqid: str, attrs: dict[str, str], visiting: set[str] | None = None) -> dict[str, str]:
        visiting = set() if visiting is None else visiting
        merged: dict[str, str] = {}
        current_id = attrs.get("ID", "")
        if current_id:
            if current_id in visiting:
                return merged
            visiting = visiting | {current_id}
        for parent_id in (x.strip() for x in attrs.get("Parent", "").split(",") if x.strip()):
            if parent_id in visiting:
                continue
            parent = attrs_by_id.get((seqid, parent_id), {})
            merged.update(parent_attributes(seqid, parent, visiting))
            merged.update({k: v for k, v in parent.items() if v})
        merged.update({k: v for k, v in attrs.items() if v})
        return merged
    enriched_selected = []
    for original in selected:
        row = dict(original)
        row["attrs"] = parent_attributes(row["seqid"], row["attrs"])
        enriched_selected.append(row)
    selected = enriched_selected

    # Collapse multipart records sharing an ID, while preserving distinct loci without IDs.
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in selected:
        attrs = row["attrs"]
        # Only stable feature IDs establish multipart identity. locus_tag, protein_id
        # and especially Name can repeat across distinct biological loci.
        identifier = attrs.get("ID")
        if identifier:
            group_key = (row["seqid"], f"id:{identifier}")
        else:
            group_key = (row["seqid"], f"coord:{row['start']}:{row['end']}:{row['strand']}:{row['line_number']}")
        grouped[group_key].append(row)

    for (seqid, group_key), parts in grouped.items():
        attrs: dict[str, str] = {}
        for row in parts:
            for key, value in row["attrs"].items():
                if value and (key not in attrs or not attrs[key]):
                    attrs[key] = value
        starts = [r["start"] for r in parts]
        ends = [r["end"] for r in parts]
        strands = {r["strand"] for r in parts if r["strand"] in {"+", "-"}}
        strand = next(iter(strands)) if len(strands) == 1 else ("?" if any(r["strand"] == "?" for r in parts) else ".")
        coord_id = f"{seqid}:{min(starts)}-{max(ends)}"
        feature_id = attrs.get("locus_tag") or attrs.get("ID") or attrs.get("protein_id") or attrs.get("Name") or coord_id
        gene = Feature(
            genome=genome.name, contig=seqid, start=min(starts), end=max(ends), strand=strand,
            feature_type=feature_type, feature_id=feature_id, locus_tag=attrs.get("locus_tag", ""),
            gene_name=attrs.get("gene", "") or attrs.get("Name", ""), product=attrs.get("product", ""),
            attributes=attrs,
        )
        genome.contigs[seqid].genes.append(gene)
        genome.genes.append(gene)

    for contig in genome.contigs.values():
        contig.genes.sort(key=lambda g: (g.start, g.end, g.feature_id.casefold()))
        if not contig.length:
            contig.length = max((g.end for g in contig.genes), default=0)
            contig.length_source = "max_feature_end"
        elif contig.genes and max(g.end for g in contig.genes) > contig.length:
            genome.warnings.append(
                f"Contig {contig.name}: gene end ({max(g.end for g in contig.genes)}) exceeds recorded length ({contig.length}); length retained"
            )
    return genome


def load_fasta_lengths(genome: Genome, fasta_dir: Path) -> tuple[dict[str, int], str | None]:
    target_name = _normalise_genome_name(genome.name)
    candidates = [p for p in fasta_dir.rglob("*") if p.is_file() and any(p.name.lower().endswith(s) for s in FASTA_SUFFIXES)]
    matching = [p for p in candidates if _normalise_genome_name(_stem_without_suffix(p, FASTA_SUFFIXES)) == target_name]
    if not matching:
        return {}, f"No FASTA filename matched genome {genome.name!r} in {fasta_dir}"
    if len(matching) > 1:
        return {}, (
            f"Multiple FASTA files match genome {genome.name!r}: "
            + ", ".join(str(p) for p in sorted(matching, key=lambda p: str(p).casefold()))
            + "; no FASTA lengths were applied"
        )
    path = matching[0]
    lengths: dict[str, int] = {}
    active = None
    length = 0
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if active is not None:
                    if active in lengths:
                        raise ValueError(f"Duplicate FASTA contig identifier {active!r} in {path}")
                    lengths[active] = length
                header = line[1:].strip()
                active = header.split()[0] if header else ""
                if not active:
                    raise ValueError(f"Empty FASTA header in {path}")
                contig = genome.contigs.get(active)
                if contig is None:
                    contig = Contig(name=active, input_order=len(genome.contigs))
                    genome.contigs[active] = contig
                contig.fasta_header = header
                # Use only explicit header fields. A complete-record claim does not
                # establish circular topology; it is retained as a classification hint.
                lower = header.casefold()
                header_type = "plasmid" if " plasmid" in lower else ("chromosome" if " chromosome" in lower else "unknown")
                if contig.replicon_type != "unknown" and header_type != "unknown" and contig.replicon_type != header_type:
                    genome.warnings.append(f"{active}: conflicting replicon type in FASTA header and GFF; retaining GFF type")
                elif contig.replicon_type == "unknown" and header_type != "unknown":
                    contig.replicon_type = header_type
                    contig.replicon_id = active
                    contig.replicon_name = header.split(maxsplit=1)[1] if len(header.split(maxsplit=1)) > 1 else active
                    contig.metadata_source = "FASTA_header"
                    contig.classification_status = "explicit_header"
                if "circular" in lower:
                    if contig.replicon_topology not in {"unknown", "circular"}:
                        genome.warnings.append(f"{active}: conflicting topology in FASTA header and GFF; retaining unknown")
                        contig.replicon_topology, contig.topology_source = "unknown", "conflict"
                    elif contig.replicon_topology == "unknown":
                        contig.replicon_topology, contig.topology_source = "circular", "FASTA_header"
                elif "linear" in lower and contig.replicon_topology == "unknown":
                    contig.replicon_topology, contig.topology_source = "linear", "FASTA_header"
                length = 0
            else:
                length += len(line.strip())
        if active is not None:
            if active in lengths:
                raise ValueError(f"Duplicate FASTA contig identifier {active!r} in {path}")
            lengths[active] = length
    return lengths, None


def _metadata_value(value: object) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"", "-", "na", "n/a", "unknown", "."} else text


def _normalize_replicon_type(value: object) -> str:
    text = _metadata_value(value).casefold().replace("_", "-")
    if not text:
        return ""
    if "unplaced" in text or "unlocalized" in text:
        return ""
    if "plasmid" in text:
        return "plasmid"
    if "chromosome" in text or text in {"src", "chromosomal"}:
        return "chromosome"
    if text in {"other", "mitochondrion", "organelle"}:
        return "other"
    return ""


def _normalize_topology(value: object) -> str:
    text = _metadata_value(value).casefold()
    if text in {"true", "yes", "1", "circular"}:
        return "circular"
    if text in {"false", "no", "0", "linear"}:
        return "linear"
    return ""


def resolve_replicon_metadata(contig_id: str, assembly_report: dict | None = None,
                              fasta_header: str | dict | None = None,
                              gff_region: dict | list[dict] | None = None) -> dict:
    """Resolve replicon metadata for one sequence record without mutating inputs.

    Evidence is normalized only for comparison; all original source values are
    returned verbatim. NCBI identity uses assigned molecule plus assembly unit.
    A GFF feature ID or FASTA record ID alone is not treated as replicon identity.
    """
    if isinstance(assembly_report, list) and not assembly_report:
        assembly_report = None
    if isinstance(gff_region, list) and not gff_region:
        gff_region = None
    if isinstance(assembly_report, list):
        results = [resolve_replicon_metadata(contig_id, row, fasta_header, gff_region)
                   for row in assembly_report]
        def merge_field(field):
            vals = {r[field] for r in results if r[field] not in {"", "unknown"}}
            return (next(iter(vals)), len(vals) > 1) if len(vals) == 1 else ("unknown", len(vals) > 1)
        typ, tc = merge_field("replicon_type")
        topology, toc = merge_field("replicon_topology")
        identity, ic = merge_field("replicon_id")
        conflicts = {x for r in results for x in r["conflict_status"].split(",") if x and x != "none"}
        conflicts.update(k for k, v in (("type", tc), ("topology", toc), ("identity", ic)) if v)
        base = results[0] if results else resolve_replicon_metadata(contig_id, None, fasta_header, gff_region)
        evidence = dict(base["metadata_evidence"])
        evidence["assembly_report"] = assembly_report
        status = "conflict" if conflicts else ("classified" if typ != "unknown" else base["classification_status"])
        return {**base, "replicon_id": identity, "replicon_type": typ, "replicon_topology": topology,
                "classification_status": status, "conflict_status": ",".join(sorted(conflicts)) or "none",
                "metadata_evidence": evidence,
                "metadata_sources": sorted({s for r in results for s in r["metadata_sources"]}),
                "resolution_notes": "; ".join(r["resolution_notes"] for r in results if r["resolution_notes"])}
    if isinstance(gff_region, list):
        results = [resolve_replicon_metadata(contig_id, assembly_report, fasta_header, region)
                   for region in gff_region]
        # Every row carries identical NCBI/header evidence; combine the GFF evidence
        # and only retain a property when every explicit result agrees.
        def merge(field):
            values = {r[field] for r in results if r[field] not in {"", "unknown"}}
            return (next(iter(values)), False) if len(values) == 1 else ("unknown", len(values) > 1)
        typ, tc = merge("replicon_type")
        topology, toc = merge("replicon_topology")
        identity, ic = merge("replicon_id")
        conflicts = set(x for r in results for x in r["conflict_status"].split(",") if x and x != "none")
        conflicts.update(k for k, v in (("type", tc), ("topology", toc), ("identity", ic)) if v)
        evidence = dict(results[0]["metadata_evidence"]) if results else {
            "assembly_report": dict(assembly_report or {}), "fasta_header": fasta_header}
        evidence["gff_region"] = gff_region
        status = "conflict" if conflicts else ("classified" if typ != "unknown" else "unknown")
        return {**(results[0] if results else {}), "contig_id": contig_id,
                "replicon_id": identity, "replicon_type": typ, "replicon_topology": topology,
                "classification_status": status, "conflict_status": ",".join(sorted(conflicts)) or "none",
                "metadata_sources": sorted({source for r in results for source in r["metadata_sources"]}),
                "metadata_evidence": evidence,
                "resolution_notes": "; ".join(r["resolution_notes"] for r in results if r["resolution_notes"])}
    ar = dict(assembly_report or {})
    gff = dict(gff_region or {})
    header_text = fasta_header if isinstance(fasta_header, str) else ""
    fh = dict(fasta_header) if isinstance(fasta_header, dict) else {}
    evidence = {"assembly_report": ar, "fasta_header": fasta_header, "gff_region": gff}
    source_values = {}
    type_values = []
    topology_values = []
    identity_values = []

    ar_role = _metadata_value(ar.get("Sequence-Role"))
    ar_location = _metadata_value(ar.get("Assigned-Molecule-Location/Type"))
    ar_molecule = _metadata_value(ar.get("Assigned-Molecule"))
    ar_unit = _metadata_value(ar.get("Assembly-Unit"))
    ar_status = str(ar.get("classification_status", "")).casefold()
    ar_type = "" if ar_status == "unplaced_scaffold" or "unplaced" in ar_role.casefold() else (
        _normalize_replicon_type(ar.get("replicon_type")) or _normalize_replicon_type(ar_location))
    if ar_type:
        type_values.append(("assembly_report", ar_type))
    if ar_role:
        source_values["sequence_role"] = {"assembly_report": ar_role}
    ar_identity = ""
    if "unplaced" not in ar_role.casefold() and ar_molecule and ar_type in {"chromosome", "plasmid"}:
        ar_identity = f"{ar_unit or 'assembly_unit_unknown'}:{ar_type}:{ar_molecule}"
        identity_values.append(("assembly_report", ar_identity))
    ar_topology = _normalize_topology(ar.get("replicon_topology"))
    if ar_topology:
        topology_values.append(("assembly_report", ar_topology))

    header_type = _normalize_replicon_type(fh.get("replicon_type"))
    if not header_type:
        lower = header_text.casefold()
        header_type = "plasmid" if " plasmid" in lower else ("chromosome" if " chromosome" in lower else "")
    if header_type:
        type_values.append(("fasta_header", header_type))
    header_topology = _normalize_topology(fh.get("replicon_topology") or fh.get("topology"))
    if not header_topology:
        header_topology = "circular" if "circular" in header_text.casefold() else (
            "linear" if "linear" in header_text.casefold() else "")
    if header_topology:
        topology_values.append(("fasta_header", header_topology))
    header_identity = _metadata_value(fh.get("replicon_id"))
    if header_identity:
        identity_values.append(("fasta_header", header_identity))

    gff_type = _normalize_replicon_type(gff.get("replicon_type") or gff.get("gbkey") or
                                        gff.get("molecule_type") or gff.get("Name"))
    if gff_type:
        type_values.append(("gff_region", gff_type))
    gff_topology = _normalize_topology(gff.get("Is_circular") or gff.get("is_circular") or gff.get("replicon_topology"))
    if gff_topology:
        topology_values.append(("gff_region", gff_topology))
    gff_identity = _metadata_value(gff.get("replicon_id") or gff.get("replicon_accession") or gff.get("accession"))
    if gff_identity:
        identity_values.append(("gff_region", gff_identity))

    def resolve(values):
        distinct = {value for _, value in values}
        return (next(iter(distinct)), False) if len(distinct) == 1 else ("unknown", len(distinct) > 1)

    resolved_type, type_conflict = resolve(type_values)
    resolved_topology, topology_conflict = resolve(topology_values)
    resolved_identity, identity_conflict = resolve(identity_values)
    # An unlocalized scaffold can have a known assigned molecule but not be
    # equated with an assembled replicon fragment without other identity evidence.
    if ar_role.casefold() == "unlocalized-scaffold" and ar_identity and not gff_identity and not header_identity:
        resolved_identity = ar_identity
    conflicts = [name for name, flag in (("type", type_conflict), ("topology", topology_conflict),
                                         ("identity", identity_conflict)) if flag]
    status = "conflict" if conflicts else ("classified" if resolved_type != "unknown" else "unknown")
    sources = sorted({source for values in (type_values, topology_values, identity_values) for source, _ in values})
    notes = []
    if conflicts:
        notes.append("Conflicting explicit evidence for: " + ", ".join(conflicts) + ".")
    if not resolved_identity or resolved_identity == "unknown":
        resolved_identity = "unknown"
        notes.append("No sufficiently reliable shared replicon identity was available.")
    if resolved_topology == "unknown":
        notes.append("Topology remains unknown because no unambiguous explicit topology evidence was available.")
    if ar_role.casefold() == "unplaced-scaffold":
        resolved_identity = "unknown"
        if not conflicts:
            status = "unplaced_scaffold"
        notes.append("NCBI marks this sequence as unplaced; no replicon assignment was made.")
    return {"contig_id": contig_id, "replicon_id": resolved_identity,
            "replicon_name": ar_molecule or _metadata_value(gff.get("replicon_name")) or _metadata_value(fh.get("replicon_name")),
            "replicon_type": resolved_type, "replicon_topology": resolved_topology,
            "sequence_role": ar_role or _metadata_value(gff.get("Sequence-Role")) or _metadata_value(fh.get("sequence_role")),
            "classification_status": status,
            "conflict_status": ",".join(conflicts) if conflicts else "none",
            "metadata_sources": sources, "metadata_evidence": evidence,
            "resolution_notes": " ".join(notes)}


def parse_assembly_report(path: Path) -> dict[str, dict[str, str]]:
    """Parse a local NCBI assembly report indexed by exact GenBank/RefSeq accession.

    Each accession maps to a structured row containing the original NCBI columns,
    normalized replicon classification, source path, and an explicit unknown
    topology (assembly reports do not normally establish circularity).
    """
    required = ("Sequence-Name", "Sequence-Role", "Assigned-Molecule",
                "Assigned-Molecule-Location/Type", "GenBank-Accn", "RefSeq-Accn",
                "Assembly-Unit", "Sequence-Length")
    records: dict[str, dict[str, str]] = {}
    by_sequence_name: dict[str, dict[str, str]] = {}
    found_header = False
    with _open_text(path) as handle:
        header = None
        for line_no, line in enumerate(handle, 1):
            if line.lstrip("\ufeff").startswith("# Sequence-Name\t"):
                header = line.lstrip("\ufeff")[2:].rstrip("\r\n").split("\t")
                missing = [name for name in required if name not in header]
                if missing:
                    raise ValueError(f"Assembly report {path} is missing columns: {', '.join(missing)}")
                if len(header) != len(set(header)):
                    raise ValueError(f"Assembly report {path} has duplicate column names")
                found_header = True
                continue
            if not line.strip() or line.startswith("#"):
                continue
            if header is None:
                raise ValueError(f"Assembly report {path} has data before its column header")
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(header):
                raise ValueError(f"Assembly report {path}:{line_no} has {len(values)} fields; expected {len(header)}")
            row = dict(zip(header, values))
            record = {key: row.get(key, "") for key in required}
            if not record["Sequence-Name"]:
                raise ValueError(f"Assembly report {path}:{line_no} has an empty Sequence-Name")
            if not record["Sequence-Length"].isdigit():
                raise ValueError(f"Assembly report {path}:{line_no} has invalid Sequence-Length {record['Sequence-Length']!r}")
            record["Sequence-Length"] = str(int(record["Sequence-Length"]))
            record["source"] = str(path)
            role = record["Sequence-Role"].casefold()
            location = record["Assigned-Molecule-Location/Type"].casefold()
            if "unlocalized" in role:
                record["replicon_type"] = "unknown"
                record["classification_status"] = "unlocalized_scaffold"
            elif "unplaced" in role:
                record["replicon_type"] = "unknown"
                record["classification_status"] = "unplaced_scaffold"
            elif "plasmid" in location:
                record["replicon_type"] = "plasmid"
                record["classification_status"] = "classified"
            elif "chromosome" in location:
                record["replicon_type"] = "chromosome"
                record["classification_status"] = "classified"
            elif record["Sequence-Role"]:
                record["replicon_type"] = "other"
                record["classification_status"] = "other_sequence"
            else:
                record["replicon_type"] = "unknown"
                record["classification_status"] = "unknown"
            record["replicon_topology"] = "unknown"
            previous_name = by_sequence_name.get(record["Sequence-Name"])
            if previous_name and previous_name != record:
                raise ValueError(f"Contradictory assembly-report records for Sequence-Name {record['Sequence-Name']!r} in {path}")
            by_sequence_name[record["Sequence-Name"]] = record
            for accession in (record["GenBank-Accn"], record["RefSeq-Accn"]):
                if accession and accession.strip() not in {"-", "na"}:
                    if accession in records:
                        kind = "contradictory" if records[accession] != record else "duplicate"
                        raise ValueError(f"{kind.capitalize()} assembly-report record for accession {accession!r} in {path}")
                    records[accession] = record
    if not found_header:
        raise ValueError(f"Assembly report {path} has no '# Sequence-Name' column header")
    if not records:
        raise ValueError(f"Assembly report {path} contains no records with GenBank or RefSeq accessions")
    return records


def _load_assembly_report_mapping(mapping_tsv: Path | None, genomes: list[Genome]) -> dict[str, Path]:
    if mapping_tsv is None:
        return {}
    if not mapping_tsv.is_file():
        raise FileNotFoundError(f"Assembly-report mapping TSV does not exist: {mapping_tsv}")
    known = {_normalise_genome_name(g.name): g.name for g in genomes}
    result: dict[str, Path] = {}
    with mapping_tsv.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"genome", "assembly_report"}.issubset(reader.fieldnames):
            raise ValueError("Assembly-report mapping TSV must have genome and assembly_report columns")
        for line_no, row in enumerate(reader, 2):
            genome_raw = (row.get("genome") or "").strip()
            report_raw = (row.get("assembly_report") or "").strip()
            key = _normalise_genome_name(genome_raw)
            if not key or not report_raw:
                raise ValueError(f"Assembly-report mapping TSV row {line_no} has an empty genome or assembly_report")
            if key not in known:
                raise ValueError(f"Assembly-report mapping names unknown genome {genome_raw!r} at row {line_no}")
            if known[key] in result:
                raise ValueError(f"Duplicate assembly-report mapping for genome {known[key]!r}")
            path = Path(report_raw)
            if not path.is_absolute():
                path = mapping_tsv.parent / path
            if not path.is_file():
                raise FileNotFoundError(f"Assembly report for genome {known[key]!r} does not exist: {path}")
            result[known[key]] = path.resolve()
    return result


def collect_contig_replicon_metadata(genomes: list[Genome], fasta_dir: Path,
                                     assembly_report_mapping: Path | None = None) -> list[dict]:
    """Return one resolved replicon metadata row per non-empty FASTA record.

    FASTA/GFF association uses exact genome stems and sequence IDs. Assembly-report
    association is explicitly genome-mapped, then accession-matched against exact
    header tokens. GFF-only records are intentionally omitted.
    """
    if not fasta_dir.is_dir():
        raise NotADirectoryError(f"FASTA input must be an existing directory: {fasta_dir}")
    report_paths = _load_assembly_report_mapping(assembly_report_mapping, genomes)
    results = []
    fasta_suffixes = tuple(s.lower() for s in FASTA_SUFFIXES)
    for genome in sorted(genomes, key=lambda g: g.name.casefold()):
        matches = [p for p in fasta_dir.rglob("*") if p.is_file() and p.name.lower().endswith(fasta_suffixes)
                   and _normalise_genome_name(_stem_without_suffix(p, FASTA_SUFFIXES)) == _normalise_genome_name(genome.name)]
        if len(matches) != 1:
            raise ValueError(f"Expected exactly one FASTA matching genome {genome.name!r}; found {len(matches)}")
        report_path = report_paths.get(genome.name)
        report_index = parse_assembly_report(report_path) if report_path else {}
        accession_pattern = re.compile(r"(?<![A-Za-z0-9])([A-Za-z]{1,6}_?\d+(?:\.\d+)?)(?![A-Za-z0-9])")
        with _open_text(matches[0]) as handle:
            header = None
            seq_id = ""
            length = 0
            def emit_record():
                if header is None or length <= 0:
                    return
                tokens = set(accession_pattern.findall(header))
                matched = []
                matched_genbank = set()
                matched_refseq = set()
                seen = set()
                for token in sorted(tokens):
                    rec = report_index.get(token)
                    if rec is not None:
                        if token == rec.get("GenBank-Accn"):
                            matched_genbank.add(token)
                        if token == rec.get("RefSeq-Accn"):
                            matched_refseq.add(token)
                        signature = tuple(sorted(rec.items()))
                        if signature not in seen:
                            seen.add(signature)
                            matched.append(rec)
                contig = genome.contigs.get(seq_id)
                gff_regions = contig.gff_region_evidence if contig else []
                resolved = resolve_replicon_metadata(seq_id, matched, header, gff_regions)
                evidence = dict(resolved["metadata_evidence"])
                evidence["assembly_report_mapping"] = str(report_path) if report_path else ""
                evidence["matched_accessions"] = sorted(tokens & set(report_index))
                return {"genome": genome.name, "contig": seq_id, "fasta_header": header,
                        "contig_length": length, "gff_annotation_status": "annotated" if contig and contig.has_gff_annotation else "unannotated",
                        "matched_genbank_accessions": sorted(matched_genbank),
                        "matched_refseq_accessions": sorted(matched_refseq),
                        "assembly_report_source": str(report_path) if report_path else "",
                        **{**resolved, "metadata_evidence": evidence}}
            for line in handle:
                if line.startswith(">"):
                    row = emit_record()
                    if row:
                        results.append(row)
                    header = line[1:].strip()
                    seq_id = header.split(maxsplit=1)[0] if header else ""
                    if not seq_id:
                        raise ValueError(f"Empty FASTA header in {matches[0]}")
                    length = 0
                else:
                    length += len(line.strip())
            row = emit_record()
            if row:
                results.append(row)
    return results


def apply_fasta_lengths(genomes: list[Genome], fasta_dir: Path | None,
                        assembly_report_dir: Path | None = None) -> list[str]:
    warnings: list[str] = []
    if fasta_dir is None:
        return warnings
    if not fasta_dir.exists() or not fasta_dir.is_dir():
        raise NotADirectoryError(f"--fasta-dir must be an existing directory: {fasta_dir}")
    for genome in genomes:
        gff_contigs = {name for name, contig in genome.contigs.items() if contig.has_gff_annotation or contig.genes}
        lengths, warning = load_fasta_lengths(genome, fasta_dir)
        if warning:
            warnings.append(warning)
            continue
        for contig_name, contig in genome.contigs.items():
            if contig_name in lengths:
                contig.length = lengths[contig_name]
                contig.length_source = "fasta"
                contig.has_fasta_record = True
            else:
                warnings.append(f"{genome.name}: contig {contig_name!r} absent from matching FASTA; retained GFF-derived length")
        extra_fasta_contigs = sorted(set(lengths) - gff_contigs, key=str.casefold)
        if extra_fasta_contigs:
            for contig_name in extra_fasta_contigs:
                genome.contigs[contig_name] = Contig(
                    name=contig_name, length=lengths[contig_name],
                    length_source="fasta_no_gff_features",
                    has_fasta_record=True,
                )
            warnings.append(
                f"{genome.name}: FASTA contains {len(extra_fasta_contigs)} contig(s) without GFF features; "
                "included in assembly contig counts and length/N50 summaries, but not family-anchor comparisons "
                f"(examples: {', '.join(extra_fasta_contigs[:10])})"
            )
    if assembly_report_dir is not None:
        report_paths = sorted(assembly_report_dir.rglob("assembly_report.txt")) if assembly_report_dir.is_dir() else [assembly_report_dir]
        report_index = {}
        for report_path in report_paths:
            report_index.update(parse_assembly_report(report_path))
        for genome in genomes:
            for contig in genome.contigs.values():
                report = report_index.get(contig.name)
                if not report:
                    continue
                contig.assembly_report_source = report["source"]
                contig.matched_assembly_accession = contig.name
                contig.sequence_role = report["Sequence-Role"]
                contig.assigned_molecule = report["Assigned-Molecule"]
                contig.assigned_molecule_location_type = report["Assigned-Molecule-Location/Type"]
                contig.assembly_unit = report["Assembly-Unit"]
                contig.genbank_accession = report["GenBank-Accn"]
                contig.refseq_accession = report["RefSeq-Accn"]
                contig.replicon_name = report["Assigned-Molecule"] if report["Assigned-Molecule"] != "-" else report["Sequence-Name"]
                role = report["Sequence-Role"].casefold()
                location = report["Assigned-Molecule-Location/Type"].casefold()
                if "unplaced" in role or "unlocalized" in role:
                    contig.replicon_type = "other"
                    contig.replicon_id = ""
                    contig.classification_status = "unplaced_sequence" if "unplaced" in role else "unlocalized_scaffold"
                else:
                    contig.replicon_id = report["Assigned-Molecule"] if report["Assigned-Molecule"] != "-" else ""
                    contig.replicon_type = "plasmid" if "plasmid" in location else ("chromosome" if "chromosome" in location else "other")
                    contig.classification_status = "assembly_report"
                contig.metadata_source = "NCBI_assembly_report"
    return warnings


def _open_delimited(path: Path) -> tuple[TextIO, csv.Dialect]:
    handle = path.open("r", encoding="utf-8-sig", newline="")
    suffix = path.suffix.casefold()
    if suffix == ".csv":
        return handle, csv.excel
    if suffix in {".tsv", ".tab"}:
        return handle, csv.excel_tab
    sample = handle.read(8192)
    handle.seek(0)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
    except csv.Error:
        dialect = csv.excel_tab
    return handle, dialect


def _find_column(fieldnames: list[str], candidates: Iterable[str]) -> str | None:
    lookup = {name.strip().casefold(): name for name in fieldnames if name is not None}
    for candidate in candidates:
        if candidate.casefold() in lookup:
            return lookup[candidate.casefold()]
    return None


def _register_unscoped_family(family_map: FamilyMap, gene_id: str, family_id: str) -> None:
    """Register an unscoped key only while it remains unambiguous across genomes."""
    if not gene_id or gene_id in family_map.ambiguous_gene_keys:
        return
    current = family_map.by_gene.get(gene_id)
    if current is None:
        family_map.by_gene[gene_id] = family_id
    elif current != family_id:
        family_map.by_gene.pop(gene_id, None)
        family_map.ambiguous_gene_keys.add(gene_id)


def _register_scoped_family(family_map: FamilyMap, genome: str, gene_id: str, family_id: str) -> None:
    key = (_normalise_genome_name(genome), gene_id)
    existing = family_map.by_genome_gene.get(key)
    if existing is not None and existing != family_id:
        raise ValueError(
            f"Conflicting family assignments for gene {gene_id!r} in genome {genome!r}: "
            f"{existing!r} and {family_id!r}"
        )
    family_map.by_genome_gene[key] = family_id


def load_family_map(path: Path, gene_column: str | None = None, family_column: str | None = None,
                    genome_column: str | None = None, delimiter: str = "|") -> FamilyMap:
    if not path.exists():
        raise FileNotFoundError(f"Family map does not exist: {path}")
    handle, dialect = _open_delimited(path)
    result = FamilyMap(source=str(path))
    with handle:
        reader = csv.DictReader(handle, dialect=dialect)
        headers = reader.fieldnames or []
        if not headers:
            raise ValueError(f"Family map has no header row: {path}")
        gene_col = gene_column or _find_column(headers, (
            "gene_id", "sequence_id", "member_id", "member", "gene_identifier", "identifier",
            "input_id", "sequence", "gene", "seq_id", "id",
        ))
        family_col = family_column or _find_column(headers, (
            "family_id", "family", "gene_family", "cluster_id", "cluster", "group_id", "group",
            "representative", "representative_id", "gene_group",
        ))
        genome_col = genome_column or _find_column(headers, ("genome", "isolate", "sample", "strain", "genome_id"))

        # PyamilySeq's wide gene_presence_absence.csv layout can be used directly:
        # one family per row and genome columns contain gene identifiers.
        if gene_col is None or family_col is None:
            first_col = headers[0]
            if len(headers) >= 2 and first_col.strip().casefold() in {"gene", "family", "family_id", "gene_family"}:
                result.is_presence_absence = True
                result.gene_column = first_col
                result.family_column = first_col
                metadata_cols = {
                    "gene", "family", "family_id", "gene_family", "non-unique gene name", "annotation",
                    "no. isolates", "no. sequences", "avg sequences per isolate", "genome fragment",
                    "order within fragment", "accessory fragment", "accessory order with fragment", "qc",
                    "min group size nuc", "max group size nuc", "avg group size nuc", "description",
                }
                for row in reader:
                    fam = (row.get(first_col) or "").strip()
                    if not fam:
                        continue
                    for col in headers[1:]:
                        if col.strip().casefold() in metadata_cols:
                            continue
                        value = (row.get(col) or "").strip()
                        if not value or value in {"-", ".", "NA", "na"}:
                            continue
                        # GPA cells list multiple paralogs separated by spaces;
                        # some third-party files instead use comma separators.
                        separators = "," if "," in value else None
                        members = [x.strip() for x in value.split(separators) if x.strip()] if separators else value.split()
                        for member in members:
                            scoped_genome = _normalise_genome_name(col)
                            _register_scoped_family(result, scoped_genome, member, fam)
                            composite = f"{col}{delimiter}{member}" if delimiter else member
                            result.by_gene[composite] = fam
                            _register_unscoped_family(result, member, fam)
                return result
            raise ValueError(
                f"Could not infer gene and family columns from {headers}. Pass --map-gene-column and --map-family-column."
            )

        result.gene_column, result.family_column = gene_col, family_col
        result.genome_column = genome_col or ""
        for row in reader:
            gene_value = (row.get(gene_col) or "").strip()
            family_value = (row.get(family_col) or "").strip()
            if not gene_value or not family_value:
                continue
            if genome_col and (row.get(genome_col) or "").strip():
                gname = (row.get(genome_col) or "").strip()
                _register_scoped_family(result, gname, gene_value, family_value)
            if delimiter and delimiter in gene_value:
                prefix, member = gene_value.split(delimiter, 1)
                _register_scoped_family(result, prefix, member, family_value)
                # The raw composite is safe as an exact match unless the file itself conflicts.
                if gene_value not in result.by_gene or result.by_gene[gene_value] == family_value:
                    result.by_gene[gene_value] = family_value
                else:
                    result.by_gene.pop(gene_value, None)
                    result.ambiguous_gene_keys.add(gene_value)
                _register_unscoped_family(result, member, family_value)
            else:
                if genome_col and (row.get(genome_col) or "").strip():
                    _register_unscoped_family(result, gene_value, family_value)
                else:
                    _register_unscoped_family(result, gene_value, family_value)
    return result


def _gene_lookup_keys(gene: Feature) -> list[str]:
    values = [gene.feature_id, gene.locus_tag, gene.attributes.get("ID", ""),
              gene.attributes.get("Parent", ""), gene.attributes.get("protein_id", ""),
              gene.attributes.get("Name", ""), gene.attributes.get("gene", "")]
    values.extend(x.strip() for x in gene.attributes.get("Alias", "").split(",") if x.strip())
    deduped = []
    seen = set()
    for value in values:
        if value and value not in seen:
            deduped.append(value)
            seen.add(value)
    return deduped


def assign_families(genomes: list[Genome], family_map: FamilyMap | None, delimiter: str = "|") -> int:
    if family_map is None:
        for genome in genomes:
            for gene in genome.genes:
                gene.family = ""
        return 0
    mapped = 0
    for genome in genomes:
        for gene in genome.genes:
            gene.family = family_map.lookup(genome.name, _gene_lookup_keys(gene), delimiter)
            if gene.family:
                mapped += 1
    return mapped


def set_labels(genomes: list[Genome], label_by: str) -> None:
    for genome in genomes:
        for gene in genome.genes:
            if label_by == "family":
                candidate = gene.family
            elif label_by == "locus_tag":
                candidate = gene.locus_tag
            elif label_by == "Name":
                candidate = gene.attributes.get("Name", "")
            elif label_by == "gene":
                candidate = gene.attributes.get("gene", "") or gene.gene_name
            elif label_by == "product":
                candidate = gene.product
            elif label_by == "ID":
                candidate = gene.attributes.get("ID", "") or gene.feature_id
            else:  # auto
                candidate = gene.family or gene.locus_tag or gene.attributes.get("gene", "") or gene.attributes.get("Name", "") or gene.feature_id
            gene.label = candidate or gene.feature_id


def _inversion_count(values: list[int]) -> int:
    """Count inversions in O(n log n), using a Fenwick tree over reference ranks."""
    if len(values) < 2:
        return 0
    ranks = {v: i + 1 for i, v in enumerate(sorted(set(values)))}
    tree = [0] * (len(ranks) + 1)

    def add(index: int) -> None:
        while index < len(tree):
            tree[index] += 1
            index += index & -index

    def prefix(index: int) -> int:
        total = 0
        while index > 0:
            total += tree[index]
            index -= index & -index
        return total

    inversions = 0
    inserted = 0
    for value in values:
        rank = ranks[value]
        inversions += inserted - prefix(rank)
        add(rank)
        inserted += 1
    return inversions


def kendall_tau_from_ranks(reference_ranks_in_query_order: list[int]) -> float | None:
    n = len(reference_ranks_in_query_order)
    if n < 2 or len(set(reference_ranks_in_query_order)) != n:
        return None
    pairs = n * (n - 1) / 2
    inversions = _inversion_count(reference_ranks_in_query_order)
    return (pairs - 2 * inversions) / pairs


def _single_copy_families(genome: Genome) -> dict[str, tuple[Contig, Feature]]:
    counts = Counter(g.family for g in genome.genes if g.family)
    result: dict[str, tuple[Contig, Feature]] = {}
    for contig in genome.contigs.values():
        for gene in contig.genes:
            if gene.family and counts[gene.family] == 1:
                result[gene.family] = (contig, gene)
    return result


def infer_contig_orientations(genomes: list[Genome], reference_name: str | None,
                              min_anchors: int = 5, min_abs_tau: float = 0.5) -> tuple[str | None, list[dict]]:
    if not genomes or not any(g.family for genome in genomes for g in genome.genes):
        return None, []
    if not reference_name:
        raise ValueError("Reference genome must be explicitly specified; no alphabetical default is selected")
    lookup = {g.name.casefold(): g for g in genomes}
    reference = lookup.get(reference_name.casefold())
    if reference is None:
        choices = ", ".join(g.name for g in genomes)
        raise ValueError(f"Reference genome {reference_name!r} not found among inputs. Available: {choices}")
    reference_unique = _single_copy_families(reference)
    reference_counts = Counter(g.family for g in reference.genes if g.family)
    orientation_rows: list[dict] = []

    # The reference is defined to retain its input direction; its contig order is still size-ranked.
    for contig in reference.contigs.values():
        contig.orientation = "forward"
        contig.orientation_status = "reference"
        contig.reference_contig = contig.name
        contig.orientation_anchors = sum(1 for gene in contig.genes if gene.family and reference_counts[gene.family] == 1)
        contig.kendall_tau_raw = 1.0 if contig.orientation_anchors >= 2 else None
        for gene in contig.genes:
            gene.orientation = "reference"
            gene.orientation_reference_contig = contig.name

    for genome in genomes:
        if genome.name == reference.name:
            continue
        genome_unique = _single_copy_families(genome)
        genome_family_counts = Counter(g.family for g in genome.genes if g.family)
        reference_by_contig: dict[str, dict[str, int]] = {}
        for ref_contig in reference.contigs.values():
            family_to_rank: dict[str, int] = {}
            for idx, gene in enumerate(ref_contig.genes, 1):
                if gene.family and reference_counts[gene.family] == 1:
                    family_to_rank[gene.family] = idx
            reference_by_contig[ref_contig.name] = family_to_rank

        for contig in genome.contigs.values():
            candidate_records = []
            query_genes = [g for g in contig.genes if g.family and genome_family_counts[g.family] == 1 and g.family in reference_unique]
            query_genes.sort(key=lambda g: (g.start, g.end, g.feature_id.casefold()))
            for ref_contig in reference.contigs.values():
                ref_ranks = reference_by_contig[ref_contig.name]
                anchors = [g for g in query_genes if g.family in ref_ranks]
                ranks = [ref_ranks[g.family] for g in anchors]
                tau = kendall_tau_from_ranks(ranks)
                candidate_records.append((len(anchors), abs(tau) if tau is not None else -1.0,
                                          ref_contig.name.casefold(), ref_contig, tau, anchors))
            if not candidate_records:
                contig.orientation = "undetermined"
                contig.orientation_status = "no_reference_contig_candidates"
                contig.reference_contig = ""
                contig.orientation_anchors = 0
                contig.kendall_tau_raw = None
            else:
                # Prefer the most informative reference contig, then the strongest monotonic relationship.
                candidate_records.sort(key=lambda x: (-x[0], -x[1], x[2]))
                n_anchors, _abs_tau, _name, ref_contig, tau, anchors = candidate_records[0]
                contig.reference_contig = ref_contig.name
                contig.orientation_anchors = n_anchors
                contig.kendall_tau_raw = tau
                tied_candidates = [candidate for candidate in candidate_records[1:]
                                   if candidate[0] == n_anchors and abs(candidate[1] - _abs_tau) < 1e-12]
                if tied_candidates:
                    contig.reference_contig = ""
                    contig.orientation = "undetermined"
                    contig.orientation_status = "ambiguous_reference_contig"
                elif n_anchors < min_anchors or tau is None or abs(tau) < min_abs_tau:
                    contig.orientation = "undetermined"
                    contig.orientation_status = "insufficient_or_ambiguous_anchors"
                elif tau >= 0:
                    contig.orientation = "forward"
                    contig.orientation_status = "inferred_from_family_order"
                else:
                    contig.orientation = "reverse"
                    contig.orientation_status = "inferred_from_family_order"
            for gene in contig.genes:
                gene.orientation = contig.orientation
                gene.orientation_reference_contig = contig.reference_contig
            orientation_rows.append({
                "genome": genome.name, "contig": contig.name, "contig_length": contig.length,
                "contig_rank_by_length": "", "reference_genome": reference.name,
                "reference_contig": contig.reference_contig, "orientation": contig.orientation,
                "orientation_status": contig.orientation_status, "unique_shared_families": contig.orientation_anchors,
                "kendall_tau_raw": contig.kendall_tau_raw,
            })

    # Keep the reference contigs in the orientation report too.
    for contig in reference.contigs.values():
        orientation_rows.append({
            "genome": reference.name, "contig": contig.name, "contig_length": contig.length,
            "contig_rank_by_length": "", "reference_genome": reference.name,
            "reference_contig": contig.name, "orientation": "forward",
            "orientation_status": "reference", "unique_shared_families": contig.orientation_anchors,
            "kendall_tau_raw": contig.kendall_tau_raw,
        })
    return reference.name, orientation_rows


def _contig_plot_gap(genome: Genome) -> int:
    n_contigs = len(genome.contigs)
    if n_contigs <= 1:
        return 0
    total = sum(contig.length for contig in genome.contigs.values())
    return max(100, min(5000, int(total * 0.002 / max(1, n_contigs - 1))))


def prepare_ordered_genomes(genomes: list[Genome], label_by: str) -> None:
    for genome in genomes:
        sorted_contigs = sorted(genome.contigs.values(), key=lambda c: (-c.length, c.name.casefold(), c.input_order))
        offset = 0
        gap = _contig_plot_gap(genome)
        genome_order = 0
        for contig_rank, contig in enumerate(sorted_contigs, 1):
            contig.genes.sort(key=lambda g: (g.start, g.end, g.feature_id.casefold()))
            previous_end = None
            for gene_order, gene in enumerate(contig.genes, 1):
                gene.gene_order_on_contig = gene_order
                gene_order += 1
                genome_order += 1
                gene.gene_order_genome = genome_order
                gene.normalized_strand = gene.strand
                if contig.orientation == "reverse":
                    gene.normalized_strand = {"+": "-", "-": "+"}.get(gene.strand, gene.strand)
                gene.normalized_gene_order = (len(contig.genes) - gene.gene_order_on_contig + 1
                                              if contig.orientation == "reverse" else gene.gene_order_on_contig)
                gene.bp_from_previous_gene = None if previous_end is None else gene.start - previous_end - 1
                previous_end = max(previous_end or 0, gene.end)
                if contig.orientation == "reverse":
                    gene.cumulative_start = offset + (contig.length - gene.end + 1)
                    gene.cumulative_end = offset + (contig.length - gene.start + 1)
                else:
                    gene.cumulative_start = offset + gene.start
                    gene.cumulative_end = offset + gene.end
                gene.orientation = contig.orientation
            # Leave a scale-aware gap between contigs, but no trailing gap after the last contig.
            offset += contig.length
            if contig_rank < len(sorted_contigs):
                offset += gap
        set_labels([genome], label_by)


def _contig_n50(lengths: list[int]) -> tuple[int, int]:
    lengths = sorted((x for x in lengths if x > 0), reverse=True)
    if not lengths:
        return 0, 0
    half = sum(lengths) / 2
    cumulative = 0
    for idx, length in enumerate(lengths, 1):
        cumulative += length
        if cumulative >= half:
            return length, idx
    return lengths[-1], len(lengths)


def build_breakpoint_reports(genomes: list[Genome], reference: Genome | None, family_map: FamilyMap | None,
                             window_genes: int = 10, max_end_anchor_genes: int = 5,
                             max_end_anchor_bp: int = 10000) -> tuple[list[dict], list[dict]]:
    """Map oriented query-contig ends onto a reference gene axis and summarise recurrent regions.

    Only contigs with an inferred forward/reverse orientation and single-copy family anchors
    are used. Each observed contig contributes a left and a right endpoint observation.
    """
    if reference is None or family_map is None or not any(g.family for x in genomes for g in x.genes):
        return [], []
    reference_counts = Counter(g.family for g in reference.genes if g.family)
    reference_gene_counts: dict[str, int] = {c.name: len(c.genes) for c in reference.contigs.values()}
    reference_family_ranks: dict[str, dict[str, int]] = {}
    for ref_contig in reference.contigs.values():
        rank_map: dict[str, int] = {}
        for rank, gene in enumerate(sorted(ref_contig.genes, key=lambda g: (g.start, g.end, g.feature_id.casefold())), 1):
            if gene.family and reference_counts[gene.family] == 1:
                rank_map[gene.family] = rank
        reference_family_ranks[ref_contig.name] = rank_map

    breakpoint_rows: list[dict] = []
    endpoint_rows: list[dict] = []
    for genome in genomes:
        if genome.name == reference.name:
            continue
        query_counts = Counter(g.family for g in genome.genes if g.family)
        for contig in genome.contigs.values():
            if contig.orientation not in {"forward", "reverse"} or not contig.reference_contig:
                continue
            ref_contig = reference.contigs.get(contig.reference_contig)
            if ref_contig is None:
                continue
            ref_ranks = reference_family_ranks.get(ref_contig.name, {})
            ordered_genes = sorted(contig.genes, key=lambda g: (g.start, g.end, g.feature_id.casefold()))
            if contig.orientation == "reverse":
                ordered_genes = ordered_genes[::-1]
            anchors = [(idx, gene, ref_ranks[gene.family]) for idx, gene in enumerate(ordered_genes)
                       if gene.family and query_counts[gene.family] == 1 and gene.family in ref_ranks]
            if len(anchors) < 2:
                continue
            left_idx, left_gene, left_rank = anchors[0]
            right_idx, right_gene, right_rank = anchors[-1]
            rank_values = [x[2] for x in anchors]
            ref_gene_count = reference_gene_counts.get(ref_contig.name, len(ref_contig.genes))
            row = {
                "genome": genome.name, "contig": contig.name, "contig_length_bp": contig.length,
                "orientation": contig.orientation, "reference_genome": reference.name,
                "reference_contig": ref_contig.name, "unique_shared_families": len(anchors),
                "reference_gene_count": ref_gene_count,
                "reference_order_min_anchor": min(rank_values), "reference_order_max_anchor": max(rank_values),
                "left_boundary_anchor_family": left_gene.family, "left_boundary_anchor_gene_id": left_gene.feature_id,
                "left_boundary_reference_gene_order": left_rank,
                "genes_from_left_contig_end_to_anchor": left_idx,
                "bp_from_left_contig_end_to_anchor": (left_gene.start - 1 if contig.orientation == "forward"
                                                       else max(0, contig.length - left_gene.end)),
                "left_boundary_support_eligible": (contig.length_source != "max_feature_end" and
                    left_idx + 1 <= max_end_anchor_genes and
                    (left_gene.start - 1 if contig.orientation == "forward" else max(0, contig.length - left_gene.end)) <= max_end_anchor_bp),
                "right_boundary_anchor_family": right_gene.family, "right_boundary_anchor_gene_id": right_gene.feature_id,
                "right_boundary_reference_gene_order": right_rank,
                "genes_from_right_contig_end_to_anchor": len(ordered_genes) - 1 - right_idx,
                "bp_from_right_contig_end_to_anchor": (max(0, contig.length - right_gene.end) if contig.orientation == "forward"
                                                        else right_gene.start - 1),
                "right_boundary_support_eligible": (contig.length_source != "max_feature_end" and
                    (len(ordered_genes) - right_idx) <= max_end_anchor_genes and
                    (max(0, contig.length - right_gene.end) if contig.orientation == "forward" else right_gene.start - 1) <= max_end_anchor_bp),
                "distance_left_anchor_to_reference_start_genes": max(0, left_rank - 1),
                "distance_right_anchor_to_reference_end_genes": max(0, ref_gene_count - right_rank),
                "near_reference_start_within_window": left_rank <= window_genes,
                "near_reference_end_within_window": (ref_gene_count - right_rank) < window_genes,
            }
            breakpoint_rows.append(row)
            endpoint_rows.append({"genome": genome.name, "contig": contig.name, "reference_contig": ref_contig.name,
                                  "boundary": "left", "reference_gene_order": left_rank,
                                  "family": left_gene.family, "anchor_gene_id": left_gene.feature_id,
                                  "genes_from_contig_end_to_anchor": left_idx,
                                  "bp_from_contig_end_to_anchor": (left_gene.start - 1 if contig.orientation == "forward"
                                                                    else max(0, contig.length - left_gene.end)),
                                  "eligible": (contig.length_source != "max_feature_end" and left_idx + 1 <= max_end_anchor_genes and
                                      (left_gene.start - 1 if contig.orientation == "forward" else max(0, contig.length - left_gene.end)) <= max_end_anchor_bp),
                                  "near_reference_contig_end": left_rank <= window_genes})
            endpoint_rows.append({"genome": genome.name, "contig": contig.name, "reference_contig": ref_contig.name,
                                  "boundary": "right", "reference_gene_order": right_rank,
                                  "family": right_gene.family, "anchor_gene_id": right_gene.feature_id,
                                  "genes_from_contig_end_to_anchor": len(ordered_genes) - 1 - right_idx,
                                  "bp_from_contig_end_to_anchor": (max(0, contig.length - right_gene.end) if contig.orientation == "forward"
                                                                    else right_gene.start - 1),
                                  "eligible": (contig.length_source != "max_feature_end" and
                                      (len(ordered_genes) - right_idx) <= max_end_anchor_genes and
                                      (max(0, contig.length - right_gene.end) if contig.orientation == "forward" else right_gene.start - 1) <= max_end_anchor_bp),
                                  "near_reference_contig_end": (ref_gene_count - right_rank) < window_genes})

    # Cluster nearby contig-end anchors into local reference-gene windows. Complete-linkage
    # by bounded rank span prevents a long chain of endpoints becoming one huge interval.
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in endpoint_rows:
        # Only anchors sufficiently close to a physical contig edge contribute to
        # recurrence counts. Both edge types are clustered together because a
        # shared breakpoint can be a left edge in one assembly and a right edge in another.
        if row["eligible"]:
            grouped[row["reference_contig"]].append(row)
    support_rows: list[dict] = []
    for ref_contig, rows in sorted(grouped.items()):
        rows.sort(key=lambda r: (r["reference_gene_order"], r["genome"].casefold(), r["contig"].casefold()))
        clusters: list[list[dict]] = []
        for row in rows:
            if not clusters or row["reference_gene_order"] - clusters[-1][0]["reference_gene_order"] > window_genes:
                clusters.append([row])
            else:
                clusters[-1].append(row)
        for cluster in clusters:
            ranks = [r["reference_gene_order"] for r in cluster]
            genomes_in_cluster = sorted({r["genome"] for r in cluster}, key=str.casefold)
            support_rows.append({
                "reference_genome": reference.name, "reference_contig": ref_contig,
                "boundary_sides": ",".join(sorted({r["boundary"] for r in cluster})),
                "reference_gene_order_min": min(ranks), "reference_gene_order_max": max(ranks),
                "reference_gene_order_median": statistics.median(ranks),
                "n_contig_end_observations": len(cluster), "n_distinct_genomes": len(genomes_in_cluster),
                "genomes": ",".join(genomes_in_cluster),
                "candidate_recurrent_region": len(genomes_in_cluster) >= 2,
                "near_reference_contig_end": any(r["near_reference_contig_end"] for r in cluster),
                "interpretation": "Candidate recurrent reference-projected assembly-end region; descriptive screen, not a significance test or biological breakpoint.",
            })
    return breakpoint_rows, support_rows


def build_physical_endpoint_rows(genomes: list[Genome]) -> list[dict]:
    """Emit exactly two FASTA-coordinate endpoints for every non-empty contig record."""
    rows = []
    for genome in genomes:
        for contig in genome.contigs.values():
            if contig.length <= 0 or not contig.has_fasta_record:
                continue
            for side, coordinate in (("left", 1), ("right", contig.length)):
                rows.append({"genome": genome.name, "contig": contig.name,
                    "contig_length": contig.length, "endpoint_side": side,
                    "physical_coordinate": coordinate, "replicon_id": contig.replicon_id,
                    "replicon_name": contig.replicon_name, "replicon_type": contig.replicon_type,
                    "replicon_topology": contig.replicon_topology,
                    "metadata_source": contig.metadata_source, "topology_source": contig.topology_source,
                    "classification_status": contig.classification_status,
                    "reference": "", "mapped_coordinate": "",
                    "mapping_method": "none", "mapping_confidence": 0,
                    "anchor_count": 0, "mapping_status": "unplaced"})
    return rows


def _write_tsv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        if fieldnames:
            writer.writeheader()
            writer.writerows(rows)


def _read_mapping_tsv(path: Path, required: tuple[str, ...]) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Mapping file does not exist: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not set(required).issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain TSV columns: {', '.join(required)}")
        return [{key: (row.get(key) or "").strip() for key in reader.fieldnames}
                for row in reader]


def build_reference_reports(genomes: list[Genome], replicon_metadata: list[dict], outdir: Path,
                            fasta_dir: Path | None, subcohort_map: Path | None = None,
                            reference_map: Path | None = None, assembly_metadata: Path | None = None,
                            similarity_matrix: Path | None = None,
                            legacy_reference: str | None = None, write_outputs: bool = True) -> dict:
    """Validate explicit sub-cohort/reference maps and write transparent candidate reports."""
    by_genome = {_normalise_genome_name(g.name): g for g in genomes}
    assignments = {key: "default" for key in by_genome}
    if subcohort_map:
        assignments = {}
        for line, row in enumerate(_read_mapping_tsv(subcohort_map, ("genome", "subcohort")), 2):
            genome, cohort = row["genome"], row["subcohort"]
            key = _normalise_genome_name(genome)
            if not genome or not cohort:
                raise ValueError(f"Sub-cohort map row {line} has an empty genome or subcohort")
            if key not in by_genome:
                raise ValueError(f"Sub-cohort map names unknown genome {genome!r} at row {line}")
            if key in assignments:
                raise ValueError(f"Duplicate sub-cohort assignment for genome {genome!r}")
            assignments[key] = cohort
        missing = [g.name for key, g in by_genome.items() if key not in assignments]
        if missing:
            raise ValueError("Sub-cohort map is missing genome assignments: " + ", ".join(missing))

    cohort_names = sorted(set(assignments.values()), key=str.casefold)
    selected: dict[str, str] = {}
    selected_presets: dict[str, str] = {}
    if legacy_reference and reference_map:
        raise ValueError("Use either --reference-genome or --reference-map, not both")
    if legacy_reference:
        if len(cohort_names) != 1:
            raise ValueError("--reference-genome cannot define references for multiple subcohorts; use --reference-map")
        key = _normalise_genome_name(legacy_reference)
        if key not in by_genome:
            raise ValueError(f"Reference genome {legacy_reference!r} is absent from the input cohort")
        selected[cohort_names[0]] = by_genome[key].name
    if reference_map:
        for line, row in enumerate(_read_mapping_tsv(reference_map, ("subcohort", "reference_genome")), 2):
            cohort, reference = row["subcohort"], row["reference_genome"]
            if not cohort or not reference:
                raise ValueError(f"Reference map row {line} has an empty subcohort or reference_genome")
            if cohort not in cohort_names:
                raise ValueError(f"Reference map names unknown subcohort {cohort!r} at row {line}")
            if cohort in selected:
                raise ValueError(f"Duplicate reference assignment for subcohort {cohort!r}")
            key = _normalise_genome_name(reference)
            if key not in by_genome:
                raise ValueError(f"Reference genome {reference!r} is absent from the input cohort")
            if assignments[key] != cohort:
                raise ValueError(f"Reference genome {reference!r} is assigned to {assignments[key]!r}, not {cohort!r}")
            if not fasta_dir:
                raise ValueError(f"Reference genome {reference!r} requires --fasta-dir")
            g = by_genome[key]
            fasta_matches = [p for p in fasta_dir.rglob("*") if p.is_file() and p.name.lower().endswith(FASTA_SUFFIXES)
                             and _normalise_genome_name(_stem_without_suffix(p, FASTA_SUFFIXES)) == key]
            if len(fasta_matches) != 1:
                raise ValueError(f"Reference genome {reference!r} must have exactly one matching FASTA assembly")
            if not g.genes or not g.selected_feature_type:
                raise ValueError(f"Reference genome {reference!r} has no usable GFF annotation for gene-family landmarks")
            rows = [r for r in replicon_metadata if _normalise_genome_name(r["genome"]) == key]
            if not rows or len({r["contig"] for r in rows}) != len(rows):
                raise ValueError(f"Reference genome {reference!r} has missing or inconsistent FASTA sequence metadata")
            conflicts = [r["contig"] for r in rows if r.get("conflict_status", "none") != "none"]
            if conflicts:
                raise ValueError(f"Reference genome {reference!r} has conflicting replicon metadata for: {', '.join(conflicts)}")
            selected[cohort] = g.name
            cohort_preset = row.get("minimap2_preset", "")
            if cohort_preset:
                if cohort_preset not in {"asm5", "asm10", "asm20"}:
                    raise ValueError(
                        f"Unsupported minimap2_preset {cohort_preset!r} for subcohort {cohort!r}; "
                        "choose asm5, asm10, or asm20"
                    )
                selected_presets[cohort] = cohort_preset
        missing = sorted(set(cohort_names) - set(selected), key=str.casefold)
        if missing:
            raise ValueError("Reference map is missing sub-cohort assignments: " + ", ".join(missing))

    extra_metadata: dict[str, dict[str, str]] = {}
    if assembly_metadata:
        for line, row in enumerate(_read_mapping_tsv(assembly_metadata, ("genome",)), 2):
            key = _normalise_genome_name(row["genome"])
            if key not in by_genome:
                raise ValueError(f"Assembly metadata names unknown genome {row['genome']!r} at row {line}")
            if key in extra_metadata:
                raise ValueError(f"Duplicate assembly metadata for genome {row['genome']!r}")
            extra_metadata[key] = row

    similarity: dict[tuple[str, str], float] = {}
    if similarity_matrix:
        with similarity_matrix.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.reader(handle))
        if len(rows) < 2 or len(rows[0]) < 2:
            raise ValueError("Similarity matrix must be a labeled square CSV with a header row and one row per genome")
        labeled_rows = rows[0][0].strip() == ""
        if labeled_rows:
            labels = [x.strip() for x in rows[0][1:]]
            data_rows = rows[1:]
            if len(data_rows) != len(labels):
                raise ValueError("Similarity matrix must be square")
        else:
            labels = [x.strip() for x in rows[0]]
            data_rows = rows[1:]
            if len(data_rows) != len(labels):
                raise ValueError("Similarity matrix must be square")
        if len(set(map(_normalise_genome_name, labels))) != len(labels):
            raise ValueError("Similarity matrix has duplicate genome labels")
        for index, row in enumerate(data_rows):
            expected = len(labels) + 1 if labeled_rows else len(labels)
            if len(row) != expected:
                raise ValueError("Similarity matrix rows must match the header width")
            row_label = row[0] if labeled_rows else labels[index]
            if labeled_rows and _normalise_genome_name(row_label) != _normalise_genome_name(labels[index]):
                raise ValueError("Similarity matrix row and column labels are not in matching order")
            values = row[1:] if labeled_rows else row
            for label, value in zip(labels, values):
                try:
                    similarity[(_normalise_genome_name(row_label), _normalise_genome_name(label))] = float(value)
                except ValueError as exc:
                    raise ValueError(f"Non-numeric similarity in {similarity_matrix}") from exc

    fasta_paths = {}
    if fasta_dir:
        for path in fasta_dir.rglob("*"):
            if path.is_file() and path.name.lower().endswith(FASTA_SUFFIXES):
                fasta_paths.setdefault(_normalise_genome_name(_stem_without_suffix(path, FASTA_SUFFIXES)), []).append(path)
    metadata_by_genome: dict[str, list[dict]] = defaultdict(list)
    for row in replicon_metadata:
        metadata_by_genome[_normalise_genome_name(row["genome"])].append(row)

    candidate_rows = []
    stats = {}
    for key, genome in by_genome.items():
        fasta_contigs = [c for c in genome.contigs.values() if c.has_fasta_record]
        lengths = sorted((c.length for c in fasta_contigs), reverse=True)
        n50, _l50 = _contig_n50(lengths)
        total = sum(lengths)
        rep_rows = metadata_by_genome.get(key, [])
        known = sum(all(r.get(field) not in (None, "", "unknown") for field in
                        ("replicon_id", "replicon_type", "replicon_topology")) for r in rep_rows)
        types = {r.get("replicon_type") for r in rep_rows}
        cohort = assignments[key]
        peers = [other for other, group in assignments.items() if group == cohort and other != key]
        vals = [similarity[(key, peer)] for peer in peers if (key, peer) in similarity]
        warnings = []
        if not fasta_contigs: warnings.append("FASTA assembly unavailable")
        if not genome.genes: warnings.append("GFF annotation unavailable")
        if not rep_rows: warnings.append("replicon metadata unavailable")
        if extra_metadata.get(key, {}).get("assembly_level") in (None, ""):
            warnings.append("assembly level not supplied")
        if extra_metadata.get(key, {}).get("completeness") in (None, ""):
            warnings.append("completeness metadata not supplied")
        stats[key] = {"length": total, "n50": n50, "ncontigs": len(fasta_contigs), "known_replicons": known,
                      "metadata_n": len(rep_rows), "annotations": len(genome.genes), "conflicts": sum(r.get("conflict_status", "none") != "none" for r in rep_rows)}
        candidate_rows.append({
            "genome": genome.name, "subcohort": cohort, "n_fasta_contigs": len(fasta_contigs),
            "n_gff_annotated_contigs": sum(c.has_gff_annotation and c.has_fasta_record for c in genome.contigs.values()),
            "assembly_length_bp": total, "n50_bp": n50,
            "assembly_level": extra_metadata.get(key, {}).get("assembly_level", ""),
            "completeness": extra_metadata.get(key, {}).get("completeness", ""),
            "has_chromosome_annotation": "yes" if "chromosome" in types else "no" if rep_rows else "unknown",
            "has_plasmid_annotation": "yes" if "plasmid" in types else "no" if rep_rows else "unknown",
            "replicon_metadata_completeness": f"{known}/{len(rep_rows)}" if rep_rows else "unknown",
            "replicon_metadata_conflict_count": stats[key]["conflicts"],
            "sourmash_within_subcohort_n": len(vals), "sourmash_mean_similarity": statistics.mean(vals) if vals else "",
            "sourmash_min_similarity": min(vals) if vals else "", "sourmash_max_similarity": max(vals) if vals else "",
            "reference_suitability_warnings": "; ".join(warnings),
        })
    # A transparent lexicographic review order; this is not a completeness score and never selects a reference.
    for cohort in cohort_names:
        group = [r for r in candidate_rows if r["subcohort"] == cohort]
        group.sort(key=lambda r: (stats[_normalise_genome_name(r["genome"])]["conflicts"],
                                  -stats[_normalise_genome_name(r["genome"])]["known_replicons"],
                                  -r["n_gff_annotated_contigs"], -r["n50_bp"], r["n_fasta_contigs"], r["genome"].casefold()))
        for rank, row in enumerate(group, 1):
            row["candidate_review_rank"] = rank
            row["recommendation"] = "candidate for manual review; not selected"
            row["candidate_review_basis"] = "conflict count ascending; known replicon records descending; GFF-annotated contigs descending; N50 descending; FASTA contig count ascending; genome ID ascending"
    candidate_rows.sort(key=lambda r: (r["subcohort"].casefold(), r["candidate_review_rank"]))
    if write_outputs:
        _write_tsv(outdir / "reference_candidates.tsv", candidate_rows)

    manifest_rows = []
    for cohort in cohort_names:
        ref_name = selected.get(cohort, "")
        key = _normalise_genome_name(ref_name) if ref_name else ""
        genome = by_genome.get(key) if key else None
        rep_rows = metadata_by_genome.get(key, []) if key else []
        replicons = {}
        for row in rep_rows:
            rid = row.get("replicon_id", "unknown")
            group_id = rid if rid not in ("", "unknown") else f"unknown:{row['contig']}"
            replicons.setdefault(group_id, {"replicon_id": rid, "replicon_type": row.get("replicon_type", "unknown"),
                                       "replicon_topology": row.get("replicon_topology", "unknown"), "contigs": []})["contigs"].append(row["contig"])
        warnings = [] if ref_name else ["No reference was specified; candidate ranking is advisory and no genome was selected."]
        if genome and not genome.genes: warnings.append("No GFF annotation available for family landmarks.")
        assembly = stats[key] if key else {}
        fasta = fasta_paths.get(key, []) if key else []
        if ref_name and not fasta:
            warnings.append("Reference FASTA path is unavailable; selection metadata only.")
        manifest_rows.append({"subcohort": cohort, "selected_reference_genome": ref_name,
            "selection_method": "user-specified" if ref_name else "not_selected",
            "reference_fasta_path": str(fasta[0]) if len(fasta) == 1 else "",
            "reference_gff_path": str(genome.path) if genome else "",
            "assembly_length_bp": assembly.get("length", ""), "n50_bp": assembly.get("n50", ""),
            "n_contigs": assembly.get("ncontigs", ""), "reference_replicons": json.dumps(list(replicons.values()), sort_keys=True),
            "available_annotations": genome.selected_feature_type if genome and genome.genes else "",
            "selection_warnings": "; ".join(warnings)})
    if write_outputs:
        _write_tsv(outdir / "reference_manifest.tsv", manifest_rows)
    return {"candidate_count": len(candidate_rows), "subcohort_count": len(cohort_names),
            "selected_references": selected, "assignments": assignments,
            "presets": selected_presets}


def build_contig_end_gene_candidates(genomes: list[Genome], max_end_anchor_genes: int = 5,
                                    max_end_anchor_bp: int = 10000) -> list[dict]:
    """Report annotated genes near physical contig ends as review candidates.

    Proximity and annotation flags do not establish truncation; they identify
    loci that deserve sequence-level inspection.
    """
    rows = []
    for genome in genomes:
        for contig in genome.contigs.values():
            ordered = sorted(contig.genes, key=lambda gene: (gene.start, gene.end, gene.feature_id.casefold()))
            count = len(ordered)
            for index, gene in enumerate(ordered):
                left_bp = max(0, gene.start - 1)
                right_bp = max(0, contig.length - gene.end)
                near_left = index + 1 <= max_end_anchor_genes and left_bp <= max_end_anchor_bp
                near_right = (contig.length_source != "max_feature_end" and
                              count - index <= max_end_anchor_genes and right_bp <= max_end_anchor_bp)
                flags = []
                for key in ("partial", "pseudo", "pseudogene", "start_range", "end_range", "exception",
                            "frameshift", "ribosomal_slippage"):
                    value = gene.attributes.get(key, "")
                    normalized = value.casefold().strip()
                    # Pyrodigal encodes the two partial-boundary states as
                    # binary digits (00 = complete, 10/01/11 = partial).
                    if key == "partial" and normalized.isdigit():
                        if any(digit != "0" for digit in normalized):
                            flags.append(f"{key}={value}")
                        continue
                    if normalized and normalized not in {"false", "0", "no", "none", "00"}:
                        flags.append(f"{key}={value}")
                if not (near_left or near_right or flags):
                    continue
                rows.append({
                    "genome": genome.name, "contig": contig.name, "contig_length_bp": contig.length,
                    "length_source": contig.length_source,
                    "gene_id": gene.feature_id, "family": gene.family, "product": gene.product,
                    "start": gene.start, "end": gene.end, "strand": gene.strand, "gene_length_bp": gene.length,
                    "genes_from_left_contig_end": index,
                    "bp_from_left_contig_end": left_bp,
                    "genes_from_right_contig_end": count - index - 1,
                    "bp_from_right_contig_end": right_bp,
                    "near_left_end": near_left, "near_right_end": near_right,
                    "annotation_flags": ";".join(flags),
                    "candidate_reason": ";".join(
                        reason for reason, present in (("near_contig_end", near_left or near_right),
                                                       ("partial_or_disrupted_annotation_flag", bool(flags))) if present
                    ),
                    "interpretation": "Candidate only; end proximity or annotation flags do not prove truncation or gene loss.",
                })
    return rows


def build_reports(genomes: list[Genome], outdir: Path, family_map: FamilyMap | None,
                  reference_genome: str | None, orientation_rows: list[dict], breakpoint_window_genes: int = 10,
                  max_end_anchor_genes: int = 5, max_end_anchor_bp: int = 10000,
                  replicon_metadata: list[dict] | None = None) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    genome_rows: list[dict] = []
    contig_rows: list[dict] = []
    gene_rows: list[dict] = []
    contig_order_rows: list[dict] = []
    comparison_rows: list[dict] = []
    order_metric_rows: list[dict] = []
    reference = next((g for g in genomes if reference_genome and g.name.casefold() == reference_genome.casefold()), None)
    replicon_lookup = {(row["genome"], row["contig"]): row for row in (replicon_metadata or [])}
    ref_counts = Counter(g.family for g in reference.genes if g.family) if reference else Counter()
    ref_unique = _single_copy_families(reference) if reference else {}

    for genome in genomes:
        contigs = sorted(genome.contigs.values(), key=lambda c: (-c.length, c.name.casefold(), c.input_order))
        lengths = [c.length for c in contigs]
        n50, l50 = _contig_n50(lengths)
        genome_length = sum(lengths)
        families = [g.family for g in genome.genes if g.family]
        family_counts = Counter(families)
        single_copy = sum(1 for fam in family_counts.values() if fam == 1)
        multi_copy = sum(1 for fam in family_counts.values() if fam > 1)
        main_contig_genes = len(contigs[0].genes) if contigs else 0
        gene_density = (sum(g.length for g in genome.genes) / genome_length) if genome_length else 0.0
        genome_rows.append({
            "genome": genome.name, "gff": str(genome.path), "selected_feature_type": genome.selected_feature_type,
            "n_contigs": len(contigs), "n_fasta_contigs": sum(c.has_fasta_record for c in contigs),
            "n_gff_annotated_contigs": sum(c.has_gff_annotation for c in contigs),
            "n_gff_only_contigs": sum(c.has_gff_annotation and not c.has_fasta_record for c in contigs),
            "n_genes": len(genome.genes), "total_contig_length_bp": genome_length,
            "length_basis": ";".join(sorted({c.length_source for c in contigs})),
            "longest_contig_bp": max(lengths, default=0), "n50_bp": n50, "l50_contigs": l50,
            "median_contig_length_bp": statistics.median(lengths) if lengths else 0,
            "genes_on_longest_contig": main_contig_genes,
            "fraction_genes_on_longest_contig": main_contig_genes / len(genome.genes) if genome.genes else 0,
            "summed_gene_span_bp": sum(g.length for g in genome.genes),
            "gene_span_fraction_of_assembly": gene_density,
            "family_mapped_genes": sum(1 for g in genome.genes if g.family),
            "unique_families": single_copy, "multicopy_families": multi_copy,
            "n_gff_warnings": len(genome.warnings),
        })
        offset = 0
        cumulative_contig_length = 0
        genome_gene_idx = 0
        for rank, contig in enumerate(contigs, 1):
            genes = sorted(contig.genes, key=lambda g: (g.start, g.end, g.feature_id.casefold()))
            strand_counts = Counter(g.strand for g in genes)
            contig_rows.append({
                "genome": genome.name, "contig": contig.name, "contig_rank_by_length": rank,
                "input_contig_order": contig.input_order, "contig_length_bp": contig.length,
                "length_source": contig.length_source, "n_genes": len(genes),
                "has_fasta_record": contig.has_fasta_record,
                "has_gff_annotation": contig.has_gff_annotation,
                "first_gene_start": min((g.start for g in genes), default=""),
                "last_gene_end": max((g.end for g in genes), default=""),
                "plus_strand_genes": strand_counts.get("+", 0), "minus_strand_genes": strand_counts.get("-", 0),
                "unknown_strand_genes": sum(strand_counts.get(x, 0) for x in (".", "?")),
                "family_mapped_genes": sum(1 for g in genes if g.family),
                "orientation": contig.orientation, "orientation_status": contig.orientation_status,
                "reference_genome": reference.name if reference else "",
                "reference_contig": contig.reference_contig, "unique_shared_families": contig.orientation_anchors,
                "kendall_tau_raw": contig.kendall_tau_raw if contig.kendall_tau_raw is not None else "",
            })
            contig_order_rows.append({
                "genome": genome.name, "contig_rank_by_length": rank, "contig": contig.name,
                "contig_length_bp": contig.length, "length_source": contig.length_source,
                "fraction_of_genome_length": contig.length / genome_length if genome_length else 0.0,
                "cumulative_contig_length_bp": cumulative_contig_length + contig.length, "n_genes": len(genes),
            })
            cumulative_contig_length += contig.length
            for gene in genes:
                genome_gene_idx += 1
                gene_rows.append({
                    "genome": genome.name, "contig": contig.name, "contig_rank_by_length": rank,
                    "gene_order_on_contig": gene.gene_order_on_contig, "gene_order_genome_by_contig_rank": genome_gene_idx,
                    "gene_id": gene.feature_id, "locus_tag": gene.locus_tag, "gene_name": gene.gene_name,
                    "product": gene.product, "family": gene.family, "label": gene.label,
                    "feature_type": gene.feature_type, "start": gene.start, "end": gene.end,
                    "midpoint": gene.midpoint, "length_bp": gene.length, "strand": gene.strand,
                    "normalized_strand": gene.normalized_strand, "normalized_gene_order": gene.normalized_gene_order,
                    "bp_from_previous_gene": gene.bp_from_previous_gene,
                    "contig_length_bp": contig.length, "length_source": contig.length_source,
                    "orientation": contig.orientation, "orientation_reference_contig": contig.reference_contig,
                    "cumulative_plot_start": gene.cumulative_start, "cumulative_plot_end": gene.cumulative_end,
                })

    _write_tsv(outdir / "genome_summary.tsv", genome_rows)
    _write_tsv(outdir / "contig_summary.tsv", contig_rows)
    _write_tsv(outdir / "contig_order.tsv", contig_order_rows)
    _write_tsv(outdir / "gene_order.tsv", gene_rows)
    _write_tsv(outdir / "contig_orientation.tsv", orientation_rows)

    if reference and family_map and any(g.family for genome in genomes for g in genome.genes):
        for genome in genomes:
            if genome.name == reference.name:
                continue
            # Each query contig is compared only with the single best-matching reference contig.
            ranked_contigs = sorted(genome.contigs.values(), key=lambda c: (-c.length, c.name.casefold(), c.input_order))
            contig_rank_lookup = {c.name: i for i, c in enumerate(ranked_contigs, 1)}
            query_family_counts = Counter(g.family for g in genome.genes if g.family)
            for contig in genome.contigs.values():
                ref_contig_name = contig.reference_contig
                ref_contig = reference.contigs.get(ref_contig_name) if ref_contig_name else None
                if not ref_contig:
                    continue
                ref_rank_by_family = {}
                for idx, rg in enumerate(sorted(ref_contig.genes, key=lambda g: (g.start, g.end, g.feature_id.casefold())), 1):
                    if rg.family and ref_counts[rg.family] == 1:
                        ref_rank_by_family[rg.family] = idx
                qgenes = sorted(contig.genes, key=lambda g: (g.start, g.end, g.feature_id.casefold()))
                shared = [g for g in qgenes if g.family in ref_rank_by_family and query_family_counts[g.family] == 1
                          and family_map.lookup(reference.name, _gene_lookup_keys(ref_unique[g.family][1]))]
                shared_families = [g.family for g in shared]
                tau_raw = kendall_tau_from_ranks([ref_rank_by_family[g.family] for g in shared])
                oriented_shared = shared[::-1] if contig.orientation == "reverse" else shared
                oriented_ranks = [ref_rank_by_family[g.family] for g in oriented_shared]
                oriented_tau = kendall_tau_from_ranks(oriented_ranks)
                adjacent_pairs = list(zip(oriented_ranks, oriented_ranks[1:]))
                adjacent_anchor_concordance = (
                    sum(1 for left, right in adjacent_pairs if right > left) / len(adjacent_pairs)
                    if adjacent_pairs else ""
                )
                order_metric_rows.append({
                    "genome": genome.name, "contig": contig.name, "reference_genome": reference.name,
                    "reference_contig": ref_contig.name, "orientation": contig.orientation,
                    "shared_unique_families": len(shared), "kendall_tau_raw": tau_raw if tau_raw is not None else "",
                    "kendall_tau_after_orientation_normalisation": oriented_tau if oriented_tau is not None else "",
                    "consecutive_anchor_order_concordance": adjacent_anchor_concordance,
                    "metric_note": "Anchor-order metrics use families single-copy within each genome; contigs are compared separately.",
                })
                for g in qgenes:
                    if not g.family or g.family not in ref_rank_by_family or query_family_counts[g.family] != 1:
                        continue
                    ref_feature = ref_unique.get(g.family, (None, None))[1]
                    if ref_feature is None:
                        continue
                    comparison_rows.append({
                        "genome": genome.name, "contig": contig.name, "contig_rank_by_length": contig_rank_lookup.get(contig.name, ""),
                        "gene_id": g.feature_id, "family": g.family, "query_gene_order_on_contig": g.gene_order_on_contig,
                        "query_normalized_gene_order": g.normalized_gene_order,
                        "reference_genome": reference.name, "reference_contig": ref_contig.name,
                        "reference_gene_id": ref_feature.feature_id,
                        "reference_gene_order_on_contig": next((i for i, rg in enumerate(sorted(ref_contig.genes, key=lambda x: (x.start, x.end, x.feature_id.casefold())), 1) if rg is ref_feature), ""),
                        "orientation": contig.orientation,
                        "query_strand": g.strand, "query_normalized_strand": g.normalized_strand,
                        "query_start": g.start, "query_end": g.end,
                    })
        _write_tsv(outdir / "gene_order_metrics_vs_reference.tsv", order_metric_rows)
        _write_tsv(outdir / "family_order_matches.tsv", comparison_rows)
    else:
        _write_tsv(outdir / "gene_order_metrics_vs_reference.tsv", [], [
            "genome", "contig", "reference_genome", "reference_contig", "orientation", "shared_unique_families",
            "kendall_tau_raw", "kendall_tau_after_orientation_normalisation", "consecutive_anchor_order_concordance", "metric_note",
        ])
        _write_tsv(outdir / "family_order_matches.tsv", [], [
            "genome", "contig", "contig_rank_by_length", "gene_id", "family", "query_gene_order_on_contig",
            "query_normalized_gene_order", "reference_genome", "reference_contig", "reference_gene_id",
            "reference_gene_order_on_contig", "orientation", "query_strand", "query_normalized_strand", "query_start", "query_end",
        ])

    breakpoint_rows, breakpoint_support_rows = build_breakpoint_reports(
        genomes, reference, family_map if family_map and any(g.family for x in genomes for g in x.genes) else None,
        window_genes=breakpoint_window_genes, max_end_anchor_genes=max_end_anchor_genes,
        max_end_anchor_bp=max_end_anchor_bp,
    )
    physical_endpoint_rows = build_physical_endpoint_rows(genomes)
    for endpoint in physical_endpoint_rows:
        metadata = replicon_lookup.get((endpoint["genome"], endpoint["contig"]), {})
        for field in ("replicon_id", "replicon_name", "replicon_type", "replicon_topology",
                      "sequence_role", "classification_status", "conflict_status"):
            endpoint[field] = metadata.get(field, "unknown" if field in {"replicon_type", "replicon_topology"} else "")
        endpoint["metadata_sources"] = json.dumps(metadata.get("metadata_sources", []), sort_keys=True)
        endpoint["metadata_source"] = endpoint["metadata_sources"]
        endpoint["topology_source"] = json.dumps({
            source: values for source, values in metadata.get("metadata_evidence", {}).items()
            if source in {"fasta_header", "gff_region", "assembly_report"}
        }, sort_keys=True) if metadata.get("replicon_topology") != "unknown" else ""
        endpoint["metadata_evidence"] = json.dumps(metadata.get("metadata_evidence", {}), sort_keys=True)
        endpoint["resolution_notes"] = metadata.get("resolution_notes", "")
        endpoint["matched_accessions"] = json.dumps({
            "genbank": metadata.get("matched_genbank_accessions", []),
            "refseq": metadata.get("matched_refseq_accessions", []),
        }, sort_keys=True)
        endpoint["assembly_report_path"] = metadata.get("assembly_report_source", "")
        endpoint["has_gff_annotation"] = metadata.get("gff_annotation_status") == "annotated"
    _write_tsv(outdir / "physical_endpoints.tsv", physical_endpoint_rows, [
        "genome", "contig", "contig_length", "endpoint_side", "physical_coordinate", "replicon_id",
        "replicon_name", "replicon_type", "replicon_topology", "metadata_source", "topology_source",
        "sequence_role", "classification_status", "conflict_status", "metadata_sources", "metadata_evidence",
        "resolution_notes", "matched_accessions", "assembly_report_path", "has_gff_annotation",
        "reference", "mapped_coordinate", "mapping_method", "mapping_confidence", "anchor_count", "mapping_status",
    ])
    manifest_rows = []
    for row in replicon_metadata or []:
        manifest_rows.append({
            "genome": row["genome"], "contig": row["contig"], "contig_length": row["contig_length"],
            "has_gff_annotation": row["gff_annotation_status"] == "annotated",
            "replicon_id": row.get("replicon_id", ""), "replicon_name": row.get("replicon_name", ""),
            "replicon_type": row.get("replicon_type", "unknown"),
            "replicon_topology": row.get("replicon_topology", "unknown"),
            "sequence_role": row.get("sequence_role", ""),
            "classification_status": row.get("classification_status", "unknown"),
            "conflict_status": row.get("conflict_status", "none"),
            "metadata_sources": json.dumps(row.get("metadata_sources", []), sort_keys=True),
            "metadata_evidence": json.dumps(row.get("metadata_evidence", {}), sort_keys=True),
            "resolution_notes": row.get("resolution_notes", ""),
            "matched_accessions": json.dumps({"genbank": row.get("matched_genbank_accessions", []),
                                               "refseq": row.get("matched_refseq_accessions", [])}, sort_keys=True),
            "assembly_report_path": row.get("assembly_report_source", ""),
        })
    manifest_fields = ["genome", "contig", "contig_length", "has_gff_annotation", "replicon_id", "replicon_name",
        "replicon_type", "replicon_topology", "sequence_role", "classification_status", "conflict_status",
        "metadata_sources", "metadata_evidence", "resolution_notes", "matched_accessions", "assembly_report_path"]
    _write_tsv(outdir / "replicon_manifest.tsv", manifest_rows, manifest_fields)
    _write_tsv(outdir / "contig_breakpoints_vs_reference.tsv", breakpoint_rows, [
        "genome", "contig", "contig_length_bp", "orientation", "reference_genome", "reference_contig",
        "unique_shared_families", "reference_gene_count", "reference_order_min_anchor", "reference_order_max_anchor",
        "left_boundary_anchor_family", "left_boundary_anchor_gene_id", "left_boundary_reference_gene_order",
        "genes_from_left_contig_end_to_anchor", "bp_from_left_contig_end_to_anchor", "left_boundary_support_eligible",
        "right_boundary_anchor_family", "right_boundary_anchor_gene_id", "right_boundary_reference_gene_order",
        "genes_from_right_contig_end_to_anchor", "bp_from_right_contig_end_to_anchor", "right_boundary_support_eligible",
        "distance_left_anchor_to_reference_start_genes", "distance_right_anchor_to_reference_end_genes",
        "near_reference_start_within_window", "near_reference_end_within_window",
    ])
    _write_tsv(outdir / "assembly_breakpoint_support.tsv", breakpoint_support_rows, [
        "reference_genome", "reference_contig", "boundary_sides", "reference_gene_order_min", "reference_gene_order_max",
        "reference_gene_order_median", "n_contig_end_observations", "n_distinct_genomes", "genomes",
        "candidate_recurrent_region", "near_reference_contig_end", "interpretation",
    ])
    contig_end_candidates = build_contig_end_gene_candidates(
        genomes, max_end_anchor_genes=max_end_anchor_genes, max_end_anchor_bp=max_end_anchor_bp
    )
    _write_tsv(outdir / "contig_end_gene_candidates.tsv", contig_end_candidates, [
        "genome", "contig", "contig_length_bp", "length_source", "gene_id", "family", "product", "start", "end", "strand",
        "gene_length_bp", "genes_from_left_contig_end", "bp_from_left_contig_end",
        "genes_from_right_contig_end", "bp_from_right_contig_end", "near_left_end", "near_right_end",
        "annotation_flags", "candidate_reason", "interpretation",
    ])

    return {
        "n_genomes": len(genomes), "n_contigs": sum(len(g.contigs) for g in genomes),
        "n_genes": sum(len(g.genes) for g in genomes), "n_family_mapped_genes": sum(1 for g in genomes for x in g.genes if x.family),
        "reference_genome": reference.name if reference else "", "n_contig_orientation_records": len(orientation_rows),
        "n_family_order_matches": len(comparison_rows), "n_contig_breakpoints_mapped": len(breakpoint_rows),
        "n_recurrent_breakpoint_regions": sum(1 for r in breakpoint_support_rows if r["candidate_recurrent_region"]),
        "n_contig_end_gene_candidates": len(contig_end_candidates),
        "n_genomes_with_warnings": sum(1 for g in genomes if g.warnings),
    }


def _require_matplotlib():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
        from matplotlib.lines import Line2D
        return matplotlib, plt, PolyCollection, Line2D
    except ImportError as exc:
        raise RuntimeError("Plotting requires matplotlib. Install PyamilySeq's analysis extra or run with --no-plots.") from exc


def make_plots(genomes: list[Genome], outdir: Path, reference_genome: str | None,
               max_genomes_plotted: int = 50, max_labels_per_genome: int = 20) -> list[str]:
    _matplotlib, plt, PolyCollection, Line2D = _require_matplotlib()
    plot_dir = outdir / "figures"
    plot_dir.mkdir(parents=True, exist_ok=True)
    notes: list[str] = []
    displayed = genomes[:max_genomes_plotted]
    if len(genomes) > len(displayed):
        notes.append(f"Figures include the first {len(displayed)} genomes in alphabetical input order; tabular reports include all {len(genomes)}.")

    # Plot 1: contig size profile by descending contig rank.
    fig_height = max(4.5, min(14.0, 0.22 * len(displayed) + 2.2))
    fig, ax = plt.subplots(figsize=(10.5, fig_height))
    max_contigs = max((len(g.contigs) for g in displayed), default=1)
    for genome in displayed:
        contigs = sorted((c for c in genome.contigs.values() if c.length > 0), key=lambda c: (-c.length, c.name.casefold(), c.input_order))
        omitted = len(genome.contigs) - len(contigs)
        if omitted:
            notes.append(f"Contig size profile omits {omitted} zero-length contig(s) for {genome.name} because the y-axis is logarithmic.")
        ax.plot(range(1, len(contigs) + 1), [c.length / 1e6 for c in contigs], marker="o", markersize=2.5,
                linewidth=1.0, alpha=0.7, label=genome.name)
    ax.set_xlabel("Contig rank (longest first)")
    ax.set_ylabel("Contig length (Mb; log scale)")
    ax.set_yscale("log")
    ax.set_xlim(1, max(2, max_contigs))
    ax.grid(True, which="both", axis="y", linestyle=":", linewidth=0.5, alpha=0.5)
    if len(displayed) <= 20:
        ax.legend(loc="best", fontsize=7, frameon=False, ncol=1)
    else:
        ax.text(0.99, 0.98, f"{len(displayed)} genomes shown", ha="right", va="top", transform=ax.transAxes, fontsize=8)
    ax.set_title("Contig length profiles (contigs ranked within each genome)")
    fig.tight_layout()
    fig.savefig(plot_dir / "contig_size_profiles.png", dpi=220, bbox_inches="tight")
    fig.savefig(plot_dir / "contig_size_profiles.svg", bbox_inches="tight")
    plt.close(fig)

    # Plot 2: arrows for each gene on concatenated, length-ranked contig tracks.
    fig_height = max(4.5, min(18.0, 0.28 * len(displayed) + 2.0))
    fig, ax = plt.subplots(figsize=(15.0, fig_height))
    plus_color, minus_color, unknown_color = "#2979A8", "#E07A27", "#9A9A9A"
    all_x_max = 0
    for row_index, genome in enumerate(displayed):
        y = len(displayed) - row_index - 1
        contigs = sorted(genome.contigs.values(), key=lambda c: (-c.length, c.name.casefold(), c.input_order))
        offset = 0
        gap = _contig_plot_gap(genome)
        polys_plus, polys_minus, polys_unknown = [], [], []
        gene_label_candidates = []
        for c_rank, contig in enumerate(contigs, 1):
            genes = sorted(contig.genes, key=lambda g: (g.start, g.end, g.feature_id.casefold()))
            for gene in genes:
                if contig.orientation == "reverse":
                    x0 = offset + (contig.length - gene.end + 1)
                    x1 = offset + (contig.length - gene.start + 1)
                    plot_strand = {"+": "-", "-": "+"}.get(gene.strand, gene.strand)
                else:
                    x0 = offset + gene.start
                    x1 = offset + gene.end
                    plot_strand = gene.normalized_strand if contig.orientation == "forward" else gene.strand
                real_width = max(1, x1 - x0 + 1)
                draw_width = max(real_width, contig.length * 0.00018)
                if plot_strand == "+":
                    left, right = x0, x0 + draw_width
                    color_target = polys_plus
                    points = [(left, y - 0.17), (right - min(draw_width * 0.22, 250), y - 0.17),
                              (right - min(draw_width * 0.22, 250), y - 0.28), (right, y),
                              (right - min(draw_width * 0.22, 250), y + 0.28),
                              (right - min(draw_width * 0.22, 250), y + 0.17), (left, y + 0.17)]
                elif plot_strand == "-":
                    left, right = x1 - draw_width, x1
                    color_target = polys_minus
                    points = [(right, y - 0.17), (left + min(draw_width * 0.22, 250), y - 0.17),
                              (left + min(draw_width * 0.22, 250), y - 0.28), (left, y),
                              (left + min(draw_width * 0.22, 250), y + 0.28),
                              (left + min(draw_width * 0.22, 250), y + 0.17), (right, y + 0.17)]
                else:
                    color_target = polys_unknown
                    points = [(x0, y - 0.12), (x0 + draw_width, y - 0.12),
                              (x0 + draw_width, y + 0.12), (x0, y + 0.12)]
                color_target.append(points)
                # Text labels are deliberately limited to keep whole-genome views readable.
                if len(gene.label) <= 16:
                    gene_label_candidates.append((x0 + real_width / 2, y + 0.28, gene.label))
            # Dotted separators between contigs are restricted to this genome row.
            if c_rank < len(contigs):
                ax.plot([offset + contig.length + gap / 2] * 2, [y - 0.38, y + 0.38], color="#777777", linestyle=":", linewidth=0.55)
            offset += contig.length
            if c_rank < len(contigs):
                offset += gap
        all_x_max = max(all_x_max, offset)
        if polys_plus:
            ax.add_collection(PolyCollection(polys_plus, facecolors=plus_color, edgecolors="none", alpha=0.88, rasterized=True))
        if polys_minus:
            ax.add_collection(PolyCollection(polys_minus, facecolors=minus_color, edgecolors="none", alpha=0.88, rasterized=True))
        if polys_unknown:
            ax.add_collection(PolyCollection(polys_unknown, facecolors=unknown_color, edgecolors="none", alpha=0.88, rasterized=True))
        if len(gene_label_candidates) > max_labels_per_genome:
            notes.append(
                f"Gene track labels for {genome.name} were limited to {max_labels_per_genome} "
                f"of {len(gene_label_candidates)} eligible labels."
            )
        if max_labels_per_genome > 0 and gene_label_candidates:
            if len(gene_label_candidates) > max_labels_per_genome:
                indices = [round(i * (len(gene_label_candidates) - 1) / (max_labels_per_genome - 1))
                           for i in range(max_labels_per_genome)] if max_labels_per_genome > 1 else [0]
                gene_label_candidates = [gene_label_candidates[i] for i in indices]
            for x, yy, label in gene_label_candidates:
                ax.text(x, yy, label, fontsize=4.5, rotation=55, ha="left", va="bottom", clip_on=True)
    ax.set_xlim(0, max(all_x_max, 1))
    ax.set_ylim(-0.7, len(displayed) - 0.3)
    ax.set_yticks(range(len(displayed)))
    ax.set_yticklabels([g.name for g in displayed][::-1], fontsize=7)
    ax.set_xlabel("Concatenated contig coordinates (contigs sorted longest-first; dashed marks separate contigs)")
    ax.set_ylabel("Genome")
    ax.grid(True, axis="x", linestyle=":", linewidth=0.5, alpha=0.35)
    ax.legend(handles=[Line2D([0], [0], color=plus_color, lw=5, label="+ strand"),
                       Line2D([0], [0], color=minus_color, lw=5, label="− strand"),
                       Line2D([0], [0], color=unknown_color, lw=5, label="Unknown strand")],
              loc="upper right", frameon=False, fontsize=8, ncol=3)
    subtitle = "reference-normalised direction" if reference_genome else "raw GFF direction"
    ax.set_title(f"Genome gene-order tracks ({subtitle})")
    fig.tight_layout()
    fig.savefig(plot_dir / "gene_order_tracks.png", dpi=220, bbox_inches="tight")
    fig.savefig(plot_dir / "gene_order_tracks.svg", bbox_inches="tight")
    plt.close(fig)

    # Plot 3: reference-relative dotplots from unique family anchors, if available.
    reference = next((g for g in genomes if reference_genome and g.name.casefold() == reference_genome.casefold()), None)
    if reference:
        ref_counts = Counter(g.family for g in reference.genes if g.family)
        plot_genomes = [g for g in displayed if g.name != reference.name]
        plot_genomes = plot_genomes[:24]
        if plot_genomes:
            ncols = 4
            nrows = math.ceil(len(plot_genomes) / ncols)
            fig, axes = plt.subplots(nrows, ncols, figsize=(4.0 * ncols, 3.1 * nrows), squeeze=False)
            for ax, genome in zip(axes.flat, plot_genomes):
                xs, ys, colors = [], [], []
                query_family_counts = Counter(g.family for g in genome.genes if g.family)
                cmap = plt.get_cmap("tab20")
                contigs = sorted(genome.contigs.values(), key=lambda c: (-c.length, c.name.casefold(), c.input_order))
                for contig_rank, contig in enumerate(contigs):
                    if not contig.reference_contig or contig.orientation == "undetermined":
                        continue
                    ref_contig = reference.contigs.get(contig.reference_contig)
                    if not ref_contig:
                        continue
                    ref_genes = sorted(ref_contig.genes, key=lambda g: (g.start, g.end, g.feature_id.casefold()))
                    ref_family_rank = {g.family: i for i, g in enumerate(ref_genes, 1) if g.family and ref_counts[g.family] == 1}
                    q_genes = sorted(contig.genes, key=lambda g: (g.start, g.end, g.feature_id.casefold()))
                    if contig.orientation == "reverse":
                        q_genes = q_genes[::-1]
                    q_count = max(len(q_genes) - 1, 1)
                    r_count = max(len(ref_genes) - 1, 1)
                    for qi, gene in enumerate(q_genes):
                        if not gene.family or gene.family not in ref_family_rank or query_family_counts[gene.family] != 1:
                            continue
                        xs.append((ref_family_rank[gene.family] - 1) / r_count)
                        ys.append(qi / q_count)
                        colors.append(cmap(contig_rank % 20))
                if xs:
                    ax.scatter(xs, ys, s=7, alpha=0.7, c=colors, linewidths=0)
                    ax.plot([0, 1], [0, 1], linestyle="--", linewidth=0.7, color="#555555", alpha=0.7)
                ax.set_title(genome.name, fontsize=8)
                ax.set_xlabel("Reference contig-relative gene position", fontsize=7)
                ax.set_ylabel("Query contig-relative position (normalised orientation)", fontsize=7)
                ax.tick_params(labelsize=6)
                ax.set_xlim(-0.02, 1.02)
                ax.set_ylim(-0.02, 1.02)
                ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.3)
            for ax in axes.flat[len(plot_genomes):]:
                ax.axis("off")
            fig.suptitle(f"Family-anchor gene-order comparison against reference: {reference.name}", y=1.01)
            fig.tight_layout()
            fig.savefig(plot_dir / "reference_gene_order_dotplots.png", dpi=220, bbox_inches="tight")
            fig.savefig(plot_dir / "reference_gene_order_dotplots.svg", bbox_inches="tight")
            plt.close(fig)
    return notes


def make_breakpoint_support_plot(outdir: Path, bin_percent: int = 2) -> list[str]:
    """Show how often query assemblies have contig ends along the main reference replicon."""
    support_path = outdir / "assembly_breakpoint_support.tsv"
    if not support_path.exists():
        return []
    with support_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        return []
    _matplotlib, plt, _PolyCollection, _Line2D = _require_matplotlib()
    if bin_percent < 1 or bin_percent > 10:
        raise ValueError("bin_percent must be between 1 and 10")
    # The largest reference contig is generally the chromosome. Plotting only it
    # answers the cohort-level question without giving many small replicons equal visual weight.
    contig_sizes: dict[str, int] = defaultdict(int)
    for row in rows:
        contig_sizes[row["reference_contig"]] = max(
            contig_sizes[row["reference_contig"]], int(row["reference_gene_order_max"])
        )
    ref_contig = max(contig_sizes, key=lambda name: (contig_sizes[name], name.casefold()))
    selected = [row for row in rows if row["reference_contig"] == ref_contig]
    gene_order_path = outdir / "gene_order.tsv"
    if not gene_order_path.exists():
        return ["Assembly-end overview skipped because gene_order.tsv is unavailable for physical reference coordinates."]
    ref_genome = selected[0]["reference_genome"]
    reference_genes: list[tuple[int, float, int]] = []
    with gene_order_path.open("r", encoding="utf-8", newline="") as handle:
        for gene in csv.DictReader(handle, delimiter="\t"):
            if gene["genome"] == ref_genome and gene["contig"] == ref_contig:
                reference_genes.append((int(gene["gene_order_on_contig"]), float(gene["midpoint"]),
                                        int(gene["contig_length_bp"])))
    reference_genes.sort()
    if not reference_genes:
        return ["Assembly-end overview skipped because reference gene coordinates were not found."]
    ref_length = reference_genes[0][2]
    gene_count = len(reference_genes)

    def rank_to_percent(rank: float) -> float:
        # Interpolate median anchor ranks to physical bp positions on the reference.
        lower = max(1, min(gene_count, math.floor(rank)))
        upper = max(1, min(gene_count, math.ceil(rank)))
        low_bp = reference_genes[lower - 1][1]
        high_bp = reference_genes[upper - 1][1]
        bp = low_bp + (rank - lower) * (high_bp - low_bp)
        return bp / ref_length * 100

    n_bins = math.ceil(100 / bin_percent)
    genomes_by_bin: list[set[str]] = [set() for _ in range(n_bins)]
    for row in selected:
        position_percent = rank_to_percent(float(row["reference_gene_order_median"]))
        bin_index = min(n_bins - 1, int(position_percent // bin_percent))
        genomes_by_bin[bin_index].update(filter(None, row["genomes"].split(",")))
    counts = [len(genomes) for genomes in genomes_by_bin]
    centers = [(i + 0.5) * bin_percent for i in range(n_bins)]
    fig, ax = plt.subplots(figsize=(11, 3.8))
    ax.bar(centers, counts, width=bin_percent * 0.88, color="#2878A5", edgecolor="white", linewidth=0.35)
    peak = max(counts, default=0)
    if peak:
        peak_bin = counts.index(peak)
        ax.annotate(f"Peak: {peak}/{len({g for row in selected for g in row['genomes'].split(',') if g})} genomes",
                    (centers[peak_bin], peak), xytext=(0, 8), textcoords="offset points",
                    ha="center", fontsize=9, color="#7A2828")
    ax.set(xlim=(0, 100), ylim=(0, max(1, peak * 1.18)),
           xlabel="Reference chromosome position (%)",
           ylabel="Genomes with a contig end in bin",
           title=f"Reference-projected contig-end locations ({ref_contig})")
    ax.set_xticks(range(0, 101, 10))
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.55)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(0.01, 0.01, f"{len({g for row in selected for g in row['genomes'].split(',') if g})} query genomes; {bin_percent}% bins by physical bp; contig ends must pass configured edge limits.",
             fontsize=8, color="#555555")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(outdir / "figures" / "assembly_breakpoint_support.png", dpi=220, bbox_inches="tight")
    fig.savefig(outdir / "figures" / "assembly_breakpoint_support.svg", bbox_inches="tight")
    plt.close(fig)
    return [f"Assembly-end overview plots the largest reference replicon ({ref_contig}); full per-contig results remain in TSV reports."]


def write_run_summary(outdir: Path, genomes: list[Genome], counts: dict, warnings: list[str],
                      family_map: FamilyMap | None, reference_genome: str | None, args: argparse.Namespace,
                      plot_notes: list[str]) -> None:
    lines = [
        "PyamilySeq Genome Order Analysis",
        "================================",
        f"Version: {VERSION}",
        f"Run UTC: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        "",
        "Input and analysis summary",
        f"Genomes parsed: {counts['n_genomes']}",
        f"Contigs parsed: {counts['n_contigs']}",
        f"Genes/loci parsed: {counts['n_genes']}",
        f"Genes mapped to a family: {counts['n_family_mapped_genes']}",
        f"Contigs with reference-projected breakpoints: {counts.get('n_contig_breakpoints_mapped', 0)}",
        f"Candidate recurrent reference-projected assembly-end regions (descriptive, not significance-tested): {counts.get('n_recurrent_breakpoint_regions', 0)}",
        f"Genes near contig ends or with partial/disruption annotation flags: {counts.get('n_contig_end_gene_candidates', 0)}",
        f"Reference genome: {reference_genome or 'not set'}",
        f"Family map: {family_map.source if family_map else 'not supplied'}",
        f"Genomes with parser warnings: {counts['n_genomes_with_warnings']}",
        "",
        "Interpretation notes",
        "- Contig ranks are assigned within each genome by descending contig length, with contig ID as a deterministic tie-breaker.",
        "- Endpoint support requires FASTA or explicit GFF sequence length. Contigs sized only by maximum observed feature end are excluded from endpoint-support counts.",
        "- A GFF-only run describes gene coordinates and strands, but cannot establish homologous contig orientation between genomes.",
        "- Reference-relative orientation uses families present once per genome and the reference contig with the most shared anchors; the configured minimum anchor count and absolute Kendall tau are required, and tied reference-contig matches remain undetermined.",
        "- Contig ranks are assembly descriptions, not evidence that contigs occupy the same biological position across genomes.",
        "- Order metrics compare each query contig with one reference contig. They do not compare arbitrary concatenated contig orders as if the assembly were contiguous.",
        "- Gene-span fraction sums annotated feature lengths; it is not the same as coding density when annotations overlap or selected features are not CDS.",
        "",
        "Files written",
        "- genome_summary.tsv: per-genome assembly and annotation summary.",
        "- contig_summary.tsv: per-contig lengths, genes, strands, and inferred orientation.",
        "- contig_order.tsv: contigs ranked by length within each genome.",
        "- gene_order.tsv: per-locus coordinates, strand, order, labels, and family mapping.",
        "- contig_orientation.tsv: evidence and decision for reference-relative orientation.",
        "- gene_order_metrics_vs_reference.tsv: strand-aware order-concordance metrics when a family map and reference are supplied.",
        "- family_order_matches.tsv: matched single-copy family anchors and their order in query/reference contigs.",
        "- contig_breakpoints_vs_reference.tsv: contig ends projected onto reference gene order; eligible endpoints require known contig lengths.",
        "- assembly_breakpoint_support.tsv: descriptive clusters of nearby reference-projected assembly ends; not a significance test.",
        "- contig_end_gene_candidates.tsv: genes near physical contig ends or carrying annotation flags; candidates require sequence-level validation.",
        "- figures/: contig size profiles, directional genome tracks, reference dotplots, and breakpoint support plots when applicable.",
    ]
    if not family_map or counts.get("n_family_mapped_genes", 0) == 0:
        lines.extend([
            "",
            "Family-based reports unavailable",
            "- gene_order_metrics_vs_reference.tsv and family_order_matches.tsv contain headers only because no gene-family map was supplied or mapped.",
            "- contig_breakpoints_vs_reference.tsv and assembly_breakpoint_support.tsv contain headers only because reference-projected breakpoints require a family map and reference genome.",
            "- Supply a PyamilySeq gene-family map and --reference-genome to populate these reports and create reference comparison plots.",
        ])
    if warnings:
        lines.extend(["", "Run-level warnings", *[f"- {w}" for w in warnings]])
    genome_warnings = [f"{g.name}: {warning}" for g in genomes for warning in g.warnings]
    if genome_warnings:
        lines.extend(["", "GFF parsing warnings", *[f"- {w}" for w in genome_warnings[:500]]])
        if len(genome_warnings) > 500:
            lines.append(f"- ... {len(genome_warnings) - 500} additional GFF warnings omitted from this text summary; warnings are reflected in genome_summary.tsv.")
    if plot_notes:
        lines.extend(["", "Figure notes", *[f"- {w}" for w in plot_notes]])
    (outdir / "run_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest = {
        "tool": "PyamilySeq Genome Order Analysis", "version": VERSION,
        "run_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "arguments": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "counts": counts,
        "genomes": [{"name": g.name, "gff": str(g.path), "feature_type": g.selected_feature_type,
                      "n_contigs": len(g.contigs), "n_genes": len(g.genes), "warnings": len(g.warnings)} for g in genomes],
        "family_map": asdict(family_map) if family_map else None,
        "run_warnings": warnings,
    }
    # Keep dictionary lookup contents out of the manifest; only record mapping provenance/settings.
    if manifest["family_map"]:
        manifest["family_map"] = {k: v for k, v in manifest["family_map"].items()
                                   if k not in {"by_gene", "by_genome_gene", "ambiguous_gene_keys"}}
    (outdir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="Track-Genome-Order",
        description="Summarise contig architecture and plot strand-aware gene order from a collection of GFF/GFF3 files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--gff-dir", required=True, type=Path, help="Directory containing one GFF/GFF3 file per genome (recursive).")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory (created if necessary).")
    parser.add_argument("--pattern", default="*", help="Recursive file-name glob applied under --gff-dir; suffix must still be a recognised GFF/GFF3 suffix.")
    parser.add_argument("--feature-type", default="auto", help="Feature type in GFF column 3. 'auto' prefers CDS, then gene, then RNA-like features.")
    parser.add_argument("--fasta-dir", type=Path, help="Optional FASTA directory for exact contig lengths; file stems should match GFF stems.")
    parser.add_argument("--assembly-report-map", type=Path,
                        help="Optional TSV with genome and assembly_report columns; relative report paths resolve from the TSV directory.")
    parser.add_argument("--subcohort-map", type=Path, help="Optional TSV with genome and subcohort columns; every input genome must be assigned exactly once.")
    parser.add_argument("--reference-map", type=Path, help="Optional TSV with subcohort and reference_genome columns; exactly one explicit reference is required per subcohort.")
    parser.add_argument("--minimap2-preset", choices=("asm5", "asm10", "asm20"), default="asm10",
                        help="Explicit minimap2 assembly preset used for nucleotide contig placement; a per-subcohort value in --reference-map overrides this.")
    parser.add_argument("--minimap2", default="minimap2", help="Path or executable name for minimap2.")
    parser.add_argument("--threads", type=int, default=1, help="Threads used by minimap2.")
    parser.add_argument("--placement-min-aligned-fraction", type=float, default=0.8,
                        help="Minimum query fraction covered for a confident placement (0 to 1].")
    parser.add_argument("--placement-min-mapq", type=int, default=20,
                        help="Minimum minimap2 MAPQ for a confident placement; MAPQ is not a calibrated probability.")
    parser.add_argument("--placement-min-block-length", type=int, default=1000,
                        help="Minimum alignment block length considered for confident placement.")
    parser.add_argument("--assembly-metadata", type=Path, help="Optional TSV with genome and optional assembly_level/completeness columns for candidate review.")
    parser.add_argument("--similarity-matrix", type=Path, help="Optional labeled square Sourmash similarity CSV used only for within-subcohort diagnostics.")
    parser.add_argument("--family-map", type=Path, help="Optional PyamilySeq gene-family assignment table or gene_presence_absence.csv.")
    parser.add_argument("--map-gene-column", help="Gene/member column in a long-format family map; inferred by default.")
    parser.add_argument("--map-family-column", help="Family/cluster column in a long-format family map; inferred by default.")
    parser.add_argument("--map-genome-column", help="Optional genome column in a long-format family map; inferred by default.")
    parser.add_argument("--sequence-delimiter", default="|", help="Delimiter used by combined gene IDs such as genome|gene.")
    parser.add_argument("--reference-genome", help="Explicit legacy single-cohort orientation/reference genome. No reference is selected automatically.")
    parser.add_argument("--min-orientation-anchors", type=int, default=5, help="Minimum unique shared families needed to infer contig orientation.")
    parser.add_argument("--min-orientation-tau", type=float, default=0.5,
                        help="Minimum absolute Kendall tau required for orientation inference (range 0 to 1); higher values are more conservative.")
    parser.add_argument("--breakpoint-window-genes", type=int, default=10, help="Maximum reference-gene-rank span for clustering nearby contig-end observations.")
    parser.add_argument("--max-end-anchor-genes", type=int, default=5, help="Maximum anchor position from a physical contig end, counting the anchor gene itself (positions 1 through N are eligible); applies to reference-anchor eligibility only, not physical endpoint inventory.")
    parser.add_argument("--max-end-anchor-bp", type=int, default=10000, help="Only count a contig endpoint for recurrence if its nearest unique-family anchor is within this many bp of the physical end.")
    parser.add_argument("--label-by", choices=("auto", "family", "locus_tag", "Name", "gene", "product", "ID"), default="auto", help="Gene label used in tabular output and figures.")
    parser.add_argument("--max-genomes-plotted", type=int, default=50, help="Maximum number of genomes drawn in overview figures; all genomes remain in reports.")
    parser.add_argument("--max-labels-per-genome", type=int, default=20, help="Label genes on a track only when that genome has no more than this many genes.")
    parser.add_argument("--no-plots", action="store_true", help="Write tabular reports only; does not require matplotlib.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return parser


def clear_stale_outputs(outdir: Path) -> None:
    """Remove only report/figure filenames owned by this tool to prevent stale plots."""
    reports = {
        "genome_summary.tsv", "contig_summary.tsv", "contig_order.tsv", "gene_order.tsv",
        "contig_orientation.tsv", "gene_order_metrics_vs_reference.tsv", "family_order_matches.tsv",
        "contig_breakpoints_vs_reference.tsv", "assembly_breakpoint_support.tsv", "run_summary.txt", "run_manifest.json",
        "contig_end_gene_candidates.tsv", "physical_endpoints.tsv", "replicon_manifest.tsv",
        "reference_candidates.tsv", "reference_manifest.tsv",
        "contig_placement_blocks.tsv", "contig_placement_summary.tsv",
        "placement_reference_sequence_manifest.tsv", "placement_query_sequence_manifest.tsv",
    }
    figures = {
        "contig_size_profiles.png", "contig_size_profiles.svg", "gene_order_tracks.png", "gene_order_tracks.svg",
        "reference_gene_order_dotplots.png", "reference_gene_order_dotplots.svg",
        "assembly_breakpoint_support.png", "assembly_breakpoint_support.svg",
    }
    outdir.mkdir(parents=True, exist_ok=True)
    for name in reports:
        path = outdir / name
        if path.is_file():
            path.unlink()
    figure_dir = outdir / "figures"
    for name in figures:
        path = figure_dir / name
        if path.is_file():
            path.unlink()
    placement_dir = outdir / "placement"
    if placement_dir.is_dir():
        shutil.rmtree(placement_dir)


def run(args: argparse.Namespace) -> int:
    gff_paths = discover_gffs(args.gff_dir, args.pattern)
    if not gff_paths:
        raise FileNotFoundError(f"No recognised GFF/GFF3 files found under {args.gff_dir} using pattern {args.pattern!r}")
    genomes: list[Genome] = []
    seen_names: set[str] = set()
    for path in gff_paths:
        genome = parse_gff(path, args.feature_type)
        if genome.name.casefold() in seen_names:
            raise ValueError(f"Genome filename stems are not unique: {genome.name!r}. Rename inputs or use one file per genome with unique stems.")
        seen_names.add(genome.name.casefold())
        genomes.append(genome)
    genomes.sort(key=lambda g: g.name.casefold())
    run_warnings: list[str] = []
    run_warnings.extend(apply_fasta_lengths(genomes, args.fasta_dir))
    if args.assembly_report_map and not args.fasta_dir:
        raise ValueError("--assembly-report-map requires --fasta-dir")
    replicon_metadata = collect_contig_replicon_metadata(
        genomes, args.fasta_dir, args.assembly_report_map
    ) if args.fasta_dir else []
    family_map = load_family_map(args.family_map, args.map_gene_column, args.map_family_column,
                                 args.map_genome_column, args.sequence_delimiter) if args.family_map else None
    mapped_count = assign_families(genomes, family_map, args.sequence_delimiter)
    if family_map and mapped_count == 0:
        run_warnings.append("A family map was supplied, but no GFF loci matched it. Check genome stems, gene identifiers, and --sequence-delimiter.")
    elif family_map and mapped_count < sum(len(g.genes) for g in genomes):
        total_loci = sum(len(g.genes) for g in genomes)
        run_warnings.append(
            f"Family-map coverage is partial: {mapped_count}/{total_loci} loci mapped "
            f"({mapped_count / total_loci:.1%}); inspect gene_order.tsv for unmatched loci."
        )
    reference_genome = args.reference_genome
    orientation_rows: list[dict] = []
    if family_map and mapped_count and args.reference_genome:
        reference_genome, orientation_rows = infer_contig_orientations(
            genomes, reference_genome, args.min_orientation_anchors, args.min_orientation_tau
        )
        contig_ranks = {
            genome.name: {contig.name: rank for rank, contig in enumerate(
                sorted(genome.contigs.values(), key=lambda x: (-x.length, x.name.casefold(), x.input_order)), 1)}
            for genome in genomes
        }
        for row in orientation_rows:
            row["contig_rank_by_length"] = contig_ranks.get(row["genome"], {}).get(row["contig"], "")
        n_undetermined = sum(1 for g in genomes if g.name != reference_genome for c in g.contigs.values() if c.orientation == "undetermined")
        if n_undetermined:
            run_warnings.append(f"{n_undetermined} contigs could not be confidently oriented against reference {reference_genome!r}; their original coordinate direction is retained.")
    else:
        reference_genome = None
        for genome in genomes:
            for contig in genome.contigs.values():
                contig.orientation = "not_assessed"
                contig.orientation_status = "no_family_map" if not family_map else "no_families_mapped"
                contig.reference_contig = ""
                for gene in contig.genes:
                    gene.orientation = contig.orientation
                    gene.orientation_reference_contig = ""
        orientation_rows = []
        for genome in genomes:
            ranked_contigs = sorted(genome.contigs.values(), key=lambda x: (-x.length, x.name.casefold(), x.input_order))
            for rank, contig in enumerate(ranked_contigs, 1):
                orientation_rows.append({
                    "genome": genome.name, "contig": contig.name, "contig_length": contig.length,
                    "contig_rank_by_length": rank, "reference_genome": "", "reference_contig": "",
                    "orientation": contig.orientation, "orientation_status": contig.orientation_status,
                    "unique_shared_families": 0, "kendall_tau_raw": "",
                })
    prepare_ordered_genomes(genomes, args.label_by)
    # Fail before clearing prior outputs when the optional plotting dependency is absent.
    if not args.no_plots:
        _require_matplotlib()
    report_args = (genomes, replicon_metadata, args.outdir, args.fasta_dir,
                   args.subcohort_map, args.reference_map, args.assembly_metadata,
                   args.similarity_matrix, args.reference_genome)
    build_reference_reports(*report_args, write_outputs=False)
    clear_stale_outputs(args.outdir)
    reference_config = build_reference_reports(*report_args)
    from .nucleotide_placement import run_nucleotide_placement
    placement_result = run_nucleotide_placement(
        genomes, args.fasta_dir, args.outdir, reference_config["assignments"],
        reference_config["selected_references"],
        {cohort: reference_config["presets"].get(cohort, args.minimap2_preset)
         for cohort in reference_config["selected_references"]},
        replicon_metadata, args.minimap2, args.threads,
        args.placement_min_aligned_fraction, args.placement_min_mapq,
        args.placement_min_block_length,
    )
    run_warnings.extend(placement_result.get("warnings", []))
    counts = build_reports(genomes, args.outdir, family_map, reference_genome, orientation_rows,
                           args.breakpoint_window_genes, args.max_end_anchor_genes, args.max_end_anchor_bp,
                           replicon_metadata)
    plot_notes: list[str] = []
    plot_failure: RuntimeError | None = None
    if not args.no_plots:
        try:
            plot_notes = make_plots(genomes, args.outdir, reference_genome, args.max_genomes_plotted, args.max_labels_per_genome)
            plot_notes.extend(make_breakpoint_support_plot(args.outdir))
        except RuntimeError as exc:
            plot_failure = exc
            plot_notes.append(f"Plot generation failed: {exc}")
    write_run_summary(args.outdir, genomes, counts, run_warnings, family_map, reference_genome, args, plot_notes)
    print(f"Genome Order Analysis complete: {counts['n_genomes']} genomes, {counts['n_contigs']} contigs, {counts['n_genes']} loci.")
    print(f"Reports: {args.outdir}")
    if reference_genome:
        print(f"Reference genome: {reference_genome}")
    if run_warnings:
        print(f"Warnings: {len(run_warnings)} (see run_summary.txt)")
    if plot_failure:
        print(f"WARNING: reports were written, but plot generation failed: {plot_failure}", file=sys.stderr)
        return 2
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.min_orientation_anchors < 2:
        parser.error("--min-orientation-anchors must be at least 2")
    if not 0 <= args.min_orientation_tau <= 1:
        parser.error("--min-orientation-tau must be between 0 and 1")
    if args.breakpoint_window_genes < 1:
        parser.error("--breakpoint-window-genes must be at least 1")
    if args.max_end_anchor_genes < 0:
        parser.error("--max-end-anchor-genes cannot be negative")
    if args.max_end_anchor_bp < 0:
        parser.error("--max-end-anchor-bp cannot be negative")
    if args.max_genomes_plotted < 1:
        parser.error("--max-genomes-plotted must be at least 1")
    if args.max_labels_per_genome < 0:
        parser.error("--max-labels-per-genome cannot be negative")
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    if not 0 < args.placement_min_aligned_fraction <= 1:
        parser.error("--placement-min-aligned-fraction must be in (0, 1]")
    if args.placement_min_mapq < 0:
        parser.error("--placement-min-mapq cannot be negative")
    if args.placement_min_block_length < 1:
        parser.error("--placement-min-block-length must be at least 1")
    try:
        return run(args)
    except (FileNotFoundError, NotADirectoryError, ValueError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
