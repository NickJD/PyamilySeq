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
        genome.contigs[contig_name] = Contig(
            name=contig_name, length=length, length_source=source, input_order=first_seen
        )

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
                active = line[1:].strip().split()[0] if line[1:].strip() else ""
                if not active:
                    raise ValueError(f"Empty FASTA header in {path}")
                length = 0
            else:
                length += len(line.strip())
        if active is not None:
            if active in lengths:
                raise ValueError(f"Duplicate FASTA contig identifier {active!r} in {path}")
            lengths[active] = length
    return lengths, None


def apply_fasta_lengths(genomes: list[Genome], fasta_dir: Path | None) -> list[str]:
    warnings: list[str] = []
    if fasta_dir is None:
        return warnings
    if not fasta_dir.exists() or not fasta_dir.is_dir():
        raise NotADirectoryError(f"--fasta-dir must be an existing directory: {fasta_dir}")
    for genome in genomes:
        lengths, warning = load_fasta_lengths(genome, fasta_dir)
        if warning:
            warnings.append(warning)
            continue
        for contig_name, contig in genome.contigs.items():
            if contig_name in lengths:
                contig.length = lengths[contig_name]
                contig.length_source = "fasta"
            else:
                warnings.append(f"{genome.name}: contig {contig_name!r} absent from matching FASTA; retained GFF-derived length")
        gff_contigs = set(genome.contigs)
        extra_fasta_contigs = sorted(set(lengths) - gff_contigs, key=str.casefold)
        if extra_fasta_contigs:
            for contig_name in extra_fasta_contigs:
                genome.contigs[contig_name] = Contig(
                    name=contig_name, length=lengths[contig_name],
                    length_source="fasta_no_gff_features",
                )
            warnings.append(
                f"{genome.name}: FASTA contains {len(extra_fasta_contigs)} contig(s) without GFF features; "
                "included in assembly contig counts and length/N50 summaries, but not family-anchor comparisons "
                f"(examples: {', '.join(extra_fasta_contigs[:10])})"
            )
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
    lookup = {g.name.casefold(): g for g in genomes}
    reference = lookup.get((reference_name or "").casefold()) if reference_name else genomes[0]
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
                "left_boundary_support_eligible": (left_idx <= max_end_anchor_genes and
                    (left_gene.start - 1 if contig.orientation == "forward" else max(0, contig.length - left_gene.end)) <= max_end_anchor_bp),
                "right_boundary_anchor_family": right_gene.family, "right_boundary_anchor_gene_id": right_gene.feature_id,
                "right_boundary_reference_gene_order": right_rank,
                "genes_from_right_contig_end_to_anchor": len(ordered_genes) - 1 - right_idx,
                "bp_from_right_contig_end_to_anchor": (max(0, contig.length - right_gene.end) if contig.orientation == "forward"
                                                        else right_gene.start - 1),
                "right_boundary_support_eligible": ((len(ordered_genes) - 1 - right_idx) <= max_end_anchor_genes and
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
                                  "eligible": (left_idx <= max_end_anchor_genes and
                                      (left_gene.start - 1 if contig.orientation == "forward" else max(0, contig.length - left_gene.end)) <= max_end_anchor_bp),
                                  "near_reference_contig_end": left_rank <= window_genes})
            endpoint_rows.append({"genome": genome.name, "contig": contig.name, "reference_contig": ref_contig.name,
                                  "boundary": "right", "reference_gene_order": right_rank,
                                  "family": right_gene.family, "anchor_gene_id": right_gene.feature_id,
                                  "genes_from_contig_end_to_anchor": len(ordered_genes) - 1 - right_idx,
                                  "bp_from_contig_end_to_anchor": (max(0, contig.length - right_gene.end) if contig.orientation == "forward"
                                                                    else right_gene.start - 1),
                                  "eligible": ((len(ordered_genes) - 1 - right_idx) <= max_end_anchor_genes and
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
                "interpretation": "Candidate recurrent contig-end region; investigate assembly and biological alternatives.",
            })
    return breakpoint_rows, support_rows


def _write_tsv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        if fieldnames:
            writer.writeheader()
            writer.writerows(rows)


def build_reports(genomes: list[Genome], outdir: Path, family_map: FamilyMap | None,
                  reference_genome: str | None, orientation_rows: list[dict], breakpoint_window_genes: int = 10,
                  max_end_anchor_genes: int = 5, max_end_anchor_bp: int = 10000) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    genome_rows: list[dict] = []
    contig_rows: list[dict] = []
    gene_rows: list[dict] = []
    contig_order_rows: list[dict] = []
    comparison_rows: list[dict] = []
    order_metric_rows: list[dict] = []
    reference = next((g for g in genomes if reference_genome and g.name.casefold() == reference_genome.casefold()), None)
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
            "n_contigs": len(contigs), "n_genes": len(genome.genes), "total_contig_length_bp": genome_length,
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

    return {
        "n_genomes": len(genomes), "n_contigs": sum(len(g.contigs) for g in genomes),
        "n_genes": sum(len(g.genes) for g in genomes), "n_family_mapped_genes": sum(1 for g in genomes for x in g.genes if x.family),
        "reference_genome": reference.name if reference else "", "n_contig_orientation_records": len(orientation_rows),
        "n_family_order_matches": len(comparison_rows), "n_contig_breakpoints_mapped": len(breakpoint_rows),
        "n_recurrent_breakpoint_regions": sum(1 for r in breakpoint_support_rows if r["candidate_recurrent_region"]),
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
        return ["Breakpoint overview skipped because gene_order.tsv is unavailable for physical reference coordinates."]
    ref_genome = selected[0]["reference_genome"]
    reference_genes: list[tuple[int, float, int]] = []
    with gene_order_path.open("r", encoding="utf-8", newline="") as handle:
        for gene in csv.DictReader(handle, delimiter="\t"):
            if gene["genome"] == ref_genome and gene["contig"] == ref_contig:
                reference_genes.append((int(gene["gene_order_on_contig"]), float(gene["midpoint"]),
                                        int(gene["contig_length_bp"])))
    reference_genes.sort()
    if not reference_genes:
        return ["Breakpoint overview skipped because reference gene coordinates were not found."]
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
           title=f"Where assemblies break along the reference chromosome ({ref_contig})")
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
    return [f"Breakpoint overview plots the largest reference replicon ({ref_contig}); full per-contig results remain in TSV reports."]


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
        f"Candidate recurrent breakpoint regions (seen in at least two genomes): {counts.get('n_recurrent_breakpoint_regions', 0)}",
        f"Reference genome: {reference_genome or 'not set'}",
        f"Family map: {family_map.source if family_map else 'not supplied'}",
        f"Genomes with parser warnings: {counts['n_genomes_with_warnings']}",
        "",
        "Interpretation notes",
        "- Contig ranks are assigned within each genome by descending contig length, with contig ID as a deterministic tie-breaker.",
        "- Length source is recorded per contig. Without a FASTA or explicit GFF sequence length, length is estimated from the maximum observed feature end and may be too short.",
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
        "- contig_breakpoints_vs_reference.tsv: inferred contig ends projected onto reference gene order.",
        "- assembly_breakpoint_support.tsv: clusters of nearby contig ends recurring across genomes.",
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
    parser.add_argument("--family-map", type=Path, help="Optional PyamilySeq gene-family assignment table or gene_presence_absence.csv.")
    parser.add_argument("--map-gene-column", help="Gene/member column in a long-format family map; inferred by default.")
    parser.add_argument("--map-family-column", help="Family/cluster column in a long-format family map; inferred by default.")
    parser.add_argument("--map-genome-column", help="Optional genome column in a long-format family map; inferred by default.")
    parser.add_argument("--sequence-delimiter", default="|", help="Delimiter used by combined gene IDs such as genome|gene.")
    parser.add_argument("--reference-genome", help="Genome filename stem to use as orientation/order reference. Defaults to the first genome alphabetically when a family map is supplied.")
    parser.add_argument("--min-orientation-anchors", type=int, default=5, help="Minimum unique shared families needed to infer contig orientation.")
    parser.add_argument("--min-orientation-tau", type=float, default=0.5,
                        help="Minimum absolute Kendall tau required for orientation inference (range 0 to 1); higher values are more conservative.")
    parser.add_argument("--breakpoint-window-genes", type=int, default=10, help="Maximum reference-gene-rank span for clustering nearby contig-end observations.")
    parser.add_argument("--max-end-anchor-genes", type=int, default=5, help="Only count a contig endpoint for recurrence if its nearest unique-family anchor is within this many genes of the physical end.")
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
    if family_map and mapped_count:
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
    clear_stale_outputs(args.outdir)
    counts = build_reports(genomes, args.outdir, family_map, reference_genome, orientation_rows,
                           args.breakpoint_window_genes, args.max_end_anchor_genes, args.max_end_anchor_bp)
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
    try:
        return run(args)
    except (FileNotFoundError, NotADirectoryError, ValueError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
