"""Canonical gene-family parsing, reconciliation, and assignment models.

The public PyamilySeq outputs are deliberately generated from :class:`FamilyLedger`
rather than rebuilding family membership independently in each writer.  This keeps
the central invariant explicit: every unique input FASTA identifier is assigned to
exactly one final family.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
import os
import re
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Mapping


class FamilyInputError(ValueError):
    """Raised when clustering inputs cannot produce an unambiguous family ledger."""


def open_text(path, mode="rt"):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    return opener(path, mode, encoding="utf-8")


def fasta_identifier(defline: str) -> str:
    """Return the canonical FASTA identifier (the first whitespace token)."""

    identifier = defline.strip().split(maxsplit=1)[0] if defline.strip() else ""
    if not identifier:
        raise FamilyInputError("Encountered an empty FASTA header")
    return identifier


def read_fasta_records(path: str) -> OrderedDict[str, str]:
    """Read FASTA records and reject identifiers that cannot be represented safely."""

    records: OrderedDict[str, str] = OrderedDict()
    current = None
    chunks: list[str] = []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks)
                current = fasta_identifier(line[1:])
                if current in records:
                    raise FamilyInputError(
                        f"Duplicate FASTA identifier {current!r} at line {line_number} in {path}"
                    )
                chunks = []
            elif current is None:
                if line.strip():
                    raise FamilyInputError(
                        f"Sequence data before the first FASTA header at line {line_number} in {path}"
                    )
            else:
                chunks.append("".join(line.split()))
    if current is not None:
        if current in records:
            raise FamilyInputError(f"Duplicate FASTA identifier {current!r} in {path}")
        records[current] = "".join(chunks)
    if not records:
        raise FamilyInputError(f"No FASTA records found in {path}")
    empty = [identifier for identifier, sequence in records.items() if not sequence]
    if empty:
        sample = ", ".join(empty[:10])
        raise FamilyInputError(f"FASTA contains {len(empty)} empty sequence(s): {sample}")
    return records


def _stable_id(prefix: str, members: Iterable[str]) -> str:
    payload = "\0".join(sorted(set(members))).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(payload).hexdigest()[:12]}"


@dataclass(frozen=True)
class ParsedCluster:
    label: str
    members: tuple[str, ...]
    representative: str


@dataclass
class Family:
    family_id: str
    members: tuple[str, ...]
    representative: str
    taxa: tuple[str, ...]
    genomes: tuple[str, ...]
    stage: str
    source_families: tuple[str, ...] = ()
    source_clusters: tuple[str, ...] = ()
    category: str = ""

    @property
    def prevalence(self) -> int:
        return len(self.taxa)

    @property
    def copy_number(self) -> Counter:
        return Counter(self.genomes)


class TaxonResolver:
    """Resolve sequence identifiers to genomes and analysis taxa.

    Mapping files may contain either ``identifier<TAB>taxon`` or
    ``identifier<TAB>genome<TAB>taxon``.  Exact sequence mappings take
    precedence over genome-prefix mappings.
    """

    def __init__(self, mode="Species", mapping_path=None, sequence_delimiter="|"):
        self.mode = mode
        self.sequence_delimiter = sequence_delimiter
        self.mapping: dict[str, tuple[str | None, str]] = {}
        if mapping_path:
            self._load_mapping(mapping_path)

    def _load_mapping(self, path):
        with open_text(path) as handle:
            sample = handle.read(4096)
            handle.seek(0)
            delimiter = "\t" if "\t" in sample else ","
            reader = csv.reader(handle, delimiter=delimiter)
            for line_number, row in enumerate(reader, start=1):
                row = [cell.strip() for cell in row]
                if not row or not row[0] or row[0].startswith("#"):
                    continue
                if line_number == 1 and row[0].lower() in {
                    "identifier", "sequence", "sequence_id", "genome", "genome_id"
                }:
                    continue
                if len(row) == 2:
                    key, taxon = row
                    genome = None
                elif len(row) >= 3:
                    key, genome, taxon = row[:3]
                else:
                    raise FamilyInputError(
                        f"Taxon map line {line_number} in {path} has fewer than two columns"
                    )
                if not taxon:
                    raise FamilyInputError(f"Empty taxon on line {line_number} in {path}")
                value = (genome or None, taxon)
                if key in self.mapping and self.mapping[key] != value:
                    raise FamilyInputError(f"Conflicting taxon mappings for {key!r} in {path}")
                self.mapping[key] = value

    def _prefix(self, sequence_id):
        if self.sequence_delimiter and self.sequence_delimiter in sequence_id:
            return sequence_id.split(self.sequence_delimiter, 1)[0]
        return sequence_id

    def genome(self, sequence_id):
        prefix = self._prefix(sequence_id)
        for key in (sequence_id, prefix):
            mapped = self.mapping.get(key)
            if mapped and mapped[0]:
                return mapped[0]
        return prefix

    def taxon(self, sequence_id):
        prefix = self._prefix(sequence_id)
        for key in (sequence_id, prefix):
            mapped = self.mapping.get(key)
            if mapped:
                return mapped[1]
        genome = self.genome(sequence_id)
        if self.mode == "Genus":
            return genome.split("_", 1)[0]
        return genome


def _deduplicate_members(members: Iterable[str]) -> tuple[str, ...]:
    return tuple(OrderedDict.fromkeys(member for member in members if member))


def parse_cdhit(path: str) -> list[ParsedCluster]:
    clusters: list[ParsedCluster] = []
    label = None
    members: list[str] = []
    representative = None
    member_re = re.compile(r"^\s*\d+\s+[^,]+,\s+>(.*?)\.\.\.\s*(.*)$")

    def finish():
        if label is None:
            return
        unique = _deduplicate_members(members)
        if not unique:
            raise FamilyInputError(f"CD-HIT cluster {label!r} in {path} has no members")
        rep = representative or min(unique)
        if rep not in unique:
            raise FamilyInputError(f"Representative {rep!r} is not in CD-HIT cluster {label!r}")
        clusters.append(ParsedCluster(label, unique, rep))

    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                finish()
                match = re.match(r"^>Cluster\s+(.+?)\s*$", line, flags=re.IGNORECASE)
                if not match:
                    raise FamilyInputError(
                        f"Expected a CD-HIT '>Cluster N' header at line {line_number} in {path}; got {line!r}"
                    )
                label = match.group(1)
                members = []
                representative = None
                continue
            if label is None:
                raise FamilyInputError(f"CD-HIT member before a cluster header at line {line_number} in {path}")
            match = member_re.match(line)
            if not match:
                raise FamilyInputError(f"Malformed CD-HIT member line {line_number} in {path}: {line!r}")
            member = match.group(1)
            members.append(member)
            if "*" in match.group(2):
                if representative is not None and representative != member:
                    raise FamilyInputError(f"CD-HIT cluster {label!r} has multiple representatives")
                representative = member
    finish()
    if not clusters:
        raise FamilyInputError(f"No CD-HIT clusters found in {path}")
    labels = [cluster.label for cluster in clusters]
    if len(labels) != len(set(labels)):
        raise FamilyInputError(f"Duplicate CD-HIT cluster labels in {path}")
    return clusters


def _edge_rows(path: str, delimiter="\t"):
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            row = next(csv.reader([line], delimiter=delimiter))
            if len(row) < 2 or not row[0].strip() or not row[1].strip():
                raise FamilyInputError(
                    f"Expected at least two {delimiter!r}-separated columns at line {line_number} in {path}"
                )
            yield line_number, [cell.strip() for cell in row]


def parse_mmseqs(path: str, delimiter="\t") -> list[ParsedCluster]:
    grouped: OrderedDict[str, OrderedDict[str, None]] = OrderedDict()
    for _, row in _edge_rows(path, delimiter=delimiter):
        representative, member = row[:2]
        grouped.setdefault(representative, OrderedDict())
        grouped[representative][representative] = None
        grouped[representative][member] = None
    if not grouped:
        raise FamilyInputError(f"No MMseqs clusters found in {path}")
    return [
        ParsedCluster(rep, tuple(members), rep)
        for rep, members in grouped.items()
    ]


def _float_cell(row, index, name, line_number, path):
    try:
        return float(row[index])
    except (IndexError, ValueError):
        raise FamilyInputError(
            f"BLAST {name} is missing or invalid at line {line_number} in {path}"
        ) from None


def parse_blast(path: str, options=None, delimiter="\t") -> list[ParsedCluster]:
    """Parse BLAST edges with optional filtering and reciprocal-best-hit linkage.

    Two-column files are treated as an explicitly filtered undirected edge list.
    Scored BLAST files default to reciprocal-best-hit edges, which prevents one
    weak transitive hit from joining otherwise distinct families.
    """

    columns_text = getattr(options, "blast_columns", None) or (
        "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore"
    )
    column_names = [name.strip().lower() for name in columns_text.split(",")]
    column_index = {name: idx for idx, name in enumerate(column_names)}
    min_identity = float(getattr(options, "blast_min_identity", 0.0) or 0.0)
    max_evalue = float(getattr(options, "blast_max_evalue", math.inf))
    min_query_coverage = float(getattr(options, "blast_min_query_coverage", 0.0) or 0.0)
    min_subject_coverage = float(getattr(options, "blast_min_subject_coverage", 0.0) or 0.0)
    graph_mode = getattr(options, "blast_graph_mode", "reciprocal-best")

    nodes: set[str] = set()
    hits: list[tuple[str, str, float, int]] = []
    scored = False
    for line_number, row in _edge_rows(path, delimiter=delimiter):
        query, subject = row[:2]
        nodes.update((query, subject))
        if len(row) == 2:
            hits.append((query, subject, 1.0, line_number))
            continue
        scored = True
        if "pident" in column_index:
            identity = _float_cell(row, column_index["pident"], "pident", line_number, path)
            if identity < min_identity:
                continue
        if "evalue" in column_index:
            evalue = _float_cell(row, column_index["evalue"], "evalue", line_number, path)
            if evalue > max_evalue:
                continue
        alignment_length = None
        if "length" in column_index:
            alignment_length = _float_cell(row, column_index["length"], "length", line_number, path)
        if min_query_coverage or min_subject_coverage:
            if alignment_length is None or "qlen" not in column_index or "slen" not in column_index:
                raise FamilyInputError(
                    "BLAST coverage filtering requires length, qlen, and slen in --blast-columns"
                )
            query_length = _float_cell(row, column_index["qlen"], "qlen", line_number, path)
            subject_length = _float_cell(row, column_index["slen"], "slen", line_number, path)
            query_coverage = 100.0 * alignment_length / query_length if query_length else 0.0
            subject_coverage = 100.0 * alignment_length / subject_length if subject_length else 0.0
            if query_coverage < min_query_coverage or subject_coverage < min_subject_coverage:
                continue
        if "bitscore" in column_index:
            score = _float_cell(row, column_index["bitscore"], "bitscore", line_number, path)
        elif "evalue" in column_index:
            score = -math.log10(max(evalue, 1e-300))
        else:
            score = 1.0
        hits.append((query, subject, score, line_number))

    adjacency: dict[str, set[str]] = {node: set() for node in nodes}
    if scored and graph_mode == "reciprocal-best":
        best_score: dict[str, float] = {}
        best_subjects: dict[str, set[str]] = defaultdict(set)
        for query, subject, score, _ in hits:
            if query == subject:
                continue
            if query not in best_score or score > best_score[query]:
                best_score[query] = score
                best_subjects[query] = {subject}
            elif score == best_score[query]:
                best_subjects[query].add(subject)
        for query, subject, _, _ in hits:
            if query == subject:
                continue
            if subject in best_subjects.get(query, ()) and query in best_subjects.get(subject, ()):
                adjacency[query].add(subject)
                adjacency[subject].add(query)
    else:
        for query, subject, _, _ in hits:
            if query != subject:
                adjacency[query].add(subject)
                adjacency[subject].add(query)

    clusters: list[ParsedCluster] = []
    visited: set[str] = set()
    weighted_degree = Counter()
    for query, subject, score, _ in hits:
        if subject in adjacency.get(query, ()):
            weighted_degree[query] += score
            weighted_degree[subject] += score
    for seed in sorted(nodes):
        if seed in visited:
            continue
        stack = [seed]
        component = []
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            component.append(node)
            stack.extend(sorted(adjacency[node] - visited, reverse=True))
        component.sort()
        representative = min(component, key=lambda node: (-weighted_degree[node], node))
        clusters.append(ParsedCluster(str(len(clusters)), tuple(component), representative))
    if not clusters:
        raise FamilyInputError(f"No BLAST nodes found in {path}")
    return clusters


def parse_cluster_file(path: str, cluster_format: str, options=None) -> list[ParsedCluster]:
    normalised = cluster_format.strip().upper().replace("_", "-")
    if normalised in {"CD-HIT", "CDHIT", "CD-HIT-EST"}:
        return parse_cdhit(path)
    if normalised.startswith("MMSEQ"):
        return parse_mmseqs(path)
    if normalised in {"TSV", "EDGE", "EDGE-LIST"}:
        return parse_mmseqs(path)
    if normalised == "CSV":
        return parse_mmseqs(path, delimiter=",")
    if normalised.startswith("BLAST"):
        return parse_blast(path, options=options)
    raise FamilyInputError(f"Unsupported clustering format: {cluster_format}")


def _check_disjoint_clusters(clusters: Iterable[ParsedCluster], source: str):
    member_family: dict[str, str] = {}
    for cluster in clusters:
        for member in cluster.members:
            previous = member_family.get(member)
            if previous is not None and previous != cluster.label:
                raise FamilyInputError(
                    f"Sequence {member!r} occurs in multiple {source} clusters: {previous!r}, {cluster.label!r}"
                )
            member_family[member] = cluster.label


class _DisjointSet:
    def __init__(self, values):
        self.parent = {value: value for value in values}

    def add(self, value):
        self.parent.setdefault(value, value)

    def find(self, value):
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left, right):
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            keep, merge = sorted((left_root, right_root))
            self.parent[merge] = keep


@dataclass
class FamilyLedger:
    sequences: OrderedDict[str, str]
    families: OrderedDict[str, Family]
    resolver: TaxonResolver
    source_families: OrderedDict[str, ParsedCluster] = field(default_factory=OrderedDict)
    assignment: dict[str, str] = field(init=False, default_factory=dict)

    def __post_init__(self):
        for family_id, family in self.families.items():
            for member in family.members:
                previous = self.assignment.get(member)
                if previous is not None:
                    raise FamilyInputError(
                        f"Sequence {member!r} is assigned to final families {previous!r} and {family_id!r}"
                    )
                self.assignment[member] = family_id
        fasta_ids = set(self.sequences)
        assigned_ids = set(self.assignment)
        missing = sorted(fasta_ids - assigned_ids)
        unknown = sorted(assigned_ids - fasta_ids)
        if missing or unknown:
            details = []
            if missing:
                details.append(f"{len(missing)} unassigned FASTA IDs ({', '.join(missing[:10])})")
            if unknown:
                details.append(f"{len(unknown)} cluster IDs absent from FASTA ({', '.join(unknown[:10])})")
            raise FamilyInputError("Final family assignment is not exact: " + "; ".join(details))

    @property
    def taxa(self):
        return tuple(sorted({self.resolver.taxon(sequence_id) for sequence_id in self.sequences}))

    @property
    def genomes(self):
        return tuple(sorted({self.resolver.genome(sequence_id) for sequence_id in self.sequences}))

    def sorted_families(self, selected_ids=None):
        selected = set(selected_ids) if selected_ids is not None else None
        values = [family for family in self.families.values() if selected is None or family.family_id in selected]
        return sorted(values, key=lambda family: (-family.prevalence, -len(family.members), family.family_id))

    def choose_representative(self, members, preferred=()):
        candidates = [member for member in preferred if member in members]
        if not candidates:
            candidates = list(members)
        return min(candidates, key=lambda member: (-len(self.sequences[member]), member))

    def selected_single_copy(self, tolerance_percent=0.0):
        if tolerance_percent < 0 or tolerance_percent > 100:
            raise FamilyInputError("Single-copy tolerance must be between 0 and 100")
        selected = []
        for family in self.families.values():
            counts = family.copy_number
            occupied = len(counts)
            allowed = math.floor(occupied * tolerance_percent / 100.0)
            if sum(count > 1 for count in counts.values()) <= allowed:
                selected.append(family.family_id)
        return selected

    def write_assignments(self, path):
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "sequence_id", "family_id", "genome", "taxon", "representative",
                "stage", "category", "source_families", "source_clusters", "assignment_status",
            ])
            for sequence_id in sorted(self.sequences):
                family = self.families[self.assignment[sequence_id]]
                writer.writerow([
                    sequence_id,
                    family.family_id,
                    self.resolver.genome(sequence_id),
                    self.resolver.taxon(sequence_id),
                    family.representative,
                    family.stage,
                    family.category,
                    ";".join(family.source_families),
                    ";".join(family.source_clusters),
                    "assigned",
                ])


def build_family_ledger(
    fasta_path: str,
    cluster_path: str,
    cluster_format: str,
    mode="Species",
    recluster_path=None,
    sequence_tag="StORF",
    taxon_map=None,
    sequence_delimiter="|",
    options=None,
) -> FamilyLedger:
    sequences = read_fasta_records(fasta_path)
    resolver = TaxonResolver(mode=mode, mapping_path=taxon_map, sequence_delimiter=sequence_delimiter)
    parsed = parse_cluster_file(cluster_path, cluster_format, options=options)
    _check_disjoint_clusters(parsed, "initial")

    fasta_ids = set(sequences)
    clustered_ids = {member for cluster in parsed for member in cluster.members}
    unknown = sorted(clustered_ids - fasta_ids)
    if unknown:
        raise FamilyInputError(
            f"Initial clusters contain {len(unknown)} identifiers absent from FASTA: {', '.join(unknown[:10])}"
        )

    initial_clusters = list(parsed)
    for sequence_id in sorted(fasta_ids - clustered_ids):
        initial_clusters.append(ParsedCluster(f"singleton:{sequence_id}", (sequence_id,), sequence_id))

    source_families: OrderedDict[str, ParsedCluster] = OrderedDict()
    sequence_source: dict[str, str] = {}
    source_labels: dict[str, str] = {}
    source_is_first: dict[str, bool] = {}
    for parsed_cluster in initial_clusters:
        source_id = _stable_id("source", parsed_cluster.members)
        if source_id in source_families and source_families[source_id].members != parsed_cluster.members:
            raise FamilyInputError(f"Internal source-family identifier collision for {source_id}")
        source_families[source_id] = parsed_cluster
        source_labels[source_id] = parsed_cluster.label
        source_is_first[source_id] = not recluster_path or any(
            member in clustered_ids for member in parsed_cluster.members
        )
        for member in parsed_cluster.members:
            if member in sequence_source:
                raise FamilyInputError(f"Sequence {member!r} occurs in more than one initial family")
            sequence_source[member] = source_id

    dsu = _DisjointSet(source_families)
    recluster_by_root: dict[str, list[ParsedCluster]] = defaultdict(list)
    parsed_reclusters: list[ParsedCluster] = []
    if recluster_path:
        parsed_reclusters = parse_cluster_file(recluster_path, cluster_format, options=options)
        # A sequence may legitimately occur as a representative on every MMseqs row,
        # but it must not belong to distinct parsed clusters after adapter grouping.
        _check_disjoint_clusters(parsed_reclusters, "reclustered")
        recluster_unknown = sorted({m for c in parsed_reclusters for m in c.members} - fasta_ids)
        if recluster_unknown:
            raise FamilyInputError(
                f"Reclustered data contain {len(recluster_unknown)} identifiers absent from FASTA: "
                f"{', '.join(recluster_unknown[:10])}"
            )
        for recluster in parsed_reclusters:
            source_ids = list(OrderedDict.fromkeys(sequence_source[member] for member in recluster.members))
            if not source_ids:
                continue
            anchor = source_ids[0]
            for source_id in source_ids[1:]:
                dsu.union(anchor, source_id)

    components: dict[str, list[str]] = defaultdict(list)
    for source_id in source_families:
        components[dsu.find(source_id)].append(source_id)
    if parsed_reclusters:
        for recluster in parsed_reclusters:
            root = dsu.find(sequence_source[recluster.members[0]])
            recluster_by_root[root].append(recluster)

    built: list[Family] = []
    for root, source_ids in components.items():
        members = tuple(sorted(
            member
            for source_id in source_ids
            for member in source_families[source_id].members
        ))
        first_sources = [
            source_id for source_id in source_ids if source_is_first[source_id]
        ]
        has_second = bool(recluster_path) and len(first_sources) != len(source_ids)
        if not first_sources:
            stage = "only_Second"
        elif len(first_sources) > 1:
            stage = "combined"
        elif has_second:
            stage = "extended"
        else:
            stage = "First"

        reclusters = recluster_by_root.get(root, [])
        preferred = [cluster.representative for cluster in reclusters]
        preferred.extend(source_families[source_id].representative for source_id in source_ids)
        representative = min(
            [member for member in preferred if member in members] or list(members),
            key=lambda member: (-len(sequences[member]), member),
        )
        family_id = _stable_id("group", members)
        genomes = tuple(resolver.genome(member) for member in members)
        taxa = tuple(sorted(set(resolver.taxon(member) for member in members)))
        built.append(Family(
            family_id=family_id,
            members=members,
            representative=representative,
            taxa=taxa,
            genomes=genomes,
            stage=stage,
            source_families=tuple(sorted(source_ids)),
            source_clusters=tuple(sorted(source_labels[source_id] for source_id in source_ids)),
        ))

    built.sort(key=lambda family: family.family_id)
    families = OrderedDict((family.family_id, family) for family in built)
    return FamilyLedger(sequences, families, resolver, source_families)
