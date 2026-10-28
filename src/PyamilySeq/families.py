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
import logging
import re
import sqlite3
import tempfile
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Mapping


class FamilyInputError(ValueError):
    """Raised when clustering inputs cannot produce an unambiguous family ledger."""


def open_text(path, mode="rt"):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    return opener(path, mode, encoding="utf-8")


def _log_progress(logger, message: str, current: int | None = None, total: int | None = None):
    if logger is None:
        return
    if current is None:
        logger.info(message)
        return
    if total:
        percent = (100.0 * current / total) if total else 0.0
        logger.info("%s: %s (%0.1f%%)", message, f"{current:,}/{total:,}", percent)
    else:
        logger.info("%s: %s", message, f"{current:,}")


def fasta_identifier(defline: str) -> str:
    """Return the canonical FASTA identifier (the first whitespace token)."""

    identifier = defline.strip().split(maxsplit=1)[0] if defline.strip() else ""
    if not identifier:
        raise FamilyInputError("Encountered an empty FASTA header")
    return identifier


class SequenceStore:
    """Disk-backed FASTA store with bounded in-memory buffering.

    The store keeps only a small chunk of sequence records in RAM at a time and
    spills them into a local SQLite database. This makes very large FASTA files
    manageable without forcing the whole file into memory.
    """

    def __init__(self, path: str | None = None, buffer_bytes: int = 32 * 1024 * 1024, temp_dir: str | None = None):
        if path is None:
            fd, temp_path = tempfile.mkstemp(prefix="pyamilyseq_sequences_", suffix=".sqlite3", dir=temp_dir)
            os.close(fd)
            self.path = temp_path
            self._owns_path = True
        else:
            self.path = path
            self._owns_path = False
        self.buffer_bytes = max(1024 * 1024, int(buffer_bytes))
        self._conn = sqlite3.connect(self.path)
        self._conn.execute("PRAGMA journal_mode=OFF")
        self._conn.execute("PRAGMA synchronous=OFF")
        self._conn.execute("PRAGMA temp_store=MEMORY")
        self._conn.execute(
            "CREATE TABLE IF NOT EXISTS sequences (sequence_id TEXT PRIMARY KEY, sequence TEXT NOT NULL, length INTEGER NOT NULL)"
        )
        self._buffer: list[tuple[str, str, int]] = []
        self._buffer_size = 0
        self._count = 0
        self._closed = False

    @classmethod
    def from_existing_db(cls, path: str, logger=None) -> "SequenceStore":
        """Reopen a SQLite database previously produced by :meth:`from_fasta`.

        Used to resume after the (often very slow) FASTA-loading step has
        already completed, without re-parsing the original FASTA file.
        """
        if not os.path.exists(path):
            raise FamilyInputError(f"Sequence database not found: {path}")
        store = cls(path=path)
        has_rows = store._conn.execute("SELECT 1 FROM sequences LIMIT 1").fetchone() is not None
        if not has_rows:
            store.close()
            raise FamilyInputError(f"Sequence database at {path} contains no sequences")
        _log_progress(logger, f"Reusing existing sequence database {path}")
        return store

    @classmethod
    def from_fasta(
        cls,
        path: str,
        *,
        memory_limit_mb: int | float | None = None,
        buffer_cap_mb: int | float | None = None,
        temp_dir: str | None = None,
        logger=None,
        report_every: int = 100_000,
    ) -> "SequenceStore":
        if memory_limit_mb is None or float(memory_limit_mb) <= 0:
            memory_limit_mb = 128
        if buffer_cap_mb is None or float(buffer_cap_mb) <= 0:
            buffer_cap_mb = 256
        limit_bytes = int(float(memory_limit_mb) * 1024 * 1024)
        cap_bytes = int(float(buffer_cap_mb) * 1024 * 1024)
        buffer_bytes = max(8 * 1024 * 1024, min(limit_bytes // 10, cap_bytes))
        _log_progress(
            logger,
            f"Sequence buffering budget set to {buffer_bytes // (1024 * 1024):,} MB from -M={memory_limit_mb} "
            f"(cap {buffer_cap_mb} MB)",
        )
        store = cls(buffer_bytes=buffer_bytes, temp_dir=temp_dir)
        try:
            current = None
            chunks: list[str] = []
            count = 0
            _log_progress(logger, f"Loading FASTA records from {path}")
            with open_text(path) as handle:
                for line_number, line in enumerate(handle, start=1):
                    if line.startswith(">"):
                        if current is not None:
                            store.add(current, "".join(chunks))
                            count += 1
                            if report_every and count % report_every == 0:
                                _log_progress(logger, "  FASTA records loaded", count)
                        current = fasta_identifier(line[1:])
                        chunks = []
                    elif current is None:
                        if line.strip():
                            raise FamilyInputError(
                                f"Sequence data before the first FASTA header at line {line_number} in {path}"
                            )
                    else:
                        chunks.append("".join(line.split()))
            if current is not None:
                store.add(current, "".join(chunks))
                count += 1
            store.flush()
            if count == 0:
                raise FamilyInputError(f"No FASTA records found in {path}")
            _log_progress(logger, "Finished loading FASTA records", count)
            return store
        except Exception:
            store.close()
            if getattr(store, "_owns_path", False):
                try:
                    os.remove(store.path)
                except FileNotFoundError:
                    pass
            raise

    def add(self, sequence_id: str, sequence: str):
        if self._closed:
            raise FamilyInputError("Attempted to write to a closed SequenceStore")
        self._buffer.append((sequence_id, sequence, len(sequence)))
        self._buffer_size += len(sequence_id) + len(sequence)
        self._count += 1
        if self._buffer_size >= self.buffer_bytes:
            self.flush()

    def flush(self):
        if not self._buffer:
            return
        try:
            self._conn.executemany(
                "INSERT INTO sequences(sequence_id, sequence, length) VALUES (?, ?, ?)",
                self._buffer,
            )
        except sqlite3.IntegrityError as exc:
            raise FamilyInputError("Duplicate FASTA identifier encountered while loading sequences") from exc
        self._conn.commit()
        self._buffer.clear()
        self._buffer_size = 0

    def close(self):
        if self._closed:
            return
        try:
            self.flush()
        finally:
            self._conn.close()
            self._closed = True
            if getattr(self, "_owns_path", False):
                try:
                    os.remove(self.path)
                except FileNotFoundError:
                    pass

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __len__(self):
        row = self._conn.execute("SELECT COUNT(*) FROM sequences").fetchone()
        return int(row[0]) if row else 0

    def __contains__(self, sequence_id):
        row = self._conn.execute(
            "SELECT 1 FROM sequences WHERE sequence_id = ? LIMIT 1", (sequence_id,)
        ).fetchone()
        return row is not None

    def __getitem__(self, sequence_id):
        row = self._conn.execute(
            "SELECT sequence FROM sequences WHERE sequence_id = ?", (sequence_id,)
        ).fetchone()
        if row is None:
            raise KeyError(sequence_id)
        return row[0]

    def get(self, sequence_id, default=None):
        row = self._conn.execute(
            "SELECT sequence FROM sequences WHERE sequence_id = ?", (sequence_id,)
        ).fetchone()
        return row[0] if row is not None else default

    def sequence_length(self, sequence_id):
        row = self._conn.execute(
            "SELECT length FROM sequences WHERE sequence_id = ?", (sequence_id,)
        ).fetchone()
        if row is None:
            raise KeyError(sequence_id)
        return int(row[0])

    def keys(self):
        for (sequence_id,) in self._conn.execute("SELECT sequence_id FROM sequences"):
            yield sequence_id

    def items(self):
        for sequence_id, sequence in self._conn.execute("SELECT sequence_id, sequence FROM sequences"):
            yield sequence_id, sequence

    def values(self):
        for _, sequence in self.items():
            yield sequence

    def __iter__(self):
        return self.keys()


def read_fasta_records(path: str, logger=None, report_every: int = 100_000) -> OrderedDict[str, str]:
    """Read FASTA records into memory.

    Use :class:`SequenceStore.from_fasta` when memory pressure matters.
    """

    store = SequenceStore.from_fasta(path, logger=logger, report_every=report_every)
    try:
        records = OrderedDict(store.items())
    finally:
        store.close()
    if not records:
        raise FamilyInputError(f"No FASTA records found in {path}")
    empty = [identifier for identifier, sequence in records.items() if not sequence]
    if empty:
        sample = ", ".join(empty[:10])
        raise FamilyInputError(f"FASTA contains {len(empty)} empty sequence(s): {sample}")
    return records

    records: OrderedDict[str, str] = OrderedDict()
    current = None
    chunks: list[str] = []
    count = 0
    _log_progress(logger, f"Loading FASTA records from {path}")
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks)
                    count += 1
                    if report_every and count % report_every == 0:
                        _log_progress(logger, "  FASTA records loaded", count)
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
        count += 1
    if not records:
        raise FamilyInputError(f"No FASTA records found in {path}")
    empty = [identifier for identifier, sequence in records.items() if not sequence]
    if empty:
        sample = ", ".join(empty[:10])
        raise FamilyInputError(f"FASTA contains {len(empty)} empty sequence(s): {sample}")
    _log_progress(logger, "Finished loading FASTA records", count)
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


def iter_cdhit_clusters(path: str):
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
        return ParsedCluster(label, unique, rep)

    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                finished = finish()
                if finished is not None:
                    yield finished
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
    finished = finish()
    if finished is not None:
        yield finished


def parse_cdhit(path: str) -> list[ParsedCluster]:
    clusters = list(iter_cdhit_clusters(path))
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


def iter_mmseqs_clusters(path: str, delimiter="\t", logger=None, report_every: int = 100_000):
    grouped: OrderedDict[str, list[str]] = OrderedDict()
    seen: dict[str, set[str]] = {}
    row_count = 0
    _log_progress(logger, f"Parsing MMseqs clusters from {path}")
    for _, row in _edge_rows(path, delimiter=delimiter):
        representative, member = row[:2]
        members = grouped.setdefault(representative, [])
        member_seen = seen.setdefault(representative, set())
        if representative not in member_seen:
            members.append(representative)
            member_seen.add(representative)
        if member not in member_seen:
            members.append(member)
            member_seen.add(member)
        row_count += 1
        if report_every and row_count % report_every == 0:
            _log_progress(logger, "  MMseqs rows parsed", row_count)
    if not grouped:
        raise FamilyInputError(f"No MMseqs clusters found in {path}")
    _log_progress(logger, "Finished parsing MMseqs clusters", row_count)
    for rep, members in grouped.items():
        yield ParsedCluster(rep, tuple(members), rep)


def parse_mmseqs(path: str, delimiter="\t", logger=None, report_every: int = 100_000) -> list[ParsedCluster]:
    return list(iter_mmseqs_clusters(path, delimiter=delimiter, logger=logger, report_every=report_every))


def _float_cell(row, index, name, line_number, path):
    try:
        return float(row[index])
    except (IndexError, ValueError):
        raise FamilyInputError(
            f"BLAST {name} is missing or invalid at line {line_number} in {path}"
        ) from None


def iter_blast_clusters(path: str, options=None, delimiter="\t", logger=None, report_every: int = 100_000):
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
    scored = False
    best_score: dict[str, float] = {}
    best_subjects: dict[str, set[str]] = defaultdict(set)
    row_count = 0
    for line_number, row in _edge_rows(path, delimiter=delimiter):
        query, subject = row[:2]
        nodes.update((query, subject))
        if len(row) > 2:
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
        else:
            score = 1.0

        row_count += 1
        if scored and graph_mode == "reciprocal-best":
            if query == subject:
                continue
            if query not in best_score or score > best_score[query]:
                best_score[query] = score
                best_subjects[query] = {subject}
            elif score == best_score[query]:
                best_subjects[query].add(subject)
        else:
            if query != subject:
                best_subjects[query].add(subject)
                best_subjects[subject].add(query)
        if report_every and row_count % report_every == 0:
            _log_progress(logger, "  BLAST rows parsed", row_count)

    if not nodes:
        raise FamilyInputError(f"No BLAST nodes found in {path}")
    adjacency: dict[str, set[str]] = {node: set() for node in nodes}
    if scored and graph_mode == "reciprocal-best":
        for query, subjects in best_subjects.items():
            for subject in subjects:
                if query != subject and query in best_subjects.get(subject, ()) and subject in best_subjects.get(query, ()):
                    adjacency[query].add(subject)
                    adjacency[subject].add(query)
    else:
        for query, subjects in best_subjects.items():
            for subject in subjects:
                if query != subject:
                    adjacency[query].add(subject)
                    adjacency[subject].add(query)

    clusters: list[ParsedCluster] = []
    visited: set[str] = set()
    weighted_degree = Counter({node: len(neighbors) for node, neighbors in adjacency.items()})
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
    _log_progress(logger, f"Finished parsing BLAST clusters: {len(clusters):,} components from {row_count:,} rows")
    for cluster in clusters:
        yield cluster


def parse_blast(path: str, options=None, delimiter="\t", logger=None, report_every: int = 100_000) -> list[ParsedCluster]:
    return list(
        iter_blast_clusters(
            path,
            options=options,
            delimiter=delimiter,
            logger=logger,
            report_every=report_every,
        )
    )


def iter_cluster_file(path: str, cluster_format: str, options=None, logger=None):
    normalised = cluster_format.strip().upper().replace("_", "-")
    if normalised in {"CD-HIT", "CDHIT", "CD-HIT-EST"}:
        yield from iter_cdhit_clusters(path)
        return
    if normalised.startswith("MMSEQ"):
        yield from iter_mmseqs_clusters(path, logger=logger)
        return
    if normalised in {"TSV", "EDGE", "EDGE-LIST"}:
        yield from iter_mmseqs_clusters(path, logger=logger)
        return
    if normalised == "CSV":
        yield from iter_mmseqs_clusters(path, delimiter=",", logger=logger)
        return
    if normalised.startswith("BLAST"):
        yield from iter_blast_clusters(path, options=options, logger=logger)
        return
    raise FamilyInputError(f"Unsupported clustering format: {cluster_format}")


def parse_cluster_file(path: str, cluster_format: str, options=None, logger=None) -> list[ParsedCluster]:
    return list(iter_cluster_file(path, cluster_format, options=options, logger=logger))


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
    sequences: SequenceStore
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
        missing = []
        for sequence_id in self.sequences:
            if sequence_id not in self.assignment:
                missing.append(sequence_id)
                if len(missing) >= 10:
                    break
        unknown = []
        for sequence_id in self.assignment:
            if sequence_id not in self.sequences:
                unknown.append(sequence_id)
                if len(unknown) >= 10:
                    break
        if missing or unknown:
            details = []
            if missing:
                details.append(f"unassigned FASTA IDs ({', '.join(missing)})")
            if unknown:
                details.append(f"cluster IDs absent from FASTA ({', '.join(unknown)})")
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
            # Streamed in the sequence store's natural (physical) order rather than
            # sorted() - at billions of records, Python's sorted() would first have
            # to materialize every sequence_id into a list in memory before sorting,
            # which is the same order-of-magnitude memory cost we hit earlier with
            # the in-memory ID cache. If callers need alphabetical order, `sort` on
            # the resulting TSV is a far cheaper external merge-sort than doing it
            # here in Python.
            for sequence_id in self.sequences:
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
    logger=None,
    sequence_db_path=None,
) -> FamilyLedger:
    logger = logger or logging.getLogger("PyamilySeq.PyamilySeq")
    _log_progress(logger, "Building canonical family ledger")
    sequence_db_path = sequence_db_path or getattr(options, "sequence_db", None)
    if sequence_db_path:
        sequences = SequenceStore.from_existing_db(sequence_db_path, logger=logger)
    else:
        sequences = SequenceStore.from_fasta(
            fasta_path,
            memory_limit_mb=getattr(options, "mem", None),
            buffer_cap_mb=getattr(options, "sequence_buffer_cap_mb", None),
            temp_dir=getattr(options, "tmp_dir", None),
            logger=logger,
        )
    resolver = TaxonResolver(mode=mode, mapping_path=taxon_map, sequence_delimiter=sequence_delimiter)
    _log_progress(logger, f"Parsing initial {cluster_format} clusters from {cluster_path}")
    source_families: OrderedDict[str, ParsedCluster] = OrderedDict()
    sequence_source: dict[str, str] = {}
    source_labels: dict[str, str] = {}
    source_is_first: dict[str, bool] = {}
    parsed_count = 0
    for parsed_cluster in iter_cluster_file(cluster_path, cluster_format, options=options, logger=logger):
        parsed_count += 1
        source_id = _stable_id("source", parsed_cluster.members)
        source_families[source_id] = parsed_cluster
        source_labels[source_id] = parsed_cluster.label
        source_is_first[source_id] = True
        for member in parsed_cluster.members:
            if member not in sequences:
                raise FamilyInputError(
                    f"Initial clusters contain identifier absent from FASTA: {member}. "
                    "The initial cluster file and -original_fasta must describe the same sequence set; "
                    "a representative-only FASTA cannot be paired with the full original clustering."
                )
            if member in sequence_source:
                raise FamilyInputError(f"Sequence {member!r} occurs in more than one initial family")
            sequence_source[member] = source_id
    _log_progress(logger, f"Parsed {parsed_count:,} initial clusters")

    singleton_count = sum(1 for sequence_id in sequences if sequence_id not in sequence_source)
    if singleton_count:
        _log_progress(logger, f"Backfilling {singleton_count:,} missing FASTA singletons")
    for sequence_id in sequences:
        if sequence_id in sequence_source:
            continue
        source_id = _stable_id("source", (sequence_id,))
        parsed_cluster = ParsedCluster(f"singleton:{sequence_id}", (sequence_id,), sequence_id)
        source_families[source_id] = parsed_cluster
        source_labels[source_id] = parsed_cluster.label
        source_is_first[source_id] = not recluster_path
        sequence_source[sequence_id] = source_id

    dsu = _DisjointSet(source_families)
    recluster_by_root: dict[str, list[ParsedCluster]] = defaultdict(list)
    parsed_reclusters: list[ParsedCluster] = []
    if recluster_path:
        _log_progress(logger, f"Parsing reclustered families from {recluster_path}")
        parsed_reclusters = list(iter_cluster_file(recluster_path, cluster_format, options=options, logger=logger))
        # A sequence may legitimately occur as a representative on every MMseqs row,
        # but it must not belong to distinct parsed clusters after adapter grouping.
        recluster_unknown = sorted(
            {m for c in parsed_reclusters for m in c.members if m not in sequences}
        )
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
    _log_progress(logger, f"Built {len(families):,} final families")
    return FamilyLedger(sequences, families, resolver, source_families)
