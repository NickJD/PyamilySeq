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
from datetime import datetime, timezone
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass, field
from typing import Iterable, Mapping
from urllib.parse import quote


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
    """Disk-backed FASTA store with bounded buffering and optional persistence.

    Temporary stores are deleted when closed. Stores created with an explicit
    ``database_path`` are built under a staging filename, validated, atomically
    promoted to the requested path, and retained for later reuse.
    """

    DB_SCHEMA_VERSION = 1
    METADATA_TABLE = "pyamilyseq_metadata"

    def __init__(
        self,
        path: str | None = None,
        buffer_bytes: int = 32 * 1024 * 1024,
        temp_dir: str | None = None,
        *,
        initialize_schema: bool = True,
    ):
        if path is None:
            fd, temp_path = tempfile.mkstemp(
                prefix="pyamilyseq_sequences_",
                suffix=".sqlite3",
                dir=temp_dir,
            )
            os.close(fd)
            self.path = temp_path
            self._owns_path = True
        else:
            self.path = os.path.abspath(path)
            self._owns_path = False

        self.persistent = not self._owns_path
        self.buffer_bytes = max(1024 * 1024, int(buffer_bytes))
        self._buffer: list[tuple[str, str, int]] = []
        self._buffer_size = 0
        self._count = 0
        self._closed = False
        self._conn = sqlite3.connect(self.path, timeout=60.0)
        try:
            # WAL gives transactions recoverable commit semantics. Use FULL for
            # retained checkpoints; temporary stores can use NORMAL because the
            # source FASTA remains the recoverable source of truth.
            self._conn.execute("PRAGMA journal_mode=WAL")
            self._conn.execute(
                "PRAGMA synchronous=FULL" if self.persistent else "PRAGMA synchronous=NORMAL"
            )
            self._conn.execute("PRAGMA temp_store=MEMORY")
            self._conn.execute("PRAGMA busy_timeout=60000")
            # Negative cache_size values are kibibytes; this sets a 64 MiB page cache.
            self._conn.execute("PRAGMA cache_size=-65536")

            if initialize_schema:
                self._conn.execute(
                    "CREATE TABLE IF NOT EXISTS sequences ("
                    "sequence_id TEXT PRIMARY KEY, sequence TEXT NOT NULL, length INTEGER NOT NULL)"
                )
                self._conn.commit()
        except Exception:
            self._conn.close()
            self._closed = True
            raise

    @staticmethod
    def _remove_database_files(path: str) -> None:
        """Remove a database and SQLite sidecar files, if present."""
        for candidate in (path, path + "-wal", path + "-shm", path + "-journal"):
            try:
                os.remove(candidate)
            except FileNotFoundError:
                pass

    def _metadata_table_exists(self) -> bool:
        row = self._conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
            (self.METADATA_TABLE,),
        ).fetchone()
        return row is not None

    def _write_metadata(self, values: Mapping[str, object]) -> None:
        self._conn.execute(
            f"CREATE TABLE IF NOT EXISTS {self.METADATA_TABLE} ("
            "key TEXT PRIMARY KEY, value TEXT NOT NULL)"
        )
        self._conn.executemany(
            f"INSERT OR REPLACE INTO {self.METADATA_TABLE}(key, value) VALUES (?, ?)",
            [(str(key), str(value)) for key, value in values.items()],
        )
        self._conn.commit()

    def _read_metadata(self) -> dict[str, str]:
        if not self._metadata_table_exists():
            return {}
        return dict(
            self._conn.execute(
                f"SELECT key, value FROM {self.METADATA_TABLE}"
            ).fetchall()
        )

    def validate(self, expected_count=None, full_check=False):
        """Validate SQLite structure, schema and optionally stored lengths."""
        required_columns = {"sequence_id", "sequence", "length"}
        actual_columns = {
            row[1] for row in self._conn.execute("PRAGMA table_info(sequences)")
        }
        if not required_columns.issubset(actual_columns):
            raise FamilyInputError(
                "SQLite sequence database has an incompatible schema; "
                f"required columns are {sorted(required_columns)}, found {sorted(actual_columns)}"
            )

        check = "integrity_check" if full_check else "quick_check"
        check_rows = self._conn.execute(f"PRAGMA {check}").fetchall()
        check_errors = [str(row[0]) for row in check_rows if not row or row[0] != "ok"]
        if not check_rows or check_errors:
            detail = "; ".join(check_errors) if check_errors else "no result returned"
            raise FamilyInputError(
                f"SQLite database validation failed ({check}): {detail}"
            )

        actual_count = len(self)
        if expected_count is not None and actual_count != int(expected_count):
            raise FamilyInputError(
                f"SQLite sequence count mismatch: expected {int(expected_count):,}, "
                f"found {actual_count:,}"
            )

        if full_check:
            mismatch = self._conn.execute(
                "SELECT COUNT(*) FROM sequences WHERE length != length(sequence)"
            ).fetchone()[0]
            if mismatch:
                raise FamilyInputError(
                    f"SQLite database contains {mismatch:,} sequence record(s) whose stored "
                    "length does not match the sequence value"
                )

        return actual_count

    def checkpoint(self, truncate=True) -> None:
        """Flush pending records and checkpoint the WAL into the main database."""
        if self._closed:
            raise FamilyInputError("Cannot checkpoint a closed SequenceStore")
        self.flush()
        mode = "TRUNCATE" if truncate else "PASSIVE"
        result = self._conn.execute(f"PRAGMA wal_checkpoint({mode})").fetchone()
        if result and int(result[0]) != 0:
            raise FamilyInputError(
                "SQLite WAL checkpoint could not complete because the database is busy: "
                f"{result}"
            )

    @classmethod
    def from_existing_db(
        cls,
        path: str,
        logger=None,
        *,
        validate_integrity: bool = True,
    ) -> "SequenceStore":
        """Open and validate a retained sequence database.

        Databases created before metadata/checkpoint support are accepted as
        legacy stores after structural validation, but are explicitly reported
        as such. New databases must have a complete status marker.
        """
        path = os.path.abspath(path)
        if not os.path.isfile(path):
            raise FamilyInputError(f"Sequence database not found: {path}")

        # Preflight the schema with a read-only SQLite connection. This prevents
        # an arbitrary/invalid file from being modified by journal-mode changes
        # before we establish that it is a PyamilySeq sequence database.
        preflight = None
        try:
            uri = "file:" + quote(path, safe="/") + "?mode=ro"
            preflight = sqlite3.connect(uri, uri=True, timeout=60.0)
            table = preflight.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='sequences'"
            ).fetchone()
            if table is None:
                raise FamilyInputError(
                    f"File is not a PyamilySeq sequence database (missing sequences table): {path}"
                )
            columns = {
                row[1] for row in preflight.execute("PRAGMA table_info(sequences)")
            }
            required = {"sequence_id", "sequence", "length"}
            if not required.issubset(columns):
                raise FamilyInputError(
                    f"SQLite sequence database has an incompatible schema: {path}; "
                    f"required columns are {sorted(required)}, found {sorted(columns)}"
                )
        except sqlite3.DatabaseError as exc:
            raise FamilyInputError(
                f"Could not open or inspect SQLite sequence database {path}: {exc}"
            ) from exc
        finally:
            if preflight is not None:
                preflight.close()

        store = cls(path=path, initialize_schema=False)
        try:
            sequences_table = store._conn.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='sequences'"
            ).fetchone()
            if sequences_table is None:
                raise FamilyInputError(
                    f"File is not a PyamilySeq sequence database (missing sequences table): {path}"
                )

            has_metadata = store._metadata_table_exists()
            expected_count = None
            if has_metadata:
                metadata = store._read_metadata()
                schema_version = metadata.get("schema_version")
                status = metadata.get("status")

                if schema_version != str(cls.DB_SCHEMA_VERSION):
                    raise FamilyInputError(
                        f"Unsupported PyamilySeq sequence database schema version: {schema_version!r}"
                    )
                if status != "complete":
                    raise FamilyInputError(
                        f"Sequence database is not a completed checkpoint (status={status!r}): {path}"
                    )
                try:
                    expected_count = int(metadata["sequence_count"])
                except (KeyError, TypeError, ValueError):
                    raise FamilyInputError(
                        f"Sequence database metadata is missing a valid sequence_count: {path}"
                    ) from None
            else:
                _log_progress(
                    logger,
                    "WARNING: Reusing a legacy SQLite sequence database without completion metadata; "
                    "structural checks will be performed, but import completeness cannot be verified.",
                )

            if validate_integrity:
                count = store.validate(expected_count=expected_count)
            else:
                # Used only immediately after from_fasta() has already checked
                # the staged database before atomic promotion.
                count = len(store)
                if expected_count is not None and count != expected_count:
                    raise FamilyInputError(
                        f"SQLite sequence count mismatch: expected {expected_count:,}, found {count:,}"
                    )
            if count == 0:
                raise FamilyInputError(f"Sequence database at {path} contains no sequences")
            _log_progress(logger, f"Validated SQLite sequence database: {count:,} records")
            _log_progress(logger, f"Reusing existing sequence database {path}")
            return store
        except Exception:
            try:
                store.close()
            except Exception:
                pass
            raise

    @classmethod
    def from_fasta(
        cls,
        path: str,
        *,
        memory_limit_mb: int | float | None = None,
        buffer_cap_mb: int | float | None = None,
        temp_dir: str | None = None,
        database_path: str | None = None,
        logger=None,
        report_every: int = 100_000,
    ) -> "SequenceStore":
        """Build a sequence store from FASTA.

        If ``database_path`` is omitted, a temporary database is created and
        deleted on close. If supplied, the database is staged in the target
        directory, marked complete only after validation, and retained there.
        """
        if memory_limit_mb is None or float(memory_limit_mb) <= 0:
            memory_limit_mb = 128
        if buffer_cap_mb is None or float(buffer_cap_mb) <= 0:
            buffer_cap_mb = 256

        limit_bytes = int(float(memory_limit_mb) * 1024 * 1024)
        cap_bytes = int(float(buffer_cap_mb) * 1024 * 1024)
        buffer_bytes = max(8 * 1024 * 1024, min(limit_bytes // 10, cap_bytes))
        _log_progress(
            logger,
            f"Sequence buffering budget set to {buffer_bytes // (1024 * 1024):,} MB "
            f"from -M={memory_limit_mb} (cap {buffer_cap_mb} MB)",
        )

        source_path = os.path.abspath(path)
        final_path = os.path.abspath(database_path) if database_path else None
        if final_path and source_path == final_path:
            raise FamilyInputError(
                "The sequence database output path cannot be the same as the input FASTA path"
            )

        staging_path = None
        if final_path:
            parent = os.path.dirname(final_path) or os.getcwd()
            os.makedirs(parent, exist_ok=True)
            if any(os.path.exists(candidate) for candidate in (
                final_path, final_path + "-wal", final_path + "-shm", final_path + "-journal"
            )):
                raise FamilyInputError(
                    f"Sequence database output already exists (or has SQLite sidecars): {final_path}. "
                    "Choose a new path; existing checkpoints are never overwritten."
                )
            fd, staging_path = tempfile.mkstemp(
                prefix=f".{os.path.basename(final_path)}.",
                suffix=".building",
                dir=parent,
            )
            os.close(fd)
            store = cls(path=staging_path, buffer_bytes=buffer_bytes)
        else:
            store = cls(buffer_bytes=buffer_bytes, temp_dir=temp_dir)

        if final_path:
            _log_progress(logger, f"Building retained SQLite sequence checkpoint at {final_path}")
        else:
            _log_progress(logger, f"Using temporary SQLite sequence store: {store.path}")

        try:
            try:
                source_stat = os.stat(source_path)
                source_metadata = {
                    "source_fasta": source_path,
                    "source_size_bytes": source_stat.st_size,
                    "source_mtime_ns": source_stat.st_mtime_ns,
                }
            except OSError:
                source_metadata = {"source_fasta": source_path}

            store._write_metadata({
                "schema_version": cls.DB_SCHEMA_VERSION,
                "status": "building",
                "created_at_utc": datetime.now(timezone.utc).isoformat(),
                **source_metadata,
            })

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

            store._write_metadata({
                "status": "complete",
                "sequence_count": count,
                "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            })
            # One structural pass is sufficient: the database is still staged
            # and is not promoted to the requested checkpoint path until this
            # validation succeeds.
            store.validate(expected_count=count)
            _log_progress(logger, "SQLite sequence database validated", count)
            _log_progress(logger, "Finished loading FASTA records", count)

            if not final_path:
                return store

            # Close and checkpoint the staged database before making it visible
            # under the requested final filename.
            store.checkpoint(truncate=True)
            store.close()
            try:
                # A hard link provides atomic, no-overwrite promotion on the
                # same filesystem; staging is deliberately in the target dir.
                os.link(staging_path, final_path)
                os.unlink(staging_path)
            except FileExistsError:
                raise FamilyInputError(
                    f"Sequence database output appeared during creation: {final_path}"
                ) from None
            except OSError as exc:
                raise FamilyInputError(
                    f"Could not promote the completed sequence database to {final_path}: {exc}. "
                    "The target filesystem must support hard links for atomic checkpoint creation."
                ) from exc

            staging_path = None
            retained_store = cls.from_existing_db(
                final_path, logger=logger, validate_integrity=False
            )
            _log_progress(logger, f"Retained sequence database checkpoint at {final_path}")
            return retained_store

        except Exception:
            try:
                store.close()
            except Exception:
                pass
            if staging_path:
                cls._remove_database_files(staging_path)
            raise

    def add(self, sequence_id: str, sequence: str):
        if self._closed:
            raise FamilyInputError("Attempted to write to a closed SequenceStore")
        try:
            sequence.encode("ascii")
        except UnicodeEncodeError as exc:
            raise FamilyInputError(
                f"Sequence {sequence_id!r} contains non-ASCII characters; biological sequences must be ASCII"
            ) from exc
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
            self._conn.commit()
        except sqlite3.IntegrityError as exc:
            raise FamilyInputError(
                "Duplicate FASTA identifier encountered while loading sequences"
            ) from exc
        self._buffer.clear()
        self._buffer_size = 0

    def close(self):
        if self._closed:
            return
        close_error = None
        try:
            self.flush()
            try:
                result = self._conn.execute("PRAGMA wal_checkpoint(TRUNCATE)").fetchone()
                if result and int(result[0]) != 0:
                    raise FamilyInputError(
                        f"SQLite WAL checkpoint could not complete: {result}"
                    )
            except Exception as exc:
                close_error = exc
        except Exception as exc:
            close_error = exc
        finally:
            try:
                self._conn.close()
            finally:
                self._closed = True
                if self._owns_path:
                    self._remove_database_files(self.path)
        if close_error is not None:
            raise close_error

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    @staticmethod
    def _decode_sequence(sequence_id, payload):
        if payload is None:
            raise KeyError(sequence_id)
        if isinstance(payload, str):
            return payload
        try:
            # CAST(sequence AS BLOB) avoids sqlite3's implicit UTF-8 decoding.
            # Protein and nucleotide sequence alphabets are ASCII.
            return bytes(payload).decode("ascii")
        except UnicodeDecodeError as exc:
            raise FamilyInputError(
                f"Stored sequence for {sequence_id!r} contains non-ASCII/corrupt bytes in SQLite"
            ) from exc

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
            "SELECT CAST(sequence AS BLOB) FROM sequences WHERE sequence_id = ?",
            (sequence_id,),
        ).fetchone()
        if row is None:
            raise KeyError(sequence_id)
        return self._decode_sequence(sequence_id, row[0])

    def get(self, sequence_id, default=None):
        row = self._conn.execute(
            "SELECT CAST(sequence AS BLOB) FROM sequences WHERE sequence_id = ?",
            (sequence_id,),
        ).fetchone()
        return self._decode_sequence(sequence_id, row[0]) if row is not None else default

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
        for sequence_id, payload in self._conn.execute(
            "SELECT sequence_id, CAST(sequence AS BLOB) FROM sequences"
        ):
            yield sequence_id, self._decode_sequence(sequence_id, payload)

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
        return min(candidates, key=lambda member: (-self.sequences.sequence_length(member), member))

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
    supplemental_sequence_ids=None,
    sequence_db_out_path=None,
) -> FamilyLedger:
    logger = logger or logging.getLogger("PyamilySeq.PyamilySeq")
    _log_progress(logger, "Building canonical family ledger")
    sequence_db_path = sequence_db_path or getattr(options, "sequence_db", None)
    sequence_db_out_path = sequence_db_out_path or getattr(options, "sequence_db_out", None)
    if sequence_db_path and sequence_db_out_path:
        raise FamilyInputError(
            "Cannot use -sequence_db and --sequence-db-out together; choose either an existing database "
            "or a FASTA input from which to build a retained database."
        )
    if sequence_db_path:
        sequences = SequenceStore.from_existing_db(sequence_db_path, logger=logger)
    else:
        sequences = SequenceStore.from_fasta(
            fasta_path,
            memory_limit_mb=getattr(options, "mem", None),
            buffer_cap_mb=getattr(options, "sequence_buffer_cap_mb", None),
            temp_dir=getattr(options, "tmp_dir", None),
            database_path=sequence_db_out_path,
            logger=logger,
        )
    try:
        if supplemental_sequence_ids is not None:
            supplemental_sequence_ids = set(supplemental_sequence_ids)
            missing_supplemental = [
                sequence_id for sequence_id in supplemental_sequence_ids
                if sequence_id not in sequences
            ]
            if missing_supplemental:
                examples = ", ".join(sorted(missing_supplemental)[:10])
                raise FamilyInputError(
                    f"{len(missing_supplemental):,} declared supplementary sequence ID(s) are absent "
                    f"from the combined FASTA: {examples}"
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

        if supplemental_sequence_ids:
            overlap = supplemental_sequence_ids.intersection(sequence_source)
            if overlap:
                examples = ", ".join(sorted(overlap)[:10])
                raise FamilyInputError(
                    f"{len(overlap):,} sequence ID(s) declared supplementary are already present "
                    f"in the original cluster file: {examples}. Check new-gene ID normalisation."
                )

        singleton_count = sum(1 for sequence_id in sequences if sequence_id not in sequence_source)
        if singleton_count:
            _log_progress(logger, f"Backfilling {singleton_count:,} missing FASTA singletons")
            if supplemental_sequence_ids is not None:
                supplement_singletons = sum(
                    1 for sequence_id in sequences
                    if sequence_id not in sequence_source and sequence_id in supplemental_sequence_ids
                )
                original_singletons = singleton_count - supplement_singletons
                _log_progress(
                    logger,
                    f"Backfilled singleton provenance: {original_singletons:,} original; "
                    f"{supplement_singletons:,} supplementary",
                )
        for sequence_id in sequences:
            if sequence_id in sequence_source:
                continue

            source_id = _stable_id("source", (sequence_id,))
            parsed_cluster = ParsedCluster(
                f"singleton:{sequence_id}",
                (sequence_id,),
                sequence_id,
            )

            source_families[source_id] = parsed_cluster
            source_labels[source_id] = parsed_cluster.label

            if supplemental_sequence_ids is None:
                # Preserve existing behaviour for ordinary Partial/reclustered runs.
                source_is_first[source_id] = not recluster_path
            else:
                # In Supplementation, explicitly identify whether this sequence
                # originated from the original dataset or the new supplement.
                source_is_first[source_id] = (
                        sequence_id not in supplemental_sequence_ids
                )

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
                key=lambda member: (-sequences.sequence_length(member), member),
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
    except BaseException:
        # If ledger construction fails after the SQLite store has been opened,
        # close it deterministically. Temporary stores are removed; retained
        # checkpoints are closed but deliberately preserved for diagnosis/reuse.
        try:
            sequences.close()
        except Exception:
            logger.exception("Failed to close SQLite sequence store after ledger-construction error")
        raise