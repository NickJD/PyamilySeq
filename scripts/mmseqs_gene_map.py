#!/usr/bin/env python3
"""Build a PyamilySeq-compatible sequence_id -> family_id map directly from an
MMseqs2 cluster TSV, bypassing PyamilySeq Partial's SequenceStore/SQLite
staging entirely.

This shortcut is only valid when you are NOT merging in a second-pass
re-clustered file (PyamilySeq's `-reclustered` option) - without that step,
PyamilySeq's final families are exactly the raw MMseqs clusters plus a
singleton family for any FASTA sequence absent from the cluster file, so
family membership never depends on sequence content or length. That means
no FASTA sequences ever need to be loaded, staged, or scanned by anything
other than a cheap one-pass header read here.

Output is a two-column TSV (sequence_id, family_id) matching what
pyamilyseq_geneorder_representation.py's --gene-map expects.
"""

import argparse
import csv
import gzip
import hashlib
import sys
from collections import OrderedDict


def open_text(path, mode="rt"):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    return opener(path, mode, encoding="utf-8")


def stable_id(prefix: str, members) -> str:
    # Same convention as PyamilySeq's families.py _stable_id, for consistency
    # with any full Partial-mode run you may compare against.
    payload = "\0".join(sorted(set(members))).encode("utf-8")
    return f"{prefix}_{hashlib.sha256(payload).hexdigest()[:12]}"


def iter_mmseqs_clusters(path, delimiter="\t", report_every=1_000_000):
    grouped: "OrderedDict[str, list[str]]" = OrderedDict()
    seen: dict[str, set[str]] = {}
    row_count = 0
    print(f"Parsing MMseqs clusters from {path}", file=sys.stderr)
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        for row in reader:
            if not row or not row[0].strip() or row[0].startswith("#"):
                continue
            if len(row) < 2 or not row[1].strip():
                raise ValueError(f"Expected at least two columns, got: {row!r}")
            representative, member = row[0].strip(), row[1].strip()
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
                print(f"  MMseqs rows parsed: {row_count:,}", file=sys.stderr)
    if not grouped:
        raise ValueError(f"No MMseqs clusters found in {path}")
    print(f"Finished parsing MMseqs clusters: {row_count:,} rows, {len(grouped):,} clusters", file=sys.stderr)
    return grouped


def iter_fasta_ids(path):
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].strip().split(maxsplit=1)[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cluster-file", required=True, help="MMseqs2 cluster TSV (.tsv or .tsv.gz)")
    parser.add_argument("--fasta", required=True, help="Original FASTA - headers only are read, never sequence content")
    parser.add_argument("--output", required=True, help="Path to write the sequence_id<TAB>family_id TSV")
    parser.add_argument("--delimiter", default="\t")
    args = parser.parse_args()

    grouped = iter_mmseqs_clusters(args.cluster_file, delimiter=args.delimiter)

    assignment: dict[str, str] = {}
    for representative, members in grouped.items():
        family_id = stable_id("group", members)
        for member in members:
            if member in assignment:
                sys.exit(f"Sequence {member!r} occurs in more than one MMseqs cluster")
            assignment[member] = family_id
    print(f"Assigned {len(assignment):,} sequences across {len(grouped):,} MMseqs clusters", file=sys.stderr)

    singleton_count = 0
    seen_fasta: set[str] = set()
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sequence_id", "family_id"])
        for sequence_id in iter_fasta_ids(args.fasta):
            if sequence_id in seen_fasta:
                sys.exit(f"Duplicate FASTA identifier: {sequence_id!r}")
            seen_fasta.add(sequence_id)
            family_id = assignment.get(sequence_id)
            if family_id is None:
                family_id = stable_id("group", (sequence_id,))
                singleton_count += 1
            writer.writerow([sequence_id, family_id])

    unknown = set(assignment) - seen_fasta
    if unknown:
        sample = ", ".join(sorted(unknown)[:5])
        sys.exit(
            f"{len(unknown):,} MMseqs cluster identifiers are absent from the FASTA "
            f"(sample: {sample})"
        )

    print(f"Backfilled {singleton_count:,} FASTA-only singletons", file=sys.stderr)
    print(f"Wrote {len(seen_fasta):,} sequence_id -> family_id rows to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()