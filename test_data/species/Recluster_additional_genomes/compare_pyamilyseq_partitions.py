#!/usr/bin/env python3
"""Compare PyamilySeq gene-to-family partitions from two runs.

PyamilySeq gene_group_table.tsv is expected to contain two tab-separated
columns per line, without a header: gene_id<TAB>family_id.

The comparison ignores family IDs and compares the actual sets of gene IDs
assigned to each family. It also writes an overlap table useful for diagnosing
family splits and merges.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def load_gene_family_map(path: Path):
    gene_to_family: dict[str, str] = {}
    family_to_genes: dict[str, set[str]] = defaultdict(set)
    duplicate_gene_rows = 0

    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(
                    f"{path}:{line_number}: expected 2 tab-separated columns, "
                    f"found {len(fields)}"
                )
            gene_id, family_id = (field.strip() for field in fields)
            if not gene_id or not family_id:
                raise ValueError(f"{path}:{line_number}: blank gene or family ID")
            previous = gene_to_family.get(gene_id)
            if previous is not None:
                duplicate_gene_rows += 1
                if previous != family_id:
                    raise ValueError(
                        f"{path}:{line_number}: gene {gene_id!r} is assigned to both "
                        f"{previous!r} and {family_id!r}"
                    )
                continue
            gene_to_family[gene_id] = family_id
            family_to_genes[family_id].add(gene_id)

    return gene_to_family, dict(family_to_genes), duplicate_gene_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("full", type=Path, help="Full-mode gene_group_table.tsv")
    parser.add_argument("supplemented", type=Path, help="Supplementation-mode gene_group_table.tsv")
    parser.add_argument(
        "--overlap-report", type=Path, default=Path("family_partition_overlaps.tsv"),
        help="Output TSV of overlapping Full and Supplementation families",
    )
    args = parser.parse_args()

    full_gene_to_family, full_families, full_dupes = load_gene_family_map(args.full)
    supp_gene_to_family, supp_families, supp_dupes = load_gene_family_map(args.supplemented)

    full_genes = set(full_gene_to_family)
    supp_genes = set(supp_gene_to_family)
    shared_genes = full_genes & supp_genes
    full_only = full_genes - supp_genes
    supp_only = supp_genes - full_genes

    full_sets = {frozenset(members) for members in full_families.values()}
    supp_sets = {frozenset(members) for members in supp_families.values()}
    exact_sets = full_sets & supp_sets

    full_families_split = 0
    for members in full_families.values():
        overlapping_supp_families = {supp_gene_to_family[g] for g in members if g in supp_gene_to_family}
        if len(overlapping_supp_families) > 1:
            full_families_split += 1

    supp_families_merging_full = 0
    for members in supp_families.values():
        overlapping_full_families = {full_gene_to_family[g] for g in members if g in full_gene_to_family}
        if len(overlapping_full_families) > 1:
            supp_families_merging_full += 1

    # Every non-identical family contributes its members to this tally. These
    # totals are reported separately by run because the partitions may differ.
    full_genes_in_nonidentical_families = sum(
        len(members) for members in full_families.values()
        if frozenset(members) not in supp_sets
    )
    supp_genes_in_nonidentical_families = sum(
        len(members) for members in supp_families.values()
        if frozenset(members) not in full_sets
    )

    overlap_rows = []
    for full_id, full_members in full_families.items():
        overlaps = defaultdict(int)
        for gene in full_members:
            supp_id = supp_gene_to_family.get(gene)
            if supp_id is not None:
                overlaps[supp_id] += 1
        for supp_id, shared_count in overlaps.items():
            supp_size = len(supp_families[supp_id])
            overlap_rows.append({
                "full_family_id": full_id,
                "supplemented_family_id": supp_id,
                "shared_gene_count": shared_count,
                "full_family_size": len(full_members),
                "supplemented_family_size": supp_size,
                "fraction_of_full_family": f"{shared_count / len(full_members):.6f}",
                "fraction_of_supplemented_family": f"{shared_count / supp_size:.6f}",
            })

    overlap_rows.sort(key=lambda row: (-row["shared_gene_count"], row["full_family_id"], row["supplemented_family_id"]))
    args.overlap_report.parent.mkdir(parents=True, exist_ok=True)
    with args.overlap_report.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "full_family_id", "supplemented_family_id", "shared_gene_count",
            "full_family_size", "supplemented_family_size",
            "fraction_of_full_family", "fraction_of_supplemented_family",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(overlap_rows)

    print("PyamilySeq family-partition comparison")
    print(f"Full file:                 {args.full}")
    print(f"Supplemented file:         {args.supplemented}")
    print(f"Genes in Full:             {len(full_genes):,}")
    print(f"Genes in Supplemented:     {len(supp_genes):,}")
    print(f"Shared gene IDs:           {len(shared_genes):,}")
    print(f"Only in Full:              {len(full_only):,}")
    print(f"Only in Supplemented:      {len(supp_only):,}")
    print(f"Families in Full:          {len(full_families):,}")
    print(f"Families in Supplemented:  {len(supp_families):,}")
    print(f"Exactly identical families:{len(exact_sets):,}")
    print(f"Full families split across Supplemented families: {full_families_split:,}")
    print(f"Supplemented families containing multiple Full families: {supp_families_merging_full:,}")
    print(f"Genes in non-identical Full families: {full_genes_in_nonidentical_families:,}")
    print(f"Genes in non-identical Supplemented families: {supp_genes_in_nonidentical_families:,}")
    print(f"Repeated identical rows in Full table: {full_dupes:,}")
    print(f"Repeated identical rows in Supplemented table: {supp_dupes:,}")
    print(f"Overlap table:             {args.overlap_report}")

    if full_only or supp_only:
        print("\nWARNING: The gene ID sets differ. Resolve input/extraction differences before interpreting partition differences.")
        if full_only:
            print("Example Full-only gene IDs:", ", ".join(sorted(full_only)[:10]))
        if supp_only:
            print("Example Supplemented-only gene IDs:", ", ".join(sorted(supp_only)[:10]))
    elif full_sets == supp_sets:
        print("\nRESULT: The two runs produced exactly the same family partition.")
    else:
        print("\nRESULT: The two runs contain the same genes but produced different family partitions.")
        print("Inspect the overlap table for splits and merges; family IDs themselves are not being compared.")


if __name__ == "__main__":
    main()
