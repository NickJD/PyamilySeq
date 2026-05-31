#!/usr/bin/env python3

import argparse
import pandas as pd
from collections import defaultdict


RANK_INDEX = {
    "species": -1,
    "genus": -2,
    "family": -3,
    "order": -4,
    "class": -5,
    "phylum": -6,
    "domain": 0
}


def load_lineages(lineage_file):
    """
    species -> lineage list
    """
    species2lineage = {}

    with open(lineage_file) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            species = parts[0]
            lineage = parts[-1].split(";")
            species2lineage[species] = lineage

    return species2lineage


def load_genome_map(species_calls_file, species2lineage, rank):
    """
    genome -> taxon at chosen rank
    """
    rank_idx = RANK_INDEX[rank]
    genome2tax = {}

    df = pd.read_csv(species_calls_file, sep="\t")

    for _, row in df.iterrows():
        genome = row["Sample"]
        species = row["Species"]

        if species not in species2lineage:
            continue

        lineage = species2lineage[species]

        if abs(rank_idx) > len(lineage):
            tax = lineage[-1]
        else:
            tax = lineage[rank_idx]

        genome2tax[genome] = tax

    return genome2tax


def parse_pyamily_matrix(matrix_file):
    """
    rows = gene families
    columns = genomes or gene IDs depending on PyamilySeq output
    """
    df = pd.read_csv(matrix_file, sep="\t", index_col=0)
    return df


def build_presence_absence(df):
    """
    force binary presence/absence
    """
    return df.astype(bool).astype(int)


def filter_singletons(df):
    """
    remove gene families present in only 1 genome
    """
    return df[df.sum(axis=1) > 1]


def filter_group_specific(df, genome2tax, mode="remove"):
    """
    remove genes present in only one taxonomic group
    """
    if mode != "remove":
        return df

    tax_groups = defaultdict(set)

    for genome in df.columns:
        tax = genome2tax.get(genome, None)
        if tax:
            tax_groups[tax].add(genome)

    keep_rows = []

    for gene, row in df.iterrows():
        present_groups = set()

        for genome, val in row.items():
            if val > 0:
                tax = genome2tax.get(genome, None)
                if tax:
                    present_groups.add(tax)

        if len(present_groups) > 1:
            keep_rows.append(gene)

    return df.loc[keep_rows]


def attach_metadata(df, genome2tax):
    """
    build long-format table:
    genome | group | gene1 | gene2 ...
    """
    rows = []

    for genome in df.columns:
        group = genome2tax.get(genome, "NA")

        row = {"Genome": genome, "Group": group}

        for gene in df.index:
            row[gene] = int(df.loc[gene, genome])

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--matrix", required=True)
    ap.add_argument("--species-calls", required=True)
    ap.add_argument("--lineages", required=True)

    ap.add_argument("--rank", default="genus",
                    choices=RANK_INDEX.keys())

    ap.add_argument("--remove-singletons", action="store_true")
    ap.add_argument("--remove-group-specific", action="store_true")

    ap.add_argument("--output", required=True)
    ap.add_argument("--format", choices=["wide", "long"], default="wide")

    args = ap.parse_args()

    species2lineage = load_lineages(args.lineages)
    genome2tax = load_genome_map(args.species_calls, species2lineage, args.rank)

    df = parse_pyamily_matrix(args.matrix)
    df = build_presence_absence(df)

    if args.remove_singletons:
        df = filter_singletons(df)

    if args.remove_group_specific:
        df = filter_group_specific(df, genome2tax, mode="remove")

    if args.format == "long":
        out = attach_metadata(df, genome2tax)
    else:
        df.insert(0, "Group", [genome2tax.get(g, "NA") for g in df.columns])
        out = df

    out.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()