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


def load_genome_groups(species_calls_file, species2lineage, rank):
    rank_idx = RANK_INDEX[rank]
    genome2group = {}

    df = pd.read_csv(species_calls_file, sep="\t")

    for _, row in df.iterrows():
        genome = row["Sample"]
        species = row["Species"]

        if species not in species2lineage:
            continue

        lineage = species2lineage[species]

        if abs(rank_idx) > len(lineage):
            group = lineage[-1]
        else:
            group = lineage[rank_idx]

        genome2group[genome] = group

    return genome2group


def parse_pyamily_gene_to_family(mapping_file):
    """
    Expect two columns:
    GeneID    FamilyID
    """
    gene2fam = {}

    df = pd.read_csv(mapping_file, sep="\t", header=None)

    for _, row in df.iterrows():
        gene2fam[row[0]] = row[1]

    return gene2fam


def extract_gene_order(gene_id):
    """
    SAMEA7551658_00012 -> 12
    """
    try:
        return int(gene_id.split("_")[-1])
    except:
        return 10**12


def load_genes_from_fasta_headers(fasta_file):
    """
    Extract gene order per genome from FASTA headers
    """
    genome2genes = defaultdict(list)

    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                gene = line.strip()[1:]
                genome = gene.split("_")[0]
                genome2genes[genome].append(gene)

    for genome in genome2genes:
        genome2genes[genome].sort(key=extract_gene_order)

    return genome2genes


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--fasta", required=True,
                    help="protein FASTA with headers GenomeID_geneNumber")

    ap.add_argument("--gene-map", required=True,
                    help="GeneID -> PyamilySeq FamilyID table")

    ap.add_argument("--species-calls", required=True)
    ap.add_argument("--lineages", required=True)

    ap.add_argument("--rank", default="genus",
                    choices=RANK_INDEX.keys())

    ap.add_argument("--output", required=True)

    args = ap.parse_args()

    species2lineage = load_lineages(args.lineages)
    genome2group = load_genome_groups(args.species_calls, species2lineage, args.rank)

    gene2fam = parse_pyamily_gene_to_family(args.gene_map)
    genome2genes = load_genes_from_fasta_headers(args.fasta)

    with open(args.output, "w") as out:

        for genome, genes in genome2genes.items():

            group = genome2group.get(genome, "NA")

            families = []

            for gene in genes:
                fam = gene2fam.get(gene, "NA")
                families.append(fam)

            out.write(
                genome + "\t" +
                group + "\t" +
                "\t".join(families) + "\n"
            )


if __name__ == "__main__":
    main()