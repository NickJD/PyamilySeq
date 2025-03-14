import argparse
import csv
import sys
from collections import defaultdict

# Increase field size limit for large CSVs
csv.field_size_limit(sys.maxsize)


def read_gene_groups(csv_file, column_limit=100):
    """
    Reads a CSV file and extracts gene groups, considering genes per genome in each column.

    Parameters:
        csv_file (str): Path to the CSV file.
        column_limit (int): Maximum number of gene columns to consider.

    Returns:
        dict: A dictionary mapping sets of genes (frozensets) to their original group ID.
    """
    gene_groups = {}

    with open(csv_file, 'r') as f:
        reader = csv.reader(f, delimiter=",")  # Use tab as delimiter
        headers = next(reader)  # Read headers
        num_cols = len(headers)

        for row in reader:
            if len(row) < 15:  # Ensure enough columns
                continue

            # Extract up to `column_limit` gene columns (starting from column 15)
            genes_per_genome = row[14: min(14 + column_limit, num_cols)]

            # Flatten the lists into a single set of genes
            gene_names = frozenset(gene.strip() for genome in genes_per_genome for gene in genome.split() if gene)

            if gene_names:
                gene_groups[gene_names] = row[0]  # Store with group ID

    return gene_groups


def compare_gene_groups(baseline_file, comparison_file, column_limit=100):
    """
    Compares two gene clustering CSV files to detect gene rearrangements.

    Parameters:
        baseline_file (str): Path to the first (reference) CSV file.
        comparison_file (str): Path to the second (comparison) CSV file.
        column_limit (int): Maximum number of gene columns to consider.

    Returns:
        list: List of detected rearrangements.
    """
    baseline_groups = read_gene_groups(baseline_file, column_limit)
    comparison_groups = read_gene_groups(comparison_file, column_limit)

    rearrangements = []

    for gene_set, baseline_group in baseline_groups.items():
        if gene_set not in comparison_groups:
            rearrangements.append(f"Gene set from {baseline_group} has been reassigned or split.")

    for gene_set, comparison_group in comparison_groups.items():
        if gene_set not in baseline_groups:
            rearrangements.append(f"New gene set detected in {comparison_group}.")

    return rearrangements


def main():
    parser = argparse.ArgumentParser(description="Compare gene arrangements between two CSV files.")
    parser.add_argument("-b", "--baseline", required=True, help="Path to the first CSV file (baseline).")
    parser.add_argument("-c", "--comparison", required=True, help="Path to the second CSV file (to compare).")
    parser.add_argument("-o", "--output", required=True, help="Path to save the rearrangement report.")
    parser.add_argument("--limit", type=int, default=100, help="Number of gene columns to consider (default: 100).")

    args = parser.parse_args()

    print(
        f"Comparing {args.baseline} (baseline) against {args.comparison} (comparison), using {args.limit} gene columns...")

    rearrangements = compare_gene_groups(args.baseline, args.comparison, args.limit)

    with open(args.output, "w") as out_file:
        for line in rearrangements:
            out_file.write(line + "\n")

    print(f"Comparison complete. Results saved to {args.output}")


if __name__ == "__main__":
    main()
