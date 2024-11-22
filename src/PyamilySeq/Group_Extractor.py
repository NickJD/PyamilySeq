import argparse
import os
import csv


def parse_fasta(fasta_file):
    """
    Parses a FASTA file and returns a dictionary of gene IDs and sequences.
    """
    sequences = {}
    with open(fasta_file, 'r') as f:
        gene_id = None
        sequence = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if gene_id:  # Save the previous gene
                    sequences[gene_id] = ''.join(sequence)
                gene_id = line[1:].split()[0].split('|')[1].replace('ENSB_','')  # Extract the gene ID after ">"
                sequence = []
            else:
                sequence.append(line)
        if gene_id:  # Save the last gene
            sequences[gene_id] = ''.join(sequence)
    return sequences


def parse_csv(csv_file):
    """
    Parses a CSV file to extract group IDs and gene IDs (skipping the first line).
    """
    groups = {}
    with open(csv_file, 'r') as f:
        reader = csv.reader(f, delimiter=',')  # Assuming tab-delimited CSV
        next(reader)  # Skip the first line
        for row in reader:
            group_id = row[0]
            gene_ids = row[14:]  # Read from column 14 onward
            gene_ids = [gene.strip() for genes in gene_ids for gene in genes.split(';') if
                        gene.strip()]  # Flatten and clean
            groups[group_id] = gene_ids
    return groups


def write_group_fastas(groups, sequences, output_dir):
    """
    Writes individual FASTA files for each group with the relevant sequences.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for group_id, gene_ids in groups.items():
        group_file = os.path.join(output_dir, f"{group_id}.fasta")
        with open(group_file, 'w') as f:
            for gene_id in gene_ids:
                if gene_id in sequences:
                    f.write(f">{gene_id}\n{sequences[gene_id]}\n")
                else:
                    print(f"Warning: Gene ID {gene_id} not found in FASTA file.")


def main():
    parser = argparse.ArgumentParser(description="Process FASTA and CSV files to create grouped FASTA outputs.")
    parser.add_argument("-fasta", required=True, help="Input FASTA file containing gene sequences.")
    parser.add_argument("-csv", required=True, help="Input CSV file containing group and gene information.")
    parser.add_argument("-output_dir", required=True, help="Directory to save the grouped FASTA files.")

    args = parser.parse_args()

    # Parse the input files
    print("Parsing FASTA file...")
    sequences = parse_fasta(args.fasta)
    print("Parsing CSV file...")
    groups = parse_csv(args.csv)

    # Write the grouped FASTA files
    print("Writing grouped FASTA files...")
    write_group_fastas(groups, sequences, args.output_dir)
    print("Process completed successfully.")


if __name__ == "__main__":
    main()
