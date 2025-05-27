import argparse
import csv

def load_csv_ids(csv_file):
    """Load all IDs from a single-row CSV file into a set for fast lookup."""
    id_set = set()
    with open(csv_file, 'r') as f:
        csv_reader = csv.reader(f)
        for row in csv_reader:
            for item in row:  # Process all values in the single row
                id_set.add(item.strip())  # Strip whitespace and add to set
    return id_set



def extract_fasta_sequences(fasta_file, id_set, output_file):
    """Extract sequences from the FASTA file if their ID is in the CSV."""
    with open(fasta_file, 'r') as infile, open(output_file, 'w') as outfile:
        write_sequence = False
        for line in infile:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]  # Extract the sequence ID
                write_sequence = seq_id in id_set  # Check if ID is in the set
            if write_sequence:
                outfile.write(line)  # Write header or sequence


def main():
    parser = argparse.ArgumentParser(description="Extract FASTA sequences based on a list of IDs from a CSV file.")
    parser.add_argument("-in", dest='fasta_file', required=True, help="Input FASTA file")
    parser.add_argument("-ids", dest='csv_file', required=True, help="CSV file containing IDs")
    parser.add_argument("-out", dest='output_file', required=True, help="Output FASTA file")

    options = parser.parse_args()

    # Load IDs from the CSV file
    id_set = load_csv_ids(options.csv_file)

    # Extract matching sequences from the FASTA file
    extract_fasta_sequences(options.fasta_file, id_set, options.output_file)

    print(f"Extraction complete. Sequences saved to {options.output_file}")

if __name__ == "__main__":
    main()
