import argparse
import collections
import csv


def parse_fasta_ids(fasta_file):
    """Extract IDs from the FASTA file."""
    ids = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].strip().split()[0]  # Capture the ID after '>'
                ids.append(seq_id)
    return ids


def find_ids_in_csv(ids, csv_file):
    """Search for each ID in the CSV file and report the first column where it is found."""
    found_records = collections.defaultdict(list)
    with open(csv_file, 'r') as f:
        csv_reader = csv.reader(f)
        for row in csv_reader:
            if row:  # Ensure row is not empty

                for id in ids: # slow
                    if id in row:
                        found_records[row[0]].append(id)
    return found_records


def main():
    parser = argparse.ArgumentParser(description="Extract IDs from a FASTA file and search for them in a CSV file.")
    parser.add_argument("-in", action='store', dest='fasta_file',
                        help="Input FASTA file", required=True)
    parser.add_argument("-ids", action='store', dest='csv_file',
                        help="CSV file containing IDs to search for", required=True)
    parser.add_argument("-out", action='store', dest='output_file',
                        help="Output file to save found IDs", required=True)

    options = parser.parse_args()

    # Parse IDs from the FASTA file
    ids = parse_fasta_ids(options.fasta_file)

    # Find IDs in the CSV file
    found_records = find_ids_in_csv(ids, options.csv_file)

    # Write output
    with open(options.output_file, 'w') as output:
        output.write("ID,Found_In_First_Column\n")
        for seq_id, found_in in found_records.items():
            output.write(f"{seq_id},{found_in}\n")


if __name__ == "__main__":
    main()
