import argparse
import os
import csv


def parse_fasta_stats(fasta_file):
    """
    Parses a FASTA file and calculates sequence statistics.
    """
    lengths = []
    with open(fasta_file, 'r') as f:
        sequence = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sequence:  # Save the previous sequence length
                    lengths.append(len(''.join(sequence)))
                sequence = []  # Reset for the next sequence
            else:
                sequence.append(line)
        if sequence:  # Save the last sequence length
            lengths.append(len(''.join(sequence)))

    # Calculate statistics
    num_sequences = len(lengths)
    if num_sequences > 0:
        avg_length = sum(lengths) / num_sequences
        min_length = min(lengths)
        max_length = max(lengths)
        length_diff = max_length - min_length
        percent_diff = (length_diff / min_length * 100) if min_length > 0 else 0
    else:
        avg_length = min_length = max_length = length_diff = percent_diff = 0

    return {
        "num_sequences": num_sequences,
        "min_length": min_length,
        "max_length": max_length,
        "avg_length": avg_length,
        "length_diff": length_diff,
        "percent_diff": percent_diff
    }


def process_fasta_directory(input_dir, output_csv):
    """
    Processes a directory of FASTA files and writes statistics to a CSV file.
    """
    results = []
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            file_path = os.path.join(input_dir, filename)
            stats = parse_fasta_stats(file_path)
            results.append({
                "file_name": filename,
                "num_sequences": stats["num_sequences"],
                "min_length": stats["min_length"],
                "max_length": stats["max_length"],
                "avg_length": stats["avg_length"],
                "length_diff": stats["length_diff"],
                "percent_diff": stats["percent_diff"]
            })

    # Write results to a CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ["file_name", "num_sequences", "min_length", "max_length", "avg_length", "length_diff",
                      "percent_diff"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


def main():
    parser = argparse.ArgumentParser(description="Summarize sequence statistics for a directory of FASTA files.")
    parser.add_argument("-input_dir", required=True, help="Directory containing FASTA files.")
    parser.add_argument("-output_csv", required=True, help="Output CSV file to save statistics.")

    args = parser.parse_args()

    # Process the directory of FASTA files
    print("Processing FASTA files...")
    process_fasta_directory(args.input_dir, args.output_csv)
    print(f"Statistics saved to {args.output_csv}")


if __name__ == "__main__":
    main()
