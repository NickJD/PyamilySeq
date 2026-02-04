
import os
import csv
import logging


# Use centralised logger factory from constants
try:
    from .constants import configure_logger, LoggingArgumentParser
except Exception:
    from constants import configure_logger, LoggingArgumentParser


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
    logger = logging.getLogger("PyamilySeq.Group_Sizes")
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
    logger.info("Wrote statistics for %d FASTA files to %s", len(results), output_csv)


def main():
    # Early console-only logger so the parser.description is emitted via logger before argparse prints usage/help.
    early_logger = configure_logger("PyamilySeq.Group_Sizes", enable_file=False, log_dir=None, verbose=False)
    parser = LoggingArgumentParser(logger_name="PyamilySeq.Group_Sizes", description="Group-Sizes - A tool to summarise sequence statistics for a directory of FASTA files.")
    parser.add_argument("-input_dir", required=True, help="Directory containing FASTA files.")
    parser.add_argument("-output_csv", required=True, help="Output CSV file to save statistics.")
    parser.add_argument("--log", action="store_true", dest="log", help="Create a timestamped logfile for this run.")
    parser.add_argument("--log-dir", dest="log_dir", default=None, help="Directory for logfile (default: same dir as -output_csv).")

    args = parser.parse_args()

    out_dir = os.path.abspath(os.path.dirname(args.output_csv)) if args.output_csv else os.getcwd()
    log_dir = args.log_dir if args.log_dir else out_dir
    logger = configure_logger("PyamilySeq.Group_Sizes", enable_file=args.log, log_dir=log_dir, verbose=False)

    logger.info("Processing FASTA files in %s", args.input_dir)
    process_fasta_directory(args.input_dir, args.output_csv)
    logger.info("Statistics saved to %s", args.output_csv)


if __name__ == "__main__":
    main()
