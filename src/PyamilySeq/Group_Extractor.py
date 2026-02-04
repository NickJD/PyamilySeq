import os
import csv
import logging

# Use centralissed logger factory from constants
try:
    from .constants import configure_logger, LoggingArgumentParser
except Exception:
    from constants import configure_logger, LoggingArgumentParser


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

    logger = logging.getLogger("PyamilySeq.Group_Extractor")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for group_id, gene_ids in groups.items():
        group_file = os.path.join(output_dir, f"{group_id}.fasta")
        with open(group_file, 'w') as f:
            for gene_id in gene_ids:
                if gene_id in sequences:
                    f.write(f">{gene_id}\n{sequences[gene_id]}\n")
                else:
                    logger.warning("Warning: Gene ID %s not found in FASTA file.", gene_id)


def main():
    # Early console-only logger so the parser description is logged before argparse outputs.
    early_logger = configure_logger("PyamilySeq.Group_Extractor", enable_file=False, log_dir=None, verbose=False)
    parser = LoggingArgumentParser(logger_name="PyamilySeq.Group_Extractor", description="Running Group-Extractor - A tool to process FASTA and CSV files to create grouped FASTA outputs.")

    parser.add_argument("-fasta", required=True, help="Input FASTA file containing gene sequences.")
    parser.add_argument("-csv", required=True, help="Input CSV file containing group and gene information.")
    parser.add_argument("-output_dir", required=True, help="Directory to save the grouped FASTA files.")
    parser.add_argument("--log", action="store_true", dest="log", help="Create a timestamped logfile for this run.")
    parser.add_argument("--log-dir", dest="log_dir", default=None, help="Directory for logfile (default: output_dir).")

    args = parser.parse_args()

    # Setup logger writing to output_dir (optional file)
    log_dir = os.path.abspath(args.output_dir) if args.output_dir else os.getcwd()
    if hasattr(args, "log_dir") and args.log_dir:
        log_dir = args.log_dir
    # Only create a logfile when --log is provided; default is console (stdout) only.
    logger = configure_logger("PyamilySeq.Group_Extractor", enable_file=getattr(args, "log", False), log_dir=log_dir, verbose=False)

    logger.info("Parsing FASTA file: %s", args.fasta)
    sequences = parse_fasta(args.fasta)
    logger.info("Parsed %d sequences.", len(sequences))
    logger.info("Parsing CSV file: %s", args.csv)
    groups = parse_csv(args.csv)
    logger.info("Parsed %d groups.", len(groups))

    logger.info("Writing grouped FASTA files to %s", args.output_dir)
    write_group_fastas(groups, sequences, args.output_dir)
    logger.info("Process completed successfully.")


if __name__ == "__main__":
    main()
