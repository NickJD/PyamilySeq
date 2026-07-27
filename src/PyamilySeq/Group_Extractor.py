import os
import csv
import logging
import re

try:
    from .families import read_fasta_records
except ImportError:
    from families import read_fasta_records

# Use centralissed logger factory from constants
try:
    from .constants import configure_logger, LoggingArgumentParser
except Exception:
    from constants import configure_logger, LoggingArgumentParser


def parse_fasta(fasta_file):
    """
    Parses a FASTA file and returns a dictionary of gene IDs and sequences.
    """
    return read_fasta_records(fasta_file)


def parse_csv(csv_file):
    """
    Parses a CSV file to extract group IDs and gene IDs (skipping the first line).
    """
    groups = {}
    with open(csv_file, 'r') as f:
        reader = csv.reader(f, delimiter=',')  # Assuming tab-delimited CSV
        header = next(reader)
        genome_names = header[14:]
        for row in reader:
            group_id = row[0]
            gene_ids = []
            for genome, cell in zip(genome_names, row[14:]):
                for gene in re.split(r'[,;\s]+', cell.strip()):
                    if not gene:
                        continue
                    gene_ids.append(gene if '|' in gene else f'{genome}|{gene}')
            groups[group_id] = gene_ids
    return groups


def write_group_fastas(groups, sequences, output_dir):

    logger = logging.getLogger("PyamilySeq.Group_Extractor")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    by_gene = {}
    for sequence_id in sequences:
        gene = sequence_id.split('|', 1)[-1]
        by_gene.setdefault(gene, []).append(sequence_id)
    for group_id, gene_ids in groups.items():
        group_file = os.path.join(output_dir, f"{group_id}.fasta")
        with open(group_file, 'w') as f:
            for gene_id in gene_ids:
                sequence_id = gene_id if gene_id in sequences else None
                if sequence_id is None:
                    matches = by_gene.get(gene_id.split('|', 1)[-1], [])
                    sequence_id = matches[0] if len(matches) == 1 else None
                if sequence_id:
                    f.write(f">{sequence_id}\n{sequences[sequence_id]}\n")
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
