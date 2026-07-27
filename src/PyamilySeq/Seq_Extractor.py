
import copy
import csv
import os
import re

try:
    from .families import read_fasta_records
except ImportError:
    from families import read_fasta_records

# Use centralised logger factory
try:
    from .constants import configure_logger, LoggingArgumentParser
except Exception:
    from constants import configure_logger, LoggingArgumentParser

def find_gene_ids_in_csv(csv_file, group_name):
    """Find gene IDs associated with the specified group name in the CSV file, starting from column 14."""
    gene_ids = []
    with open(csv_file, newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        genomes = header[14:]
        for row in reader:
            if row and row[0] == group_name:
                for genome, cell in zip(genomes, row[14:]):
                    for gene in re.split(r'[,;\s]+', cell.strip()):
                        if gene:
                            gene_ids.append(gene if '|' in gene else f'{genome}|{gene}')
                break

    return gene_ids

def extract_sequences(fasta_file, gene_ids):
    """Extract sequences from the FASTA file that match the gene IDs."""
    records = read_fasta_records(fasta_file)
    by_gene = {}
    for sequence_id in records:
        by_gene.setdefault(sequence_id.split('|', 1)[-1], []).append(sequence_id)
    sequences = {}
    for requested in gene_ids:
        sequence_id = requested if requested in records else None
        if sequence_id is None:
            matches = by_gene.get(requested.split('|', 1)[-1], [])
            sequence_id = matches[0] if len(matches) == 1 else None
        if sequence_id:
            sequences[sequence_id] = [f'>{sequence_id}', records[sequence_id]]
    return sequences

def main():
    # Early console-only logger so parser.description appears in logger output before argparse prints the menu.
    early_logger = configure_logger("PyamilySeq.Seq_Extractor", enable_file=False, log_dir=None, verbose=False)
    parser = LoggingArgumentParser(logger_name="PyamilySeq.Seq_Extractor", description="Running Seq-Extractor - A tool to extract sequences for specified group name from CSV file and corresponding FASTA file.")

    parser.add_argument("-csv", action='store', dest='csv_file',
                        help="CSV file containing group data", required=True)
    parser.add_argument("-group", action='store', dest='group_name',
                        help="Group name to search for in the CSV", required=True)
    parser.add_argument("-fasta", action='store', dest='fasta_file',
                        help="Input FASTA file containing sequences", required=True)
    parser.add_argument("-out", action='store', dest='output_file',
                        help="Output FASTA file with extracted sequences", required=True)
    parser.add_argument("--log", action="store_true", dest="log", help="Create a timestamped logfile for this run.")
    parser.add_argument("--log-dir", dest="log_dir", default=None, help="Directory for logfile (default: dir of output_file).")

    options = parser.parse_args()

    # Setup logger
    out_dir = os.path.abspath(os.path.dirname(options.output_file)) if options.output_file else os.getcwd()
    log_dir = options.log_dir if getattr(options, "log_dir", None) else out_dir
    logger = configure_logger("PyamilySeq.Seq_Extractor", enable_file=getattr(options, "log", False), log_dir=log_dir, verbose=False)

    logger.info("Searching for gene IDs in CSV %s for group %s", options.csv_file, options.group_name)

    # Find gene IDs in CSV
    gene_ids = find_gene_ids_in_csv(options.csv_file, options.group_name)
    if not gene_ids:
        logger.warning("No gene IDs found for group name '%s' in the CSV.", options.group_name)
        print(f"No gene IDs found for group name '{options.group_name}' in the CSV.")
        return

    # Extract sequences from the FASTA file
    logger.info("Extracting sequences from FASTA: %s", options.fasta_file)
    sequences = extract_sequences(options.fasta_file, gene_ids)

    # Write matched sequences to the output FASTA file
    with open(options.output_file, 'w') as output:
        for gene_id, sequence_lines in sequences.items():
            output.write("\n".join(sequence_lines) + "\n")
    logger.info("Wrote %d sequences to %s", len(sequences), options.output_file)

if __name__ == "__main__":
    main()
