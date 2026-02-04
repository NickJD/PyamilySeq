
import copy
import os

# Use centralised logger factory
try:
    from .constants import configure_logger, LoggingArgumentParser
except Exception:
    from constants import configure_logger, LoggingArgumentParser

def find_gene_ids_in_csv(csv_file, group_name):
    """Find gene IDs associated with the specified group name in the CSV file, starting from column 14."""
    gene_ids = []
    with open(csv_file, 'r') as f:
        for line in f:
            cells = line.strip().split(',')
            if cells[0].replace('"','') == group_name:
                # Collect gene IDs from column 14 onward
                # for cell in cells[14:]:
                #     gene_ids.extend(cell.strip().replace('"','').split())  # Splitting by spaces if there are multiple IDs in a cell                break
                for cell in cells[14:]:
                    for gene in cell.strip().replace('"', '').split(';'):
                        if gene:
                            gene_ids.append(gene)

    return gene_ids

def extract_sequences(fasta_file, gene_ids):
    """Extract sequences from the FASTA file that match the gene IDs."""
    sequences = {}
    capture = False
    current_id = ""
    not_found = copy.deepcopy(gene_ids)
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract the ID part after '>' and check if it's in gene_ids
                current_id = line[1:].strip().split()[0].split('|')[1]
                capture = current_id in gene_ids
                if current_id in not_found:
                    not_found.remove(current_id)
                if capture:
                    sequences[current_id] = [line.strip()]  # Start with header line
            elif capture:
                sequences[current_id].append(line.strip())  # Append sequence lines
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
