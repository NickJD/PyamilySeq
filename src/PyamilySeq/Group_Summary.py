from collections import OrderedDict, defaultdict

try:
    from .constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError):
    from constants import *
    from utils import *


def categorise_percentage(percent):
    categories = {
        (20, 40): "20-40%",
        (40, 60): "40-60%",
        (60, 80): "60-80%",
        (80, 95): "80-95%",
        (95, 99): "95-99%",
        (99, 100): "99-100%"
    }
    for (low, high), label in categories.items():
        if low <= percent < high:
            return label
    return None


def read_cd_hit_output(clustering_output):
    clusters = OrderedDict()
    with open(clustering_output, 'r') as f:
        current_cluster_id = None
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_cluster_id = line.split(' ')[1]
                clusters[current_cluster_id] = []
            elif line and current_cluster_id is not None:
                parts = line.split('\t')
                if len(parts) > 1:
                    clustered_info = parts[1]
                    length = int(''.join(c for c in clustered_info.split(',')[0] if c.isdigit()))
                    clustered_header = '>' + clustered_info.split('>')[1].split('...')[0]
                    if 'at ' in clustered_info and '%' in clustered_info.split('at ')[-1]:
                        percent_identity = extract_identity(clustered_info)
                    elif line.endswith('*'):
                        percent_identity = 100.0
                    else:
                        raise ValueError("Percent identity not found in the string.")
                    clusters[current_cluster_id].append({
                        'header': clustered_header,
                        'length': length,
                        'percent_identity': percent_identity
                    })

    return clusters


def summarise_clusters(options, clusters, output):
    logger = logging.getLogger("PyamilySeq.Group_Summary")
    multicopy_groups = defaultdict(int)  # Counter for clusters with multicopy genes
    with open(output, 'w') as out_f:
        out_f.write(
            "Cluster_ID\tNum_Sequences\tNum_Genomes\tAvg_Length\tLength_Range\tAvg_Identity\tIdentity_Range\tGenomes_With_Multiple_Genes\tMulticopy_Percentage\n")
        for cluster_id, seqs in clusters.items():
            num_seqs = len(seqs)
            lengths = [seq['length'] for seq in seqs]
            identities = [seq['percent_identity'] for seq in seqs]
            avg_length = sum(lengths) / num_seqs if num_seqs > 0 else 0
            length_range = f"{min(lengths)}-{max(lengths)}" if num_seqs > 0 else "N/A"
            avg_identity = sum(identities) / num_seqs if num_seqs > 0 else 0
            identity_range = f"{min(identities):.2f}-{max(identities):.2f}" if num_seqs > 0 else "N/A"

            # Count genomes in cluster
            genome_to_gene_count = defaultdict(int)
            for seq in seqs:
                genome = seq['header'].split('|')[0].replace('>', '')
                genome_to_gene_count[genome] += 1

            num_genomes = len(genome_to_gene_count)
            num_genomes_with_multiple_genes = sum(1 for count in genome_to_gene_count.values() if count > 1)
            multicopy_percentage = (num_genomes_with_multiple_genes / options.genome_num) * 100 if options.genome_num > 0 else 0

            category = categorise_percentage(multicopy_percentage)
            if category:
                multicopy_groups[category] += 1

            # Write detailed output for each cluster
            out_f.write(
                f"{cluster_id}\t{num_seqs}\t{num_genomes}\t{avg_length:.2f}\t{length_range}\t{avg_identity:.2f}\t{identity_range}\t"
                f"{num_genomes_with_multiple_genes}\t{multicopy_percentage:.2f}\n"
            )

        # Define order for multicopy statistics output
        category_order = ["20-40%", "40-60%", "60-80%", "80-95%", "95-99%", "99-100%"]
        for category in category_order:
            logger.info("Clusters with multicopy genes in %s range: %s", category, multicopy_groups[category])


def main():
    # Initial logger setup before parsing arguments (use same logger name as summarise_clusters)
    early_logger = configure_logger("PyamilySeq.Group_Summary", enable_file=False, log_dir=None, verbose=False)
    # Use the LoggingArgumentParser so usage/help/error messages are emitted via the same logger
    parser = LoggingArgumentParser(logger_name="PyamilySeq.Group_Summary", description="Running Group-Summary - A tool to summarise CD-HIT clustering files.")

    # Required Arguments
    required = parser.add_argument_group('Required Parameters')
    required.add_argument('-input_cluster', action="store", dest="input_cluster", required=True,
                          help='Input CD-HIT .cluster file')
    required.add_argument('-output', action="store", dest="output", required=True,
                          help="Output TSV file to store cluster summaries - Will add '.tsv' if not provided by user")
    required.add_argument('-genome_num', action='store', dest='genome_num', type=int, required=True,
                          help='Total number of genomes in dataset')

    # Optional Arguments
    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-output_dir', action="store", dest="output_dir",
                          help='Default: Same as input file', required=False)

    misc = parser.add_argument_group("Misc Parameters")
    misc.add_argument("-verbose", action="store_true", dest="verbose",
                      help="Print verbose output.", required=False)
    misc.add_argument("-v", "--version", action="version",
                      version=f"PyamilySeq: Group-Summary version {PyamilySeq_Version} - Exiting",
                      help="Print out version number and exit")
    # Add optional logging flags
    parser.add_argument("--log", action="store_true", dest="log", help="Create a timestamped logfile for this run.")
    parser.add_argument("--log-dir", dest="log_dir", default=None, help="Directory for logfile (default: output_dir or input file dir).")

    options = parser.parse_args()

    # Setup logger once we know output paths/options
    # after we resolve output_path / options.output_dir:
    resolved_log_dir = options.log_dir if getattr(options, "log_dir", None) else (os.path.abspath(options.output_dir) if getattr(options, "output_dir", None) else os.getcwd())
    logger = configure_logger("PyamilySeq.Group_Summary", enable_file=getattr(options, "log", False), log_dir=resolved_log_dir, verbose=options.verbose)
    if options.verbose:
        logger.debug("Options: %s", vars(options))

    # File handling
    options.input_cluster = fix_path(options.input_cluster)
    if options.output_dir is None:
        options.output_dir = os.path.dirname(os.path.abspath(options.input_cluster))
    output_path = os.path.abspath(options.output_dir)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    output_name = options.output
    if not output_name.endswith('.tsv'):
        output_name += '.tsv'
    output_file_path = os.path.join(output_path, output_name)

    # Process clusters and generate summary
    clusters = read_cd_hit_output(options.input_cluster)
    summarise_clusters(options, clusters, output_file_path)
    logger.info("Summary written to %s", output_file_path)


if __name__ == "__main__":
    main()
