from collections import defaultdict
import hashlib
import logging
import os

# Use centralised logger factory from constants
try:
    from .constants import configure_logger, LoggingArgumentParser
except Exception:
    from constants import configure_logger, LoggingArgumentParser

def read_cd_hit_output(clstr_file):
    """
    Reads a CD-HIT .clstr file and extracts sequence clusters.
    Returns a dictionary where keys are sequence headers and values are cluster IDs.
    """
    seq_to_cluster = {}  # Maps sequence header -> cluster ID
    cluster_id = None

    with open(clstr_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                cluster_id = line.split(None, 1)[1]
            elif line:
                parts = line.split('\t')
                if len(parts) > 1:
                    seq_header = parts[1].split('>')[1].split('...')[0]  # Extract sequence header
                    seq_to_cluster[seq_header] = cluster_id

    return seq_to_cluster

def compare_cd_hit_clusters(file1, file2, output_file):
    logger = logging.getLogger("PyamilySeq.Group_Compare")
    logger.info("Comparing clusters: %s vs %s", file1, file2)
    # Read both clustering files
    clusters1 = read_cd_hit_output(file1)
    clusters2 = read_cd_hit_output(file2)

    # Reverse mappings: cluster ID -> list of sequences
    grouped_clusters1 = defaultdict(set)
    grouped_clusters2 = defaultdict(set)

    for seq, cluster_id in clusters1.items():
        grouped_clusters1[cluster_id].add(seq)
    for seq, cluster_id in clusters2.items():
        grouped_clusters2[cluster_id].add(seq)

    def signatures(grouped):
        result = {}
        for cluster_id, members in grouped.items():
            digest = hashlib.sha256('\0'.join(sorted(members)).encode('utf-8')).hexdigest()[:12]
            result[cluster_id] = digest
        return result

    signatures1 = signatures(grouped_clusters1)
    signatures2 = signatures(grouped_clusters2)
    sequence_shifts = 0
    only_in_file1 = defaultdict(list)
    only_in_file2 = defaultdict(list)

    # Prepare data for the TSV output
    tsv_data = []

    for seq in sorted(set(clusters1) | set(clusters2)):
        cluster1 = clusters1.get(seq)
        cluster2 = clusters2.get(seq)
        if cluster1 is None:
            only_in_file2[cluster2].append(seq)
            tsv_data.append([seq, 'NA', cluster2, 'Only in file2'])
        elif cluster2 is None:
            only_in_file1[cluster1].append(seq)
            tsv_data.append([seq, cluster1, 'NA', 'Only in file1'])
        elif grouped_clusters1[cluster1] != grouped_clusters2[cluster2]:
            sequence_shifts += 1
            tsv_data.append([seq, cluster1, cluster2, 'Partition mismatch'])

    equivalent_partitions = len(
        set(signatures1.values()) & set(signatures2.values())
    )

    # Print metrics
    logger.info("Clustering Comparison Metrics:")
    logger.info("Equivalent family partitions: %s", equivalent_partitions)
    logger.info("Sequence shifts (sequences assigned to different clusters): %s", sequence_shifts)
    logger.info("Sequences only in the first file: %s", sum(map(len, only_in_file1.values())))
    logger.info("Sequences only in the second file: %s", sum(map(len, only_in_file2.values())))

    # Write the results to a TSV file
    with open(output_file, 'w') as out_file:
        out_file.write("Sequence\tCluster ID (File 1)\tCluster ID (File 2)\tChange Type\n")
        for row in tsv_data:
            out_file.write("\t".join(map(str, row)) + "\n")

    logger.info("Results have been written to %s", output_file)

def main():
    # Early console-only logger so parser.description and argparse messages are emitted via logger
    early_logger = configure_logger("PyamilySeq.Group_Compare", enable_file=False, log_dir=None, verbose=False)
    parser = LoggingArgumentParser(logger_name="PyamilySeq.Group_Compare", description="Running Group-Compare - A tool to compare two CD-HIT .clstr files to check for clustering consistency.")

    parser.add_argument("-file1", required=True, help="First CD-HIT .clstr file")
    parser.add_argument("-file2", required=True, help="Second CD-HIT .clstr file")
    parser.add_argument("-output", required=True, help="Output file (TSV format)")
    parser.add_argument("--log", action="store_true", dest="log", help="Create a timestamped logfile for this run.")
    parser.add_argument("--log-dir", dest="log_dir", default=None, help="Directory for logfile (default: same dir as -output).")
    args = parser.parse_args()

    # Setup logger
    out_dir = os.path.abspath(os.path.dirname(args.output)) if args.output else os.getcwd()
    log_dir = args.log_dir if args.log_dir else out_dir
    logger = configure_logger("PyamilySeq.Group_Compare", enable_file=args.log, log_dir=log_dir, verbose=False)

    compare_cd_hit_clusters(args.file1, args.file2, args.output)

if __name__ == "__main__":
    main()
