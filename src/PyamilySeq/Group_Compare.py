import argparse
from collections import defaultdict

def read_cd_hit_output(clstr_file):
    """
    Reads a CD-HIT .clstr file and extracts sequence clusters.
    Returns a dictionary where keys are sequence headers and values are cluster IDs.
    """
    seq_to_cluster = {}  # Maps sequence header -> cluster ID
    cluster_id = 0  # Generic ID for clusters (since CD-HIT names don't matter)

    with open(clstr_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                cluster_id += 1  # Increment cluster ID
            elif line:
                parts = line.split('\t')
                if len(parts) > 1:
                    seq_header = parts[1].split('>')[1].split('...')[0]  # Extract sequence header
                    seq_to_cluster[seq_header] = cluster_id

    return seq_to_cluster

def compare_cd_hit_clusters(file1, file2, output_file):
    """
    Compares two CD-HIT .clstr files to check if clusters are the same.
    Writes the results to a TSV file.
    """
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

    # Initialize metrics counters
    cluster_name_changes = 0
    sequence_shifts = 0
    only_in_file1 = defaultdict(list)
    only_in_file2 = defaultdict(list)
    cluster_mismatches = defaultdict(list)

    # Prepare data for the TSV output
    tsv_data = []

    # Track changes
    for seq, cluster_id in clusters1.items():
        if seq not in clusters2:
            only_in_file1[cluster_id].append(seq)
            tsv_data.append([seq, cluster_id, "NA", "Only in file1"])
        elif clusters2[seq] != cluster_id:
            # Sequence shifts: sequence in different clusters between files
            sequence_shifts += 1
            cluster_mismatches[seq].append((cluster_id, clusters2[seq]))
            tsv_data.append([seq, cluster_id, clusters2[seq], "Mismatch"])

    for seq, cluster_id in clusters2.items():
        if seq not in clusters1:
            only_in_file2[cluster_id].append(seq)
            tsv_data.append([seq, "NA", cluster_id, "Only in file2"])
        elif clusters1[seq] != cluster_id:
            # Sequence shifts: sequence in different clusters between files
            sequence_shifts += 1
            cluster_mismatches[seq].append((clusters1[seq], cluster_id))
            tsv_data.append([seq, clusters1[seq], cluster_id, "Mismatch"])

    # Track cluster name changes (same sequences in different clusters)
    for cluster_id1, seqs1 in grouped_clusters1.items():
        for cluster_id2, seqs2 in grouped_clusters2.items():
            if seqs1 == seqs2 and cluster_id1 != cluster_id2:
                cluster_name_changes += 1
                for seq in seqs1:
                    tsv_data.append([seq, cluster_id1, cluster_id2, "Cluster name change"])

    # Print metrics
    print("🔢 Clustering Comparison Metrics:")
    print(f"Cluster name changes: {cluster_name_changes}")
    print(f"Sequence shifts (sequences assigned to different clusters): {sequence_shifts}")
    print(f"Sequences only in the first file: {len(only_in_file1)}")
    print(f"Sequences only in the second file: {len(only_in_file2)}")
    print()

    # Write the results to a TSV file
    with open(output_file, 'w') as out_file:
        out_file.write("Sequence\tCluster ID (File 1)\tCluster ID (File 2)\tChange Type\n")
        for row in tsv_data:
            out_file.write("\t".join(map(str, row)) + "\n")

    print(f"✅ Results have been written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Compare two CD-HIT .clstr files to check for clustering consistency.")
    parser.add_argument("-file1", required=True, help="First CD-HIT .clstr file")
    parser.add_argument("-file2", required=True, help="Second CD-HIT .clstr file")
    parser.add_argument("-output", required=True, help="Output file (TSV format)")
    args = parser.parse_args()

    compare_cd_hit_clusters(args.file1, args.file2, args.output)

if __name__ == "__main__":
    main()
