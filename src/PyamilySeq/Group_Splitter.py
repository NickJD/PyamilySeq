import subprocess
import os
import argparse
from collections import defaultdict, OrderedDict
from line_profiler_pycharm import profile

try:
    from .Constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from Constants import *
    from utils import *

def run_cd_hit(options, input_file, clustering_output, clustering_mode):
    cdhit_command = [
        clustering_mode,
        '-i', input_file,
        '-o', clustering_output,
        '-c', str(options.pident),
        '-s', str(options.len_diff),
        '-T', str(options.clustering_threads),
        '-M', str(options.clustering_memory),
        '-d', "0",
        '-sc', "1",
        '-sf', "1"
    ]
    if options.verbose:
        subprocess.run(cdhit_command)
    else:
        subprocess.run(cdhit_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def calculate_new_rep_seq(cluster_data):
    total_length = sum(entry['length'] for entry in cluster_data)
    avg_length = total_length / len(cluster_data)

    total_identity = sum(entry['percent_identity'] for entry in cluster_data)
    avg_identity = total_identity / len(cluster_data)

    # Calculate a score based on both length difference and percent identity
    def score(entry):
        length_diff = abs(entry['length'] - avg_length)
        identity_diff = abs(entry['percent_identity'] - avg_identity)
        return length_diff + (100 - identity_diff)  # You can weight these differently

    rep_entry = min(cluster_data, key=score)
    return rep_entry


def length_within_threshold(rep_length, length, len_diff):
    return abs(rep_length - length) / rep_length <= len_diff


def check_if_all_identical(clustered_sequences):
    lengths = {entry['length'] for cluster in clustered_sequences.values() for entry in cluster}
    perc_idents = {entry['percent_identity'] for cluster in clustered_sequences.values() for entry in cluster}

    return len(lengths) == 1 and len(perc_idents) == 1


def read_fasta_groups(fasta_file):
    groups = defaultdict(list)
    genome_count = defaultdict(int)
    current_group = None
    current_sequence = []

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_group is not None:
                    groups[current_group].append((current_group_header, ''.join(current_sequence)))

                current_group_header = line.strip()
                current_group = current_group_header.split('|')[0]
                genome = current_group_header.split('|')[1]
                current_sequence = []
                genome_count[genome] += 1
            else:
                current_sequence.append(line.strip())

        if current_group is not None:
            groups[current_group].append((current_group_header, ''.join(current_sequence)))

    return groups, genome_count


def write_fasta(sequences, output_file):
    with open(output_file, 'w') as f:
        for header, seq in sequences:
            f.write(f"{header}\n{seq}\n")


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
                    length = clustered_info.split(',')[0]
                    length = int(''.join(c for c in length if c.isdigit()))
                    clustered_header = clustered_info.split('>')[1].split('...')[0]
                    clustered_header = '>' + clustered_header

                    if 'at' in clustered_info:
                        percent_identity = extract_identity(line)

                    elif '*' in line:
                        percent_identity = 100.0
                    else:
                        raise ValueError("Percent identity not found in the string.")

                    clusters[current_cluster_id].append({
                        'header': clustered_header,
                        'length': length,
                        'percent_identity': percent_identity
                    })

    return clusters


def separate_groups(input_fasta, options, clustering_mode):
    groups, genome_count = read_fasta_groups(input_fasta)

    paralog_groups = defaultdict(int)  # To track number of paralog groups

    for group_header, sequences in groups.items():
        group_name = group_header.split('|')[0]  # Get the group part (e.g., '>Group_n')

        # Count genomes with more than one gene
        genome_to_gene_count = defaultdict(int)
        for header, _ in sequences:
            genome = header.split('|')[1]
            genome_to_gene_count[genome] += 1

        num_genomes_with_multiple_genes = sum(1 for count in genome_to_gene_count.values() if count > 1)
        total_genomes = len(genome_to_gene_count)

        # Check if the group meets the threshold for having paralogs
        if total_genomes == 0 or (num_genomes_with_multiple_genes / total_genomes) * 100 < options.percent_threshold:
            continue

        group_file_name = group_name.replace('>','')

        temp_fasta = f"{options.output_dir}/{group_file_name}.fasta"
        write_fasta(sequences, temp_fasta)

        # Run cd-hit on the individual group
        clustering_output = f"{options.output_dir}/{group_file_name}_clustering"

        run_cd_hit(options, temp_fasta, clustering_output, clustering_mode)

        # Read the clustering results to find subgroups
        clustered_sequences = read_cd_hit_output(clustering_output + '.clstr')

        # Detect if all sequences are identical in length and percentage identity
        all_same = check_if_all_identical(clustered_sequences)

        # **Global subgroup counter for the entire major group**
        subgroup_id = 0
        remaining_sequences = sequences.copy() # Track unprocessed sequences
        sequences_to_remove = []

        if not all_same:
            while remaining_sequences:
                # Track subgroups for this pass
                subgroup_sequences = []
                genome_seen = set()
                sequences_found = False  # Track if any sequence was added

                # Recalculate representative sequence dynamically based on remaining genes
                rep = calculate_new_rep_seq(
                    [entry for cluster in clustered_sequences.values() for entry in cluster if
                     entry['header'] in (h for h, _ in remaining_sequences)]
                )

                # Find the sequence corresponding to rep['header'] from the list of sequences
                rep_seq = next((seq for header, seq in sequences if header == rep['header']), None)

                # Process each genome to select the best matching sequence
                for genome in genome_to_gene_count:
                    best_sequence = None
                    best_score = -1  # Initialize with a very low similarity score

                    # Iterate over each sequence in the remaining sequences for this genome
                    for header, seq in remaining_sequences:
                        genome_id = header.split('|')[1]

                        if genome_id == genome:  # Ensure this sequence belongs to the current genome

                            length = len(seq)
                            if rep_seq == seq:
                                perc_ident = 100.0
                            else:
                                perc_ident = calculate_similarity(rep_seq, seq)  # Define a function to calculate similarity

                            # Calculate the length difference ratio (smaller ratio means closer length to the representative)
                            length_diff_ratio = abs(rep['length'] - length) / rep['length']

                            # Check if this sequence is more similar than the current best one
                            if length_within_threshold(rep['length'], length,
                                                       options.len_diff) and perc_ident >= options.pident:

                                # Combine percentage identity and length difference into a single score
                                # Here, you want a high identity and a small length difference
                                # Adjust the weight of length difference and similarity according to your requirements
                                score = perc_ident - (length_diff_ratio * 100)  # Weighting length diff (you can adjust the *100 factor)

                                # Check if this sequence has a higher score than the current best
                                if score > best_score:
                                    best_score = score
                                    best_sequence = (header, seq)  # Store the best matching sequence for this genome

                    # Once the best sequence is identified, add it to the subgroup
                    if best_sequence is not None:
                        sequences_found = True  # At least one sequence was added
                        new_header = f">{group_file_name}_subgroup_{subgroup_id}|{best_sequence[0].split('|')[1]}|{best_sequence[0].split('|')[2]}"
                        subgroup_sequences.append((new_header, best_sequence[1]))
                        sequences_to_remove.append(best_sequence)
                        genome_seen.add(genome)

                    # If no sequences were found for this pass, exit the loop
                    # if not sequences_found:
                    #     break

                # Write each subgroup into a separate FASTA file
                if subgroup_sequences:
                    subgroup_file = f"{options.output_dir}/{group_file_name}_subgroup_{subgroup_id}.fasta"
                    write_fasta(subgroup_sequences, subgroup_file)

                    # Remove processed sequences from the remaining list
                    remaining_sequences = [item for item in remaining_sequences if
                                           item[0] not in {h for h, _ in sequences_to_remove}]

                    # Increment subgroup ID globally for the next subgroup
                    subgroup_id += 1
                    paralog_groups[group_name] += 1  # Count this group as a paralog group


        else:
            # Condition 2: If sequences are identical, distribute genes evenly into subgroups
            num_subgroups = 1000
            subgroup_sequences = defaultdict(list)  # Store sequences for each subgroup
            genome_count = defaultdict(int)  # Count how many genes have been assigned to each genome

            # Iterate over all sequences regardless of whether the genome has been seen
            for header, seq in sequences:
                genome = header.split('|')[1]

                # Determine the next subgroup for this genome
                subgroup_id = genome_count[genome] % num_subgroups
                new_header = f"{group_file_name}_subgroup_{subgroup_id}|{genome}|{header.split('|')[2]}"
                subgroup_sequences[subgroup_id].append((new_header, seq))

                # Increment the count for this genome
                genome_count[genome] += 1

            # Write out each subgroup to a separate FASTA file
            for subgroup_id, seqs in subgroup_sequences.items():
                subgroup_file = f"{options.output_dir}/{group_file_name}_subgroup_{subgroup_id}.fasta"
                write_fasta(seqs, subgroup_file)

                # Increment subgroup ID globally for the next subgroup
                subgroup_id += 1
                paralog_groups[group_name] += 1  # Count this group as a paralog group



        # Clean up temporary fasta file if the option is set
        if options.delete_temp_files:
            if temp_fasta and os.path.exists(temp_fasta):
                os.remove(temp_fasta)
            if os.path.exists(clustering_output + '.clstr'):
                os.remove(clustering_output + '.clstr')
            if os.path.exists(clustering_output):
                os.remove(clustering_output)

    # Print metrics about paralog groups
    print(f"Identified {len(paralog_groups)} paralog groups:")
    for group_id, count in paralog_groups.items():
        print(f"Group ID: {group_id}, Number of new groups: {count}")


def main():
    parser = argparse.ArgumentParser(description='Group-Splitter: ' + PyamilySeq_Version + ': A tool to split "paralogous" groups identified by PyamilySeq.')
    ### Required Arguments
    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-input_fasta', action='store', dest='input_fasta',
                          help='Input FASTA file containing gene groups.',
                          required=True)
    required.add_argument('-sequence_type', action='store', dest='sequence_type', default='DNA',choices=['AA', 'DNA'],
                          help='Default - DNA: Are groups "DNA" or "AA" sequences?',
                          required=False)
    required.add_argument('-output_dir', action='store', dest='output_dir',
                          help='Output directory.',
                          required=True)

    optional = parser.add_argument_group('Optional Arguments')

    optional.add_argument('-pident', action='store', dest='pident', type=float, default=0.9,
                          help='Sequence identity threshold (default: 0.9)')
    optional.add_argument('-len_diff', action='store', dest='len_diff', type=float, default=0.05,
                          help='Length difference threshold (default: 0.05)')
    optional.add_argument('-clustering_threads', action='store', dest='clustering_threads', type=int, default=4,
                          help='Number of threads for clustering (default: 4)')
    optional.add_argument('-clustering_memory', action='store', dest='clustering_memory', type=int, default=2000,
                          help='Memory limit in MB for clustering (default: 2000)')
    optional.add_argument('-percent_threshold', action='store', dest='percent_threshold', type=float, default=80,
                          help='Minimum percentage of genomes with paralogs (default: 80.0)')
    optional.add_argument('-verbose', action='store_true', dest='verbose', help='Print verbose output.')
    optional.add_argument('-no_delete_temp_files', action='store_false', dest='delete_temp_files',
                          help='Default: Delete all temporary files after processing.')

    misc = parser.add_argument_group('Misc Arguments')
    misc.add_argument('-v', action='store_true', dest='version',
                      help='Print out version number and exit',
                      required=False)

    options = parser.parse_args()

    # Check for version flag
    if options.version:
        print(f"Group-Splitter version {PyamilySeq_Version}")
        exit(0)

    options = parser.parse_args()

    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)

    if options.sequence_type == 'DNA':
        clustering_mode = 'cd-hit-est'
    else:
        clustering_mode = 'cd-hit'

    separate_groups(options.input_fasta, options, clustering_mode)

    print("Done")


if __name__ == "__main__":
    main()
