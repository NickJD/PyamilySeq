import collections
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
        '-g', "1",
        '-sc', "1",
        '-sf', "1"
    ]
    if options.verbose:
        subprocess.run(cdhit_command)
    else:
        subprocess.run(cdhit_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

@profile
def calculate_new_rep_seq(cluster_data, length_weight=1.0, identity_weight=1.0):
    total_length = sum(entry['length'] for entry in cluster_data)
    avg_length = total_length / len(cluster_data)

    total_identity = sum(entry['percent_identity'] for entry in cluster_data)
    avg_identity = total_identity / len(cluster_data)

    # Normalize length and identity
    max_length = max(entry['length'] for entry in cluster_data)
    max_identity = 100  # Assuming percent_identity is out of 100

    # Calculate a score based on both length difference and percent identity
    def score(entry):
        normalized_length_diff = abs(entry['length'] - avg_length) / max_length
        normalized_identity_diff = abs(entry['percent_identity'] - avg_identity) / max_identity
        return (length_weight * normalized_length_diff) + (identity_weight * (1 - normalized_identity_diff))

    rep_entry = min(cluster_data, key=score)
    return rep_entry



def length_within_threshold(rep_length, length, len_diff):
    return abs(rep_length - length) / rep_length <= len_diff


def check_if_all_identical(clustered_sequences):
    lengths = {entry['length'] for cluster in clustered_sequences.values() for entry in cluster}
    perc_idents = {entry['percent_identity'] for cluster in clustered_sequences.values() for entry in cluster}

    return len(lengths) == 1 and len(perc_idents) == 1



def read_fasta_groups(options):
    groups = defaultdict(list)
    genome_count = defaultdict(int)
    current_group = None
    current_sequence = []

    # Parse the list of specific group numbers if provided
    selected_groups = None
    if options.groups is not None:
        selected_groups = [int(g.strip()) for g in options.groups.split(',')]

    with open(options.input_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_group is not None and (selected_groups is None or group_number in selected_groups):
                    groups[current_group].append((current_group_header, ''.join(current_sequence)))

                current_group_header = line.strip()
                current_group = current_group_header.split('|')[0]
                genome = current_group_header.split('|')[1]
                current_sequence = []
                genome_count[genome] += 1

                # Only process if group matches the selected_groups or if no specific groups were provided
                group_number = int(current_group.replace('>Group_', ''))  # Assuming format 'Group_n'
                if selected_groups is not None and group_number not in selected_groups:
                    current_group = None  # Skip this group
                    continue

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

                    if 'at ' in clustered_info and '%' in clustered_info.split('at ')[-1]:
                        percent_identity = extract_identity(line)
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

@profile
def separate_groups(options, clustering_mode):
    groups, genome_count = read_fasta_groups(options)

    paralog_groups = defaultdict(int)  # To track number of paralog groups

    for group_header, sequences in groups.items():
        if options.verbose:
            print(f"\n###\nCurrent Group: {group_header.replace('>','')}\n")

        group_name = group_header.split('|')[0]  # Get the group part (e.g., '>Group_n')

        # Count genomes with more than one gene
        genome_to_gene_count = defaultdict(int)
        for header, _ in sequences:
            genome = header.split('|')[1]
            genome_to_gene_count[genome] += 1

        num_genomes_with_multiple_genes = sum(1 for count in genome_to_gene_count.values() if count > 1)

        # Check if the group meets the threshold for having paralogs
        if options.groups == None:
            if (num_genomes_with_multiple_genes / options.genome_num) * 100 < options.group_threshold:
                continue


        group_file_name = group_name.replace('>','')

        temp_fasta = f"{options.output_dir}/{group_file_name}.fasta"
        write_fasta(sequences, temp_fasta)

        # Run cd-hit on the individual group
        clustering_output = f"{options.output_dir}/{group_file_name}_clustering"

        run_cd_hit(options, temp_fasta, clustering_output, clustering_mode)

        # Read the clustering results to find subgroups
        clustered_sequences = read_cd_hit_output(clustering_output + '.clstr')

        if len(clustered_sequences) == 1:
            # Detect if all sequences are identical in length and percentage identity
            all_same = check_if_all_identical(clustered_sequences)

        # **Global subgroup counter for the entire major group**
        subgroup_id = 0


        if not all_same:
            # Iterate through each cluster in clustered_sequences
            for cluster_key, cluster in clustered_sequences.items():

                remaining_sequences_tmp = sequences.copy()  # Track unprocessed sequences
                remaining_sequences = [entry for entry in remaining_sequences_tmp if entry[0] in
                                       {seq_entry['header'] for seq_entry in cluster}]
                sequences_to_remove = []

                while remaining_sequences:
                    # Track subgroups for this cluster pass
                    subgroup_sequences = []
                    genome_seen = set()

                    # Recalculate representative sequence dynamically for this cluster
                    rep = calculate_new_rep_seq(
                        [entry for entry in cluster if entry['header'] in (h for h, _ in remaining_sequences)]
                    )

                    # Find the sequence corresponding to rep['header'] from the list of sequences
                    rep_seq = next((seq for header, seq in sequences if header == rep['header']), None)

                    # Save previously checked seqs, so we don't have to compare them again.
                    checked = collections.defaultdict(float)

                    # Process each genome to select the best matching sequence
                    for genome in genome_to_gene_count:
                        best_sequence = None
                        best_score = None  # Initialise with a very low score, so that even negative scores can be selected

                        # Iterate over each sequence in the remaining sequences for this genome
                        for header, seq in remaining_sequences:
                            genome_id = header.split('|')[1]

                            if genome_id == genome:  # Ensure this sequence belongs to the current genome
                                if rep_seq == seq:
                                    levenshtein_distance = 0
                                else:
                                    if seq in checked:
                                        levenshtein_distance = checked[seq]
                                    else:
                                        levenshtein_distance = levenshtein_distance_calc(rep_seq,seq)
                                        checked[seq] = levenshtein_distance
                                # Lower Levenshtein distance means more 'similar' sequences
                                score = levenshtein_distance

                                # Check if this sequence has a higher score than the current best
                                if best_sequence == None:
                                    best_score = score
                                    best_sequence = (header, seq)  # Store the best matching sequence for this genome
                                elif score < best_score:
                                    best_score = score
                                    best_sequence = (header, seq)  # Store the best matching sequence for this genome

                        # Add the best sequence for this genome to the subgroup
                        if best_sequence is not None:
                            new_header = f">{group_file_name}_subgroup_{subgroup_id}|{best_sequence[0].split('|')[1]}|{best_sequence[0].split('|')[2]}"
                            subgroup_sequences.append((new_header, best_sequence[1]))
                            sequences_to_remove.append(best_sequence)
                            genome_seen.add(genome)

                    # Write each subgroup into a separate FASTA file
                    if subgroup_sequences:
                        subgroup_file = f"{options.output_dir}/{group_file_name}_subgroup_{subgroup_id}.fasta"
                        write_fasta(subgroup_sequences, subgroup_file)

                        # Remove processed sequences from the remaining list
                        remaining_sequences = [item for item in remaining_sequences if
                                               item[0] not in {h for h, _ in sequences_to_remove}]

                        # Increment subgroup ID for the next subgroup
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
    parser = argparse.ArgumentParser(description='PyamilySeq ' + PyamilySeq_Version + ': Group-Splitter - A tool to split multi-copy gene groups identified by PyamilySeq.')
    ### Required Arguments
    required = parser.add_argument_group('Required Parameters')
    required.add_argument('-input_fasta', action='store', dest='input_fasta',
                          help='Input FASTA file containing gene groups.',
                          required=True)
    required.add_argument('-sequence_type', action='store', dest='sequence_type', default='DNA',choices=['AA', 'DNA'],
                          help='Default - DNA: Are groups "DNA" or "AA" sequences?',
                          required=True)
    required.add_argument('-genome_num', action='store', dest='genome_num', type=int,
                          help='The total number of genomes must be provide',
                          required=True)
    required.add_argument('-output_dir', action='store', dest='output_dir',
                          help='Output directory.',
                          required=True)

    regrouping_params = parser.add_argument_group('Regrouping Parameters')
    regrouping_params.add_argument('-groups', action="store", dest='groups', default=None,
                          help='Default - auto: Detect groups to be split (see -group_threshold). '
                               'Provide "-groups 1,2,3,4" with group IDs to split specific groups.',
                          required=False)
    regrouping_params.add_argument('-group_threshold', action='store', dest='group_threshold', type=float, default=80,
                          help='Minimum percentage of genomes with multi-copy (default: 80.0) - Does not work with "-groups"')

    cdhit_params = parser.add_argument_group('CD-HIT Reclustering Parameters')
    cdhit_params.add_argument('-c', action='store', dest='pident', type=float, default=0.8,
                          help='Sequence identity threshold (default: 0.8) - Probably should be higher than what was used in initial clustering.')
    cdhit_params.add_argument('-s', action='store', dest='len_diff', type=float, default=0.20,
                          help="Length difference cutoff (default: 0.20) - Often the most impactful parameter to split 'multi-copy' gene groups.")
    cdhit_params.add_argument('-T', action='store', dest='clustering_threads', type=int, default=4,
                          help='Number of threads for clustering (default: 4)')
    cdhit_params.add_argument('-M', action='store', dest='clustering_memory', type=int, default=2000,
                          help='Memory limit in MB for clustering (default: 2000)')


    misc = parser.add_argument_group("Misc Parameters")
    misc.add_argument('-no_delete_temp_files', action='store_false', dest='delete_temp_files',
                      help='Default: Delete all temporary files after processing.',
                      required=False)
    misc.add_argument("-verbose", action="store_true", dest="verbose" ,
                      help="Print verbose output.",
                      required=False)
    misc.add_argument("-v", "--version", action="version",
                      version=f"PyamilySeq: Group-Splitter version {PyamilySeq_Version} - Exiting",
                      help="Print out version number and exit")


    options = parser.parse_args()
    print("Running PyamilySeq: Group-Splitter " + PyamilySeq_Version)



    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)

    if options.sequence_type == 'DNA':
        clustering_mode = 'cd-hit-est'
    else:
        clustering_mode = 'cd-hit'

    separate_groups(options, clustering_mode)


if __name__ == "__main__":

    main()
