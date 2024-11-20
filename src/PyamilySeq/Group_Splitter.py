
import argparse
from collections import defaultdict, OrderedDict


try:
    from .constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *
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
        '-g', str(options.fast_mode),
        '-sc', "1",
        '-sf', "1"
    ]
    if options.verbose == True:
        subprocess.run(cdhit_command)
    else:
        subprocess.run(cdhit_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

#'@profile
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



#def length_within_threshold(rep_length, length, len_diff):
#    return abs(rep_length - length) / rep_length <= len_diff


def check_if_all_identical(clustered_sequences):
    lengths = {entry['length'] for cluster in clustered_sequences.values() for entry in cluster}
    perc_idents = {entry['percent_identity'] for cluster in clustered_sequences.values() for entry in cluster}

    return len(lengths) == 1 and len(perc_idents) == 1



def read_fasta_groups(options, groups_to_use):
    groups = defaultdict(list)
    genome_count = defaultdict(int)
    current_group = None
    current_sequence = []

    if options.sequence_type == 'AA':
        affix = '_aa.fasta'
    else:
        affix = '_dna.fasta'

    combined_groups_fasta = options.input_directory + '/Gene_Groups_Output/combined_group_sequences' + affix

    if groups_to_use[0] == 'ids':
        selected_group_ids = [int(g.strip()) for g in groups_to_use[1].split(',')]
    elif groups_to_use[0] == 'groups':
        selected_groups = set(range(int(groups_to_use[1]), 101))
        # Scan the directory for filenames that match the criteria
        selected_group_ids = []
        for filename in os.listdir(os.path.dirname(combined_groups_fasta)):
            if 'core' in filename and filename.endswith('.fasta'):
                try:
                    group_number = int(filename.split('_')[2])
                    if group_number in selected_groups:
                        selected_group_ids.append(int(filename.split('_')[3].split('.')[0]))
                except ValueError:
                    continue


    group_number = None
    with open(combined_groups_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_group is not None and (selected_group_ids is None or group_number in selected_group_ids):
                    groups[current_group].append((current_group_header, ''.join(current_sequence)))

                current_group_header = line.strip()
                current_group = current_group_header.split('|')[0]
                genome = current_group_header.split('|')[1]
                current_sequence = []
                genome_count[genome] += 1

                # Only process if group matches the selected_groups or if no specific groups were provided
                group_number = int(current_group.replace('>Group_', ''))  # Assuming format 'Group_n'
                if selected_group_ids is not None and group_number not in selected_group_ids:
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

#@profile
def separate_groups(options, clustering_mode, groups_to_use):
    groups, genome_count = read_fasta_groups(options, groups_to_use)

    paralog_groups = defaultdict(lambda: {'count': 0, 'sizes': []})  # To track number of paralog groups and their sizes

    for group_header, sequences in groups.items():
        if options.verbose == True:
            print(f"\n###\nCurrent Group: {group_header.replace('>','')}\n")

        group_name = group_header.split('|')[0]  # Get the group part (e.g., '>Group_n')

        # Count genomes with more than one gene
        genome_to_gene_count = defaultdict(int)
        for header, _ in sequences:
            genome = header.split('|')[1]
            genome_to_gene_count[genome] += 1

        num_genomes_with_multiple_genes = sum(1 for count in genome_to_gene_count.values() if count > 1)

        # Check if the group meets the threshold for having paralogs
        #if options.groups == None:
        if (num_genomes_with_multiple_genes / options.genome_num) * 100 < options.group_threshold:
            continue


        group_file_name = group_name.replace('>','')

        temp_fasta = f"{options.gene_groups_output}/{group_file_name}.fasta"
        write_fasta(sequences, temp_fasta)

        # Run cd-hit on the individual group
        clustering_output = f"{options.gene_groups_output}/{group_file_name}_clustering"

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
                        subgroup_file = f"{options.sub_groups_output}/{group_file_name}_subgroup_{subgroup_id}.fasta"
                        write_fasta(subgroup_sequences, subgroup_file)

                        # Remove processed sequences from the remaining list
                        remaining_sequences = [item for item in remaining_sequences if
                                               item[0] not in {h for h, _ in sequences_to_remove}]

                        # Increment subgroup ID for the next subgroup
                        subgroup_id += 1
                        paralog_groups[group_name]['count'] += 1  # Count this group as a paralog group
                        paralog_groups[group_name]['sizes'].append(len(subgroup_sequences))  # Record the size of the subgroup




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
                subgroup_file = f"{options.input_directory}/{group_file_name}_subgroup_{subgroup_id}.fasta"
                write_fasta(seqs, subgroup_file)

                # Increment subgroup ID globally for the next subgroup
                subgroup_id += 1
                paralog_groups[group_name]['count'] += 1  # Count this group as a paralog group
                paralog_groups[group_name]['sizes'].append(len(seqs))  # Record the size of the subgroup



        # Clean up temporary fasta file if the option is set
        if options.delete_temp_files:
            if temp_fasta and os.path.exists(temp_fasta):
                os.remove(temp_fasta)
            if os.path.exists(clustering_output + '.clstr'):
                os.remove(clustering_output + '.clstr')
            if os.path.exists(clustering_output):
                os.remove(clustering_output)


    return paralog_groups


def main():
    parser = argparse.ArgumentParser(description='PyamilySeq ' + PyamilySeq_Version + ': Group-Splitter - A tool to split multi-copy gene groups identified by PyamilySeq.')
    ### Required Arguments
    required = parser.add_argument_group('Required Parameters')
    required.add_argument('-input_directory', action='store', dest='input_directory',
                          help='Provide the directory of a PyamilySeq run.',
                          required=True)
    required.add_argument('-sequence_type', action='store', dest='sequence_type', default='AA',choices=['AA', 'DNA'],
                          help='Default - AA: Are groups "DNA" or "AA" sequences?',
                          required=True)
    required.add_argument('-genome_num', action='store', dest='genome_num', type=int,
                          help='The total number of genomes must be provide',
                          required=True)


    ### Regrouping Arguments
    regrouping_params = parser.add_argument_group('Regrouping Parameters')
    regrouping_params.add_argument('-groups', action="store", dest='groups', type=int, default=None,
                          help='Default - 99: groups to be split by pangenome grouping (see -group_threshold). '
                               'Provide "-groups 99" to split specific groups.',
                          required=False)
    regrouping_params.add_argument('-group_ids', action="store", dest='group_ids', default=None,
                          help='Default - None: Provide "-group_ids 1,2,3,4" to split specific groups (see -group_threshold).',
                          required=False)
    regrouping_params.add_argument('-group_threshold', action='store', dest='group_threshold', type=float, default=80,
                          help='Default: 80: Minimum percentage of genomes with multi-copy in a gene group to be split.',
                          required=False)

    ### Output Arguments
    output_args = parser.add_argument_group('Output Parameters')
    output_args.add_argument('-a', action="store", dest='align_core', default=None,
                          help='Default - No output: SLOW! (Only works for Species mode) Output aligned and concatinated sequences of identified groups -'
                               'provide group levels at which to output "-a 99,95".',
                          required=False)


    ### CD-HIT Reclustering Arguments
    cdhit_params = parser.add_argument_group('CD-HIT Reclustering Parameters')
    cdhit_params.add_argument('-c', action='store', dest='pident', type=float, default=0.8,
                          help='Sequence identity threshold (default: 0.8) - Probably should be higher than what was used in initial clustering.')
    cdhit_params.add_argument('-s', action='store', dest='len_diff', type=float, default=0.20,
                          help="Length difference cutoff (default: 0.20) - Often the most impactful parameter to split 'multi-copy' gene groups.")
    cdhit_params.add_argument('-fastmode', action='store_true', dest='fast_mode',
                      help='Default False: Run CD-HIT with "-g 0" to speed up but reduce accuracy of clustering.',
                      required=False)
    cdhit_params.add_argument('-T', action='store', dest='clustering_threads', type=int, default=8,
                          help='Number of threads for clustering (default: 8)')
    cdhit_params.add_argument('-M', action='store', dest='clustering_memory', type=int, default=2000,
                          help='Memory limit in MB for clustering (default: 2000)')

    ### MAFFT Alignment Arguments
    alignment_args = parser.add_argument_group('Alignment Runtime Arguments - Optional when "-a" is provided.')

    alignment_args.add_argument("-t", action="store", dest="threads", type=int, default=8,
                          help="Default 8: Threads to be allocated for clustering and/or alignment.",
                          required=False)

    ### Misc Arguments
    misc = parser.add_argument_group("Misc Parameters")
    misc.add_argument('-no_delete_temp_files', action='store_false', dest='delete_temp_files',
                      help='Default: Delete all temporary files after processing.',
                      required=False)
    misc.add_argument("-verbose", action="store_true", dest="verbose",
                      help="Print verbose output.",
                      required=False)
    misc.add_argument("-v", "--version", action="version",
                      version=f"PyamilySeq: Group-Splitter version {PyamilySeq_Version} - Exiting",
                      help="Print out version number and exit")


    options = parser.parse_args()
    print("Running PyamilySeq: Group-Splitter " + PyamilySeq_Version)



    ###External tool checks:
    ##MAFFT
    if options.align_core == True:
        if is_tool_installed('mafft'):
            if options.verbose == True:
                print("mafft is installed. Proceeding with alignment.")
        else:
            exit("mafft is not installed. Please install mafft to proceed.")
    ##CD-HIT

    if is_tool_installed('cd-hit'):
        if options.verbose == True:
            print("cd-hit is installed. Proceeding with clustering.")
        if options.sequence_type == 'DNA':
            clustering_mode = 'cd-hit-est'
        else:
            clustering_mode = 'cd-hit'
        if options.fast_mode == True:
            options.fast_mode = 0
            if options.verbose == True:
                print("Running CD-HIT in fast mode.")
        else:
            options.fast_mode = 1
            if options.verbose == True:
                print("Running CD-HIT in slow mode.")
    else:
        exit("cd-hit is not installed. Please install cd-hit to proceed.")

    ##Alignment
    if options.align_core != None:
        if options.groups == None and options.group_ids == None:
            sys.exit('Must provide "-groups" or "-group_ids" when requesting alignment with "-a".')

    ##Output Directories
    gene_groups_output = os.path.join(options.input_directory, "Gene_Groups_Output")
    options.gene_groups_output = gene_groups_output
    sub_groups_output = os.path.join(options.input_directory, "Sub_Groups_Output")
    options.sub_groups_output = sub_groups_output
    if not os.path.exists(gene_groups_output):
        os.makedirs(gene_groups_output)
    if not os.path.exists(sub_groups_output):
        os.makedirs(sub_groups_output)

    ## Get Summary Stats
    summary_file = os.path.join(options.input_directory, 'summary_statistics.txt')

    # Save arguments to a text file
    params_out = os.path.join(options.input_directory, 'Group-Splitter_params.txt')
    with open(params_out, "w") as outfile:
        for arg, value in vars(options).items():
            outfile.write(f"{arg}: {value}\n")



    ## Group Selction - FIX THIS - currently fails if either are not provided
    if options.groups != None and options.group_ids != None:
        sys.exit('Must provide "-group_ids" or "-groups", not both.')
    elif options.group_ids != None:
        groups_to_use = ('ids', options.group_ids)
    elif options.groups != None:
        groups_to_use = ('groups', options.groups)
    else:
        groups_to_use = ('groups', 99)



    paralog_groups = separate_groups(options, clustering_mode, groups_to_use)
    ###
    # Print metrics about paralog groups
    print(f"Identified {len(paralog_groups)} paralog groups:")
    for group_id, data in paralog_groups.items():
        print(f"Group ID: {group_id}, Number of new groups: {data['count']}, Sizes: {data['sizes']}")
    ###


    # Read summary statistics
    with open(summary_file, 'r') as f:
        summary_data = f.read().splitlines()

    summary_info = {}
    for line in summary_data:
        if ':' in line:
            key, value = line.split(':')
            summary_info[key.strip()] = int(value.strip())

    genome_num = summary_info['Number of Genomes']
    core_99 = summary_info['First_core_99']
    core_95 = summary_info['First_core_95']
    core_15 = summary_info['First_core_15']
    core_0 = summary_info['First_core_0']
    total_gene_groups = summary_info['Total Number of First Gene Groups (Including Singletons)']

    # Initialise new core values
    new_core_99 = core_99
    new_core_95 = core_95
    new_core_15 = core_15
    new_core_0 = core_0

    # Recalculate each *_core_* value
    for group_id, data in paralog_groups.items():
        group_id = group_id.replace('>Group_', '')
        original_group = next((f for f in os.listdir(gene_groups_output) if f.endswith(f'_{group_id}.fasta')), None)
        original_group = int(original_group.split('_')[2])
        if original_group == 99:
            new_core_99 -= 1
        elif original_group == 95:
            new_core_95 -= 1
        elif original_group == 15:
            new_core_15 -= 1
        elif original_group == 0:
            new_core_0 -= 1

        for size in data['sizes']:
            if size >= math.floor(99 * genome_num / 100):
                new_core_99 += 1
            elif size >= math.floor(95 * genome_num / 100):
                new_core_95 += 1
            elif size >= math.floor(15 * genome_num / 100):
                new_core_15 += 1
            elif size >= math.floor(0 * genome_num / 100):
                new_core_0 += 1




    # Write out the new summary statistics - currently only works for default cores
    stats_out = summary_file.replace('.txt','_recalculated.txt')
    key_order = ['First_core_', 'extended_core_', 'combined_core_', 'Second_core_','only_Second_core_']
    with open(stats_out, 'w') as outfile:
        print("Number of Genomes: " + str(options.genome_num))
        outfile.write("Number of Genomes: " + str(options.genome_num) + "\n")
        print("Reclaculated Gene Groups:")
        outfile.write("Recalculated Gene Groups\n")
        print(f"First_core_99: {new_core_99}")
        outfile.write(f"First_core_99: {new_core_99}\n")
        print(f"First_core_95: {new_core_95}")
        outfile.write(f"First_core_95: {new_core_95}\n")
        print(f"First_core_15: {new_core_15}")
        outfile.write(f"First_core_15: {new_core_15}\n")
        print(f"First_core_0: {new_core_0}")
        outfile.write(f"First_core_0: {new_core_0}\n")
        print("Total Number of First Gene Groups (Including Singletons): " + str(total_gene_groups))
        outfile.write("Total Number of First Gene Groups (Including Singletons): " + str(total_gene_groups))

    # Alignment
    if options.align_core != None:
        print("\n\nProcessing gene group alignment")
        group_directory = options.gene_groups_output
        sub_group_directory = options.sub_groups_output
        genome_list = read_genomes_from_fasta(options.gene_groups_output + '/combined_group_sequences_dna.fasta')
        process_gene_groups(options, group_directory, sub_group_directory, paralog_groups, genome_list, 'concatenated_genes_post_splitting_aligned_dna.fasta')







if __name__ == "__main__":

    main()
