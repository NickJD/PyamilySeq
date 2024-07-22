#from line_profiler_pycharm import profile

from collections import OrderedDict,defaultdict
import copy
import math
import sys
from tempfile import NamedTemporaryFile



try:
    from .Constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from Constants import *
    from utils import *


def custom_sort_key(k, dict1, dict2):
    return (len(dict1[k]), len(dict2[k]))

def sort_keys_by_values(dict1, dict2):
    sorted_keys = sorted(dict1.keys(), key=lambda k: custom_sort_key(k, dict1, dict2), reverse=True)
    return sorted_keys

def select_longest_gene(sequences):
    """Select the longest sequence for each genome."""
    longest_sequences = {}
    for seq_id, sequence in sequences.items():
        genome = seq_id.split('|')[0]  # Assuming genome name can be derived from the sequence ID
        if genome not in longest_sequences or len(sequence) > len(longest_sequences[genome][1]):
            longest_sequences[genome] = (seq_id, sequence)
    return longest_sequences


def run_mafft_on_sequences(options, sequences, output_file):
    """Run mafft on the given sequences and write to output file."""
    # Create a temporary input file for mafft
    with NamedTemporaryFile('w', delete=False) as temp_input_file:
        for header, sequence in sequences.items():
            temp_input_file.write(f">{header}\n{sequence}\n")
        temp_input_file_path = temp_input_file.name

    # Run mafft
    try:
        with open(output_file, 'w') as output_f:
            if options.verbose == 'True':
                subprocess.run(
                    ['mafft', '--auto', temp_input_file_path],
                    stdout=output_f,
                    stderr=sys.stderr,
                    check=True
                )
            else:
                subprocess.run(
                    ['mafft', '--auto', temp_input_file_path],
                    stdout=output_f,
                    stderr=subprocess.DEVNULL,  # Suppress stderr
                    check=True
                )
    finally:
        os.remove(temp_input_file_path)  # Clean up the temporary file


def process_gene_families(options, directory, output_file):
    """Process each gene family file to select the longest sequence per genome and concatenate aligned sequences."""
    concatenated_sequences = {}
    output_file = directory.replace('Gene_Families_Output',output_file)

    # Iterate over each gene family file
    for gene_file in os.listdir(directory):
        if gene_file.endswith('.fasta'):
            gene_path = os.path.join(directory, gene_file)

            # Read sequences from the gene family file
            sequences = read_fasta(gene_path)

            # Select the longest sequence for each genome
            longest_sequences = select_longest_gene(sequences)

            # Run mafft on the longest sequences
            aligned_file = f"{gene_file}_aligned.fasta"
            run_mafft_on_sequences(options, {seq_id: seq for seq_id, seq in longest_sequences.values()}, aligned_file)

            # Read aligned sequences and concatenate them
            aligned_sequences = read_fasta(aligned_file)
            for genome, aligned_seq in aligned_sequences.items():
                genome_name = genome.split('|')[0]
                if genome_name not in concatenated_sequences:
                    concatenated_sequences[genome_name] = ""
                concatenated_sequences[genome_name] += aligned_seq

            # Clean up aligned file
            os.remove(aligned_file)

    # Write the concatenated sequences to the output file
    with open(output_file, 'w') as out:
        for genome, sequence in concatenated_sequences.items():
            out.write(f">{genome}\n")
            wrapped_sequence = wrap_sequence(sequence, 60)
            out.write(f"{wrapped_sequence}\n")

def gene_presence_absence_output(options, genome_dict, pangenome_clusters_First_sorted, pangenome_clusters_First_sequences_sorted):
    print("Outputting gene_presence_absence file")
    in_name = options.clusters.split('.')[0]
    gpa_outfile = open(in_name+'_gene_presence_absence.csv','w')
    gpa_outfile.write('"Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Genome Fragment","Order within Fragment","'
                     '"Accessory Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc","')
    gpa_outfile.write('","'.join(genome_dict.keys()))
    gpa_outfile.write('"\n')
    for cluster, sequences in pangenome_clusters_First_sequences_sorted.items():
        average_sequences_per_genome = len(sequences) / len(pangenome_clusters_First_sorted[cluster])
        gpa_outfile.write('"group_'+str(cluster)+'","","'+str(len(pangenome_clusters_First_sorted[cluster]))+'","'+str(len(sequences))+'","'+str(average_sequences_per_genome)+
                         '","","","","","","","","",""')

        full_out = ''
        for genome in genome_dict.keys():
            tmp_list = []
            for value in sequences:
                if value.split('|')[0] == genome:
                    tmp_list.append(value)
            if tmp_list:
                full_out += ',"'+''.join(tmp_list)+'"'
            gpa_outfile.write(full_out)
        gpa_outfile.write('\n')

### Below is some unfinished code
    # edge_list_outfile = open(in_name+'_edge_list.csv','w')
    # for cluster, sequences in pangenome_clusters_First_sequences_sorted.items():
    #     output = []
    #     for entry in sequences:
    #         # Split each entry at '|'
    #         genome, gene = entry.split('|')
    #         # Format the result as "gene  genome"
    #         output.append(f"{gene}\t{genome}")
    #     for line in output:
    #         edge_list_outfile.write(line + '\n')


def wrap_sequence(sequence, width=60):
    wrapped_sequence = []
    for i in range(0, len(sequence), width):
        wrapped_sequence.append(sequence[i:i + width])
    return "\n".join(wrapped_sequence)


def read_fasta(fasta_file):
    sequences = {}
    current_sequence = None
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                current_sequence = line[1:]
                sequences[current_sequence] = ''
            else:
                sequences[current_sequence] += line
    return sequences


def reorder_dict_by_keys(original_dict, sorted_keys):
    return {k: original_dict[k] for k in sorted_keys}

def get_cores(options,genome_dict):
    ##Calculate core groups
    groups = OrderedDict()
    cores = OrderedDict()
    prev_top = len(genome_dict)
    first = True
    for group in options.core_groups.split(','):
        calculated_floor = math.floor(int(group) / 100 * len(genome_dict))
        if first == False:
            # Ensure no overlap
            # if calculated_floor <= prev_top:
            #     calculated_floor = prev_top - 1

            groups[group] = (calculated_floor,prev_top)
        else:
            groups[group] = (calculated_floor, prev_top)
            first = False
        prev_top = calculated_floor
        first_core_group = 'first_core_' + group
        cores[first_core_group] = []
        if options.reclustered != None:
            extended_core_group = 'extended_core_' + group
            cores[extended_core_group] = []
            combined_core_group = 'combined_core_' + group
            cores[combined_core_group] = []
            second_core_group = 'second_core_' + group
            cores[second_core_group] = []
            only_second_core_group = 'only_second_core_' + group
            cores[only_second_core_group] = []
    return cores, groups

#@profile
def calc_First_only_core(cluster, pep_num, groups, cores):
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= pep_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['first_core_'+family_group].append(cluster)

#@profile
def calc_single_First_extended_Second_only_core(pep_num, groups, cores, second_num): # Count gene families extended with StORFs
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= pep_num+second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['extended_core_' + family_group].append(pep_num)


#@profile
def calc_multi_First_extended_Second_only_core(pep_num, groups, cores, second_num): # Count seperately those gene families extended with StORF_Reporter but combined >1 PEP
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= pep_num+second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['combined_core_' + family_group] += 1


#@profile
def calc_Second_only_core(groups, cores, second_num):
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['second_core_' + family_group] += 1

#@profile
def calc_only_Second_only_core(groups, cores, second_num): # only count the true storf onlies
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['only_second_core_' + family_group] += 1





#@profile
def combined_clustering_counting(options, pangenome_clusters_First, reps, combined_pangenome_clusters_First_Second_clustered):
    num_clustered_First = defaultdict(list)
    pangenome_clusters_Type = copy.deepcopy(pangenome_clusters_First)
    list_of_reps = list(reps.keys())
    for cluster, pep_genomes in pangenome_clusters_First.items():
        rep = list_of_reps[int(cluster)]  # get the rep of the current pep cluster
        Com_PEP_Genomes = 0
        Seconds = 0
        seen_Seconds = []
        added_Second_genomes = 0
        try:  # get the cluster from the storf clusters which contains this rep
            clustered_combined = combined_pangenome_clusters_First_Second_clustered[rep]  # Not true clusters - I put a PEP as key myself
            seen_clust_Genomes = []
            num_clustered_First[cluster].append(rep + '_' + str(len(pep_genomes)))
            for clust in clustered_combined:
                if options.sequence_tag not in clust:  # Not good enough at the moment
                    clust_Genome = clust.split('|')[0]
                    if clust_Genome not in seen_clust_Genomes:
                        seen_clust_Genomes.append(clust_Genome)
                        if clust_Genome not in pep_genomes:
                            Com_PEP_Genomes += 1
                    num_clustered_First[cluster].append(clust + '_' + str(reps[clust][1]))
                elif options.sequence_tag in clust:
                    Seconds += 1
                    clust_Genome = clust.split('|')[0]
                    if clust_Genome not in seen_Seconds:
                        seen_Seconds.append(clust_Genome)
                    if clust_Genome not in seen_clust_Genomes:
                        seen_clust_Genomes.append(clust_Genome)
                        if clust_Genome not in pep_genomes:
                            added_Second_genomes += 1
                else:
                    sys.exit("Error: looking for sequence_tag")

            size_of_pep_clusters = []
            peps = num_clustered_First[cluster]
            for pep in peps:
                pep = pep.rsplit('_', 1)
                size_of_pep_clusters.append(int(pep[1]))
            pangenome_clusters_Type[cluster] = [len(num_clustered_First[cluster]), sum(size_of_pep_clusters),
                                                size_of_pep_clusters, added_Second_genomes, Seconds, len(seen_Seconds)]

        except KeyError:
            ###Singleton
            num_pep_genomes = [len(pep_genomes)]
            pangenome_clusters_Type[cluster] = [1, len(pep_genomes), num_pep_genomes, added_Second_genomes, Seconds,
                                                len(seen_Seconds)]

    return pangenome_clusters_Type

#@profile
def single_clustering_counting(options, pangenome_clusters_First, reps):
    num_clustered_First = defaultdict(list)
    recorded_First = []
    pangenome_clusters_Type = copy.deepcopy(pangenome_clusters_First)
    list_of_reps = list(reps.keys())
    for cluster, First_genomes in pangenome_clusters_First.items():
        rep = list_of_reps[int(cluster)]  # get the rep of the current pep cluster

        try:  # get the cluster from the storf clusters which contains this rep
            num_clustered_First[cluster].append(rep + '_' + str(len(First_genomes)))
            size_of_First_clusters = []
            Firsts = num_clustered_First[cluster]
            for First in Firsts:
                First = First.rsplit('_', 1)
                size_of_First_clusters.append(int(First[1]))
                recorded_First.append(First[0])
            pangenome_clusters_Type[cluster] = [len(num_clustered_First[cluster]), sum(size_of_First_clusters),
                                                size_of_First_clusters, 0, 0, 0]

        except KeyError:
            ###Singleton
            num_pep_genomes = [len(First_genomes)]
            pangenome_clusters_Type[cluster] = [1, len(First_genomes), num_pep_genomes, 0, 0, 0]

    return pangenome_clusters_Type



#@profile
def combined_clustering_CDHIT(options, genome_dict):
    unique_genomes = []
    Second_in = open(options.reclustered, 'r')
    combined_pangenome_clusters_First = OrderedDict()
    combined_pangenome_clusters_First_sequences = OrderedDict()
    combined_pangenome_clusters_Second = OrderedDict()
    combined_pangenome_clusters_Second_sequences = OrderedDict()
    combined_pangenome_clusters_First_Second_clustered = OrderedDict()

    not_Second_only_cluster_ids = []
    already_seen_PEP = []
    Combined_clusters = OrderedDict()
    Combined_reps = OrderedDict()
    first = True
    for line in Second_in:
        if line.startswith('>'):
            if first == False:
                cluster_size = len(Combined_clusters[cluster_id])
                Combined_reps.update({rep: cluster_size})
                for pep in combined_pangenome_clusters_First_sequences[cluster_id]:
                    if pep != []:
                        if pep in already_seen_PEP:
                            continue
                        else:
                            already_seen_PEP.append(pep)
                if len(combined_pangenome_clusters_Second_sequences[cluster_id]) > 0 and len(combined_pangenome_clusters_First_sequences[cluster_id]) > 0:
                    if len(combined_pangenome_clusters_First_sequences[cluster_id]) > 1:  # If we have clustered >1 PEP family, we need to record 1 as key and all others are val
                        all_but_first = combined_pangenome_clusters_First_sequences[cluster_id][1:]
                        storfs_clustered = combined_pangenome_clusters_Second_sequences[cluster_id]
                        VALUE = all_but_first + storfs_clustered
                    else:
                        VALUE = combined_pangenome_clusters_Second_sequences[cluster_id]
                    KEY = combined_pangenome_clusters_First_sequences[cluster_id][0]
                    combined_pangenome_clusters_First_Second_clustered.update({KEY: VALUE})
            cluster_id = line.strip('>')
            cluster_id = cluster_id.strip('\n')
            cluster_id = cluster_id.split(' ')[1]
            Combined_clusters.update({cluster_id: []})
            combined_pangenome_clusters_First.update({cluster_id: []})
            combined_pangenome_clusters_First_sequences.update({cluster_id: []})
            combined_pangenome_clusters_Second.update({cluster_id: []})
            combined_pangenome_clusters_Second_sequences.update({cluster_id: []})

            first = False
        else:
            clustered = line.split('\t')[1]
            clustered = clustered.split('>')[1]
            clustered = clustered.split('...')[0]
            genome = clustered.split('|')[0]
            genome_dict[genome] += 1
            if '*' in line:
                rep = clustered
                Combined_reps.update({rep: 0})
            if first == False:
                Combined_clusters[cluster_id].append(clustered)
                clustered_genome = clustered.split('|')[0]
                if options.sequence_tag in line:
                    if clustered_genome not in combined_pangenome_clusters_Second[cluster_id]:
                        combined_pangenome_clusters_Second[cluster_id].append(clustered_genome)
                    combined_pangenome_clusters_Second_sequences[cluster_id].append(clustered)
                else:
                    if cluster_id not in not_Second_only_cluster_ids:
                        not_Second_only_cluster_ids.append(cluster_id)  # Tell us which StORF_Reporter clustered are unmatched to a PEP
                    if clustered_genome not in combined_pangenome_clusters_First[cluster_id]:
                        combined_pangenome_clusters_First[cluster_id].append(clustered_genome)
                    combined_pangenome_clusters_First_sequences[cluster_id].append(clustered)


    return combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids, combined_pangenome_clusters_Second, unique_genomes

def combined_clustering_Edge_List(options, genome_dict):
    if options.cluster_format == 'TSV':
        separator = '\t'
    elif options.cluster_format == 'CSV':
        separator = ','
    unique_genomes = []
    cluster_id = 0
    last_rep = ''
    Second_in = open(options.reclustered, 'r')
    combined_pangenome_clusters_First = OrderedDict()
    combined_pangenome_clusters_First_sequences = OrderedDict()
    combined_pangenome_clusters_Second = OrderedDict()
    combined_pangenome_clusters_Second_sequences = OrderedDict()
    combined_pangenome_clusters_First_Second_clustered = OrderedDict()

    not_Second_only_cluster_ids = []
    already_seen_PEP = []
    Combined_clusters = OrderedDict()
    Combined_reps = OrderedDict()
    first = True
    for line in Second_in:
        rep, child = line.strip().split(separator)
        child_genome = child.split('|')[0]  # Extracting the genome identifier from the child sequence

        if first == True:
            Combined_clusters.update({cluster_id: []})
            combined_pangenome_clusters_First.update({cluster_id: []})
            combined_pangenome_clusters_First_sequences.update({cluster_id: []})
            combined_pangenome_clusters_Second.update({cluster_id: []})
            combined_pangenome_clusters_Second_sequences.update({cluster_id: []})
            Combined_reps.update({rep: 0})
            first = False

        if first == False:
            if rep != last_rep and last_rep != '':
                cluster_size = len(Combined_clusters[cluster_id])
                Combined_reps.update({rep: cluster_size})
                for pep in combined_pangenome_clusters_First_sequences[cluster_id]:
                    if pep != []:
                        if pep in already_seen_PEP:
                            continue
                        else:
                            already_seen_PEP.append(pep)
                if len(combined_pangenome_clusters_Second_sequences[cluster_id]) > 0 and len(combined_pangenome_clusters_First_sequences[cluster_id]) > 0:
                    if len(combined_pangenome_clusters_First_sequences[cluster_id]) > 1:  # If we have clustered >1 PEP family, we need to record 1 as key and all others are val
                        all_but_first = combined_pangenome_clusters_First_sequences[cluster_id][1:]
                        storfs_clustered = combined_pangenome_clusters_Second_sequences[cluster_id]
                        VALUE = all_but_first + storfs_clustered
                    else:
                        VALUE = combined_pangenome_clusters_Second_sequences[cluster_id]
                    KEY = combined_pangenome_clusters_First_sequences[cluster_id][0]
                    combined_pangenome_clusters_First_Second_clustered.update({KEY: VALUE})

                cluster_id += 1
                Combined_clusters.update({cluster_id: []})
                combined_pangenome_clusters_First.update({cluster_id: []})
                combined_pangenome_clusters_First_sequences.update({cluster_id: []})
                combined_pangenome_clusters_Second.update({cluster_id: []})
                combined_pangenome_clusters_Second_sequences.update({cluster_id: []})
                Combined_reps.update({rep: 0})


        Combined_clusters[cluster_id].append(child)
        if options.sequence_tag in line:
            if child_genome not in combined_pangenome_clusters_Second[cluster_id]:
                combined_pangenome_clusters_Second[cluster_id].append(child_genome)
            combined_pangenome_clusters_Second_sequences[cluster_id].append(child)
        else:
            if cluster_id not in not_Second_only_cluster_ids:
                not_Second_only_cluster_ids.append(cluster_id)  # Tell us which StORF_Reporter clustered are unmatched to a PEP
            if child_genome not in combined_pangenome_clusters_First[cluster_id]:
                combined_pangenome_clusters_First[cluster_id].append(child_genome)
            combined_pangenome_clusters_First_sequences[cluster_id].append(child)

        last_rep = rep

    return combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids, combined_pangenome_clusters_Second, unique_genomes


def cluster_EdgeList(options):
    if options.cluster_format == 'TSV':
        separator = '\t'
    elif options.cluster_format == 'CSV':
        separator = ','
    cluster_id = 0
    last_rep = ''
    first = True
    First_in = open(options.clusters, 'r')
    pangenome_clusters_First = OrderedDict()
    pangenome_clusters_First_sequences = OrderedDict()
    genome_dict = defaultdict(int)
    reps = OrderedDict()
    for line in First_in:
        rep, child = line.strip().split(separator)
        child_genome = child.split('|')[0]  # Extracting the genome identifier from the child sequence
        # Counting occurrences of genomes
        genome_dict[child_genome] += 1
        if first == True:
            pangenome_clusters_First[0] = []
            pangenome_clusters_First_sequences[0] = []
            first = False

        if rep != last_rep and last_rep != '':
            cluster_id +=1
            pangenome_clusters_First[cluster_id] = []
            pangenome_clusters_First_sequences[cluster_id] = []
            cluster_size = len(pangenome_clusters_First_sequences[cluster_id-1])
            reps.update({last_rep: [cluster_size, len(pangenome_clusters_First[cluster_id-1])]})
            pangenome_clusters_First[cluster_id] = []
            pangenome_clusters_First_sequences[cluster_id] = []
        if child_genome not in pangenome_clusters_First[cluster_id]:
            pangenome_clusters_First[cluster_id].append(child_genome)

        pangenome_clusters_First_sequences[cluster_id].append(child)
        last_rep = rep
        cluster_size = len(pangenome_clusters_First_sequences[cluster_id])
        reps.update({rep: [cluster_size, len(pangenome_clusters_First[cluster_id])]})


    return genome_dict, pangenome_clusters_First, pangenome_clusters_First_sequences, reps



def cluster_CDHIT(options):
    First_in = open(options.clusters, 'r')
    clusters = OrderedDict()
    pangenome_clusters_First = OrderedDict()
    pangenome_clusters_First_sequences = OrderedDict()
    first = True
    genome_dict = defaultdict(int)
    reps = OrderedDict()
    ## Load in all data for easier reuse later
    for line in First_in:
        if line.startswith('>'):
            if first == False:
                cluster_size = len(clusters[cluster_id])
                reps.update({rep: [cluster_size, len(pangenome_clusters_First[cluster_id])]})
            cluster_id = line.strip('>')
            cluster_id = cluster_id.strip('\n')
            cluster_id = cluster_id.split(' ')[1]
            clusters.update({cluster_id: []})
            pangenome_clusters_First.update({cluster_id: []})
            pangenome_clusters_First_sequences.update({cluster_id: []})

            first = False
        else:
            clustered = line.split('\t')[1]
            clustered = clustered.split('>')[1]
            clustered = clustered.split('...')[0]
            genome = clustered.split('|')[0]
            genome_dict[genome] += 1
            if '*' in line:
                rep = clustered
                reps.update({rep: [0, 0]})
            if first == False:
                clusters[cluster_id].append(clustered)
                clustered_genome = clustered.split('|')[0]
                if clustered_genome not in pangenome_clusters_First[cluster_id]:
                    pangenome_clusters_First[cluster_id].append(clustered_genome)
                pangenome_clusters_First_sequences[cluster_id].append(clustered)
    return genome_dict, pangenome_clusters_First, pangenome_clusters_First_sequences, reps

#@profile
def cluster(options):

    if options.cluster_format == 'CD-HIT':
        genome_dict, pangenome_clusters_First, pangenome_clusters_First_sequences, reps = cluster_CDHIT(options)
    elif options.cluster_format in ['TSV','CSV']:
        genome_dict, pangenome_clusters_First, pangenome_clusters_First_sequences, reps = cluster_EdgeList(options)

    ######################################
    cores, groups = get_cores(options, genome_dict)
    ###

    if options.reclustered != None:
        if options.cluster_format == 'CD-HIT':
            combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids,combined_pangenome_clusters_Second,\
                unique_genomes = combined_clustering_CDHIT(options, genome_dict)
        if options.cluster_format == 'TSV':
            combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids,combined_pangenome_clusters_Second,\
                unique_genomes = combined_clustering_Edge_List(options, genome_dict)
        pangenome_clusters_Type = combined_clustering_counting(options, pangenome_clusters_First, reps, combined_pangenome_clusters_First_Second_clustered)
    else:
        pangenome_clusters_Type = single_clustering_counting(options, pangenome_clusters_First, reps)



    Number_Of_Second_Extending_But_Same_Genomes = 0

    sorted_first_keys = sort_keys_by_values(pangenome_clusters_First, pangenome_clusters_First_sequences)
    pangenome_clusters_First_sorted = reorder_dict_by_keys(pangenome_clusters_First, sorted_first_keys)
    pangenome_clusters_First_sequences_sorted = reorder_dict_by_keys(pangenome_clusters_First_sequences, sorted_first_keys)
    pangenome_clusters_Type_sorted = reorder_dict_by_keys(pangenome_clusters_Type, sorted_first_keys)

    print("Calculating Groups")
    for cluster, numbers in pangenome_clusters_Type_sorted.items():
    ############################### Calculate First only
        #if numbers[0] == 1 and numbers[1] >=2:
        calc_First_only_core(cluster, numbers[1],groups,cores)

        # elif numbers[0] >1 and numbers[1] >=2:
        #     calc_First_only_core(cluster, numbers[2][0],groups,cores)


    if options.reclustered != None:
        ############################# Calculate First and Reclustered-Second
        if numbers[0] == 1 and numbers[3] >= 1:  # If Seconds did not combine First reps
            calc_single_First_extended_Second_only_core(numbers[1], groups, cores, numbers[3])
        elif numbers[0] > 1 and numbers[3] >= 1:  # If unique Secondss combined multiple Firsts
            calc_multi_First_extended_Second_only_core(numbers[1], groups, cores, numbers[3])
        elif numbers[4] >= 1:
            Number_Of_Second_Extending_But_Same_Genomes += 1
        combined_pangenome_clusters_ONLY_Second_Type = defaultdict(list)
        combined_pangenome_clusters_Second_Type = defaultdict(list)
        for cluster, genomes in combined_pangenome_clusters_Second.items():
            if cluster in not_Second_only_cluster_ids:
                combined_pangenome_clusters_Second_Type[cluster] = [cluster, len(genomes)]
            else:
                combined_pangenome_clusters_ONLY_Second_Type[cluster] = [cluster, len(genomes)]
        for cluster, data in combined_pangenome_clusters_Second_Type.items():
            calc_Second_only_core(groups, cores, data[1])
        for cluster, data in combined_pangenome_clusters_ONLY_Second_Type.items():
            if data[1] >= 2:
                calc_only_Second_only_core(groups, cores, data[1])
    ###########################
    key_order = ['first_core_', 'extended_core_', 'combined_core_', 'second_core_','only_second_core_']
    print("Gene Groups:")
    for key_prefix in key_order:
        for key, value in cores.items():
            if key.startswith(key_prefix):
                print(f"{key}: {len(value)}")
    print("Total Number of Gene Groups (Including Singletons): " + str(len(pangenome_clusters_First_sequences_sorted)))

    if options.gene_presence_absence_out != None:
        gene_presence_absence_output(options,genome_dict, pangenome_clusters_First_sorted, pangenome_clusters_First_sequences_sorted)

    if options.write_families != None and options.fasta != None:
        sequences = read_fasta(options.fasta)
        input_dir = os.path.dirname(os.path.abspath(options.clusters))
        output_dir = os.path.join(input_dir, 'Gene_Families_Output')

        # Create output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for key_prefix in key_order:
            for key, values in cores.items():
                if any(part in options.write_families.split(',') for part in key.split('_')):
                    if key.startswith(key_prefix):
                        for value in values:
                            output_filename = f"{key}_{value}.fasta"
                            sequences_to_write = pangenome_clusters_First_sequences_sorted[value]
                            # Write sequences to output file that are in the sequences dictionary
                            with open(os.path.join(output_dir, output_filename), 'w') as outfile:
                                for header in sequences_to_write:
                                    if header in sequences:
                                        outfile.write(f">{header}\n")
                                        wrapped_sequence = wrap_sequence(sequences[header])
                                        outfile.write(f"{wrapped_sequence}\n")

    if options.con_core != None and options.fasta != None and options.write_families != None:
        process_gene_families(options, os.path.join(input_dir, 'Gene_Families_Output'), 'concatonated_genes_aligned.fasta')


        # groups_dir = os.path.join(input_dir, 'Gene_Families_Output')
        # """Run mafft on all .fasta files in the given directory."""
        # for filename in os.listdir(groups_dir):
        #     if filename.endswith('.fasta'):
        #         input_path = os.path.join(groups_dir, filename)
        #         output_filename = filename.replace('.fasta', '_mafft.aln')
        #         output_path = os.path.join(groups_dir, output_filename)
        #
        #         # Call mafft command
        #         try:
        #             with open(output_path, 'w') as output_file:
        #                 subprocess.run(
        #                     ['mafft', '--auto', input_path],
        #                     stdout=output_file,
        #                     stderr=subprocess.DEVNULL,  # Suppress stderr
        #                     check=True
        #                 )
        #             print(f"Processed {input_path} -> {output_path}")
        #         except subprocess.CalledProcessError as e:
        #             print(f"Failed to process {input_path}: {e}")

        ##This could be run once and not above AND here..
        # output_dir = os.path.dirname(os.path.abspath(options.clusters))
        # sequences = read_fasta(options.fasta)
        # concatenated_sequences = {genome: '' for genome in genome_dict.keys()}
        #
        #
        # for key_prefix in key_order:
        #     for key, values in cores.items():
        #         if any(part in options.con_core.split(',') for part in key.split('_')):
        #             if key.startswith(key_prefix):
        #                 for value in values:
        #                     length_capture = {genome: [] for genome in genome_dict.keys()}
        #                     sequences_to_write = pangenome_clusters_First_sequences_sorted[value]
        #                     for header in sequences_to_write:
        #                         if header in sequences:
        #                             length_capture[header.split('|')[0]].append([header,len(sequences[header])])
        #                     if all(bool(values) for values in length_capture.values()): # If a GF is not present in 'ALL' genomes, do not add to concat
        #                         for genome, lengths in length_capture.items():
        #                             max_value = float('-inf')
        #                             max_item = None
        #                             for length in lengths:
        #                                 current_value = length[1]
        #                                 if current_value > max_value:
        #                                     max_value = current_value
        #                                     max_item = length[0]
        #                             concatenated_sequences[genome.split('|')[0]] += sequences[max_item]
        #
        #
        # with open(os.path.join(output_dir, 'core_concat.fasta'), 'w') as outfile:
        #     for genome, sequence in concatenated_sequences.items():
        #         outfile.write(f">{genome}\n")
        #         wrapped_sequence = wrap_sequence(sequence)
        #         outfile.write(f"{wrapped_sequence}\n")


        # for core_gene_family in core_gene_families:
        #     found_sequences = {genome: False for genome in genomes}
        #
        #     for fasta_file in fasta_files:f
        #         sequences = read_fasta(fasta_file)
        #         for header, sequence in sequences.items():
        #             genome = header.split('|')[0]
        #             if genome in genomes and core_gene_family in header:
        #                 concatenated_sequences[genome] += sequence
        #                 found_sequences[genome] = True
        #
        #     for genome in genomes:
        #         if not found_sequences[genome]:
        #             concatenated_sequences[genome] += '-' * len(next(iter(sequences.values())))




