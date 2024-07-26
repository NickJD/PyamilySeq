#from line_profiler_pycharm import profile

import copy
import math
import sys



try:
    from .Constants import *
    from .clusterings import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from Constants import *
    from clusterings import *
    from utils import *


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
    output_dir = os.path.abspath(options.output_dir)
    in_name = options.clusters.split('.')[0].split('/')[-1]
    gpa_outfile = os.path.join(output_dir, in_name)
    gpa_outfile = open(gpa_outfile+'_gene_presence_absence.csv','w')
    gpa_outfile.write('"Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Genome Fragment","Order within Fragment","'
                     '"Accessory Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc","')
    gpa_outfile.write('","'.join(genome_dict.keys()))
    gpa_outfile.write('"\n')
    for cluster, sequences in pangenome_clusters_First_sequences_sorted.items():
        average_sequences_per_genome = len(sequences) / len(pangenome_clusters_First_sorted[cluster])
        gpa_outfile.write('"group_'+str(cluster)+'","","'+str(len(pangenome_clusters_First_sorted[cluster]))+'","'+str(len(sequences))+'","'+str(average_sequences_per_genome)+
                         '","","","","","","","","",""')


        for genome in genome_dict.keys():
            full_out = ''
            tmp_list = []
            for value in sequences:
                if value.split('|')[0] == genome:
                    tmp_list.append(value)
            if tmp_list:
                full_out += ',"'+''.join(tmp_list)+'"'
            else:
                full_out = ',""'
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




def get_cores(options,genome_dict):
    ##Calculate core groups
    groups = OrderedDict()
    cores = OrderedDict()
    prev_top = len(genome_dict)
    first = True
    for group in options.core_groups.split(','):
        calculated_floor = math.floor(int(group) / 100 * len(genome_dict))
        if first == False:
            groups[group] = (calculated_floor,prev_top)
        else:
            groups[group] = (calculated_floor, prev_top)
            first = False
        prev_top = calculated_floor
        first_core_group = 'First_core_' + group
        cores[first_core_group] = []
        if options.reclustered != None:
            extended_core_group = 'extended_core_' + group
            cores[extended_core_group] = []
            combined_core_group = 'combined_core_' + group
            cores[combined_core_group] = []
            second_core_group = 'Second_core_' + group
            cores[second_core_group] = []
            only_second_core_group = 'only_Second_core_' + group
            cores[only_second_core_group] = []
    return cores, groups

#@profile
def calc_First_only_core(cluster, First_num, groups, cores):
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= First_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['First_core_'+family_group].append(cluster)

#@profile
def calc_single_First_extended_Second_only_core(cluster, First_num, groups, cores, Second_num): # Count gene families extended with StORFs
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= First_num+Second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['extended_core_' + family_group].append(cluster)


#@profile
def calc_multi_First_extended_Second_only_core(cluster, First_num, groups, cores, Second_num): # Count seperately those gene families extended with StORF_Reporter but combined >1 PEP
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= First_num+Second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['combined_core_' + family_group].append(cluster)


#@profile
def calc_Second_only_core(cluster, Second_num, groups, cores):
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= Second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['Second_core_' + family_group].append(cluster)

#@profile
def calc_only_Second_only_core(cluster, Second_num, groups, cores): # only count the true storf onlies
    groups_as_list = list(groups.values())
    for idx in (idx for idx, (sec, fir) in enumerate(groups_as_list) if sec <= Second_num <= fir):
        res = idx
    family_group = list(groups)[res]
    cores['only_Second_core_' + family_group].append(cluster)



#@profile
def cluster(options):

    if options.cluster_format == 'CD-HIT':
        genome_dict, pangenome_clusters_First, pangenome_clusters_First_sequences, reps = cluster_CDHIT(options, '|')
    elif options.cluster_format in ['TSV','CSV']:
        genome_dict, pangenome_clusters_First, pangenome_clusters_First_sequences, reps = cluster_EdgeList(options, '|')

    ###
    cores, groups = get_cores(options, genome_dict)
    ###

    if options.reclustered != None:
        if options.cluster_format == 'CD-HIT':
            combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids,combined_pangenome_clusters_Second = combined_clustering_CDHIT(options, genome_dict, '|')
        if options.cluster_format == ['TSV','CSV']:
            combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids,combined_pangenome_clusters_Second = combined_clustering_Edge_List(options, '|')

        pangenome_clusters_Type = combined_clustering_counting(options, pangenome_clusters_First, reps, combined_pangenome_clusters_First_Second_clustered, '|')
    else:
        pangenome_clusters_Type = single_clustering_counting(pangenome_clusters_First, reps)



    Number_Of_Second_Extending_But_Same_Genomes = 0

    sorted_first_keys = sort_keys_by_values(pangenome_clusters_First, pangenome_clusters_First_sequences)
    pangenome_clusters_First_sorted = reorder_dict_by_keys(pangenome_clusters_First, sorted_first_keys)
    pangenome_clusters_First_sequences_sorted = reorder_dict_by_keys(pangenome_clusters_First_sequences, sorted_first_keys)
    pangenome_clusters_Type_sorted = reorder_dict_by_keys(pangenome_clusters_Type, sorted_first_keys)

    print("Calculating Groups")
    for cluster, numbers in pangenome_clusters_Type_sorted.items():
    ############################### Calculate First only
        calc_First_only_core(cluster, numbers[1],groups,cores)

    if options.reclustered != None:
        ############################# Calculate First and Reclustered-Second
        if numbers[0] == 1 and numbers[3] >= 1:  # If Seconds did not combine First reps
            calc_single_First_extended_Second_only_core(cluster, numbers[1], groups, cores, numbers[3])
        elif numbers[0] > 1 and numbers[3] >= 1:  # If unique Secondss combined multiple Firsts
            calc_multi_First_extended_Second_only_core(cluster, numbers[1], groups, cores, numbers[3])
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
            if data[1] >= 1:
                calc_Second_only_core(cluster, data[1], groups, cores)
        for cluster, data in combined_pangenome_clusters_ONLY_Second_Type.items():
            if data[1] >= 1:
                calc_only_Second_only_core(cluster, data[1], groups, cores)
    ###########################
    ### Output
    output_path = os.path.abspath(options.output_dir)
    stats_out = os.path.join(output_path,'summary_statistics.txt')
    key_order = ['First_core_', 'extended_core_', 'combined_core_', 'Second_core_','only_Second_core_']
    with open(stats_out, 'w') as outfile:
        print("Gene Groups:")
        outfile.write("Gene Groups:\n")
        for key_prefix in key_order:
            for key, value in cores.items():
                if key.startswith(key_prefix):
                    print(f"{key}: {len(value)}")
                    outfile.write(f"{key}: {len(value)}\n")
        print("Total Number of Gene Groups (Including Singletons): " + str(len(pangenome_clusters_First_sequences_sorted)))
        outfile.write("Total Number of Gene Groups (Including Singletons): " + str(len(pangenome_clusters_First_sequences_sorted)))

    if options.gene_presence_absence_out != None:
        gene_presence_absence_output(options,genome_dict, pangenome_clusters_First_sorted, pangenome_clusters_First_sequences_sorted)

    if options.write_families != None and options.fasta != None:
        sequences = read_fasta(options.fasta)
        output_dir = os.path.dirname(os.path.abspath(options.output_dir))
        output_dir = os.path.join(output_dir, 'Gene_Families_Output')

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
        process_gene_families(options, os.path.join(output_dir, 'Gene_Families_Output'), 'concatonated_genes_aligned.fasta')






