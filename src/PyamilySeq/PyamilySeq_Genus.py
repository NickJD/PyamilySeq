#from line_profiler_pycharm import profile


try:
    from .constants import *
    from .clusterings import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *
    from clusterings import *
    from utils import *


def gene_presence_absence_output(options, genus_dict, pangenome_clusters_First_sorted, pangenome_clusters_First_sequences_sorted):
    print("Outputting gene_presence_absence file")
    output_dir = os.path.abspath(options.output_dir)
    #in_name = options.clusters.split('.')[0].split('/')[-1]
    gpa_outfile = os.path.join(output_dir, 'gene_presence_absence.csv')
    gpa_outfile = open(gpa_outfile, 'w')
    genus_dict = OrderedDict(sorted(genus_dict.items()))
    gpa_outfile.write('"Gene","Non-unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Genome Fragment","Order within Fragment",'
                     '"Accessory Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc","')
    gpa_outfile.write('","'.join(genus_dict.keys()))
    gpa_outfile.write('"\n')
    for cluster, sequences in pangenome_clusters_First_sequences_sorted.items():
        average_sequences_per_genome = len(sequences) / len(pangenome_clusters_First_sorted[cluster])
        gpa_outfile.write('"group_'+str(cluster)+'","","","'+str(len(pangenome_clusters_First_sorted[cluster]))+'","'+str(len(sequences))+'","'+str(average_sequences_per_genome)+
                         '","","","","","","","",""')


        for genus in genus_dict.keys():
            full_out = ''
            tmp_list = []
            for value in sequences:
                if value.split('_')[0] == genus:
                    tmp_list.append(value)
            if tmp_list:
                full_out += ',"'+'  '.join(tmp_list)+'"'
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




def get_cores(options):
    ##Calculate core groups
    groups = OrderedDict()
    cores = OrderedDict()
    for group in options.core_groups.split(','):
        first_core_group = 'First_genera_' + group
        cores[first_core_group] = []
        if options.reclustered != None:
            extended_core_group = 'extended_genera_' + group
            cores[extended_core_group] = []
            combined_core_group = 'combined_genera_' + group
            cores[combined_core_group] = []
            second_core_group = 'Second_genera_' + group
            cores[second_core_group] = []
            only_second_core_group = 'only_Second_genera_' + group
            cores[only_second_core_group] = []
    return cores, groups

#@profile
def calc_First_only_core(cluster, First_num,  cores):
    try:
        cores['First_genera_' + str(First_num)].append(cluster)
    except KeyError:
        cores['First_genera_>'].append(cluster)
#@profile
def calc_single_First_extended_Second_only_core(cluster, First_num, cores, Second_num): # Count gene families extended with StORFs
    group = First_num + Second_num
    try:
        cores['extended_genera_' + str(group)].append(cluster)
    except KeyError:
        cores['extended_genera_>'].append(cluster)
#@profile
def calc_multi_First_extended_Second_only_core(cluster, First_num, cores, Second_num): # Count seperately those gene families extended with StORF_Reporter but combined >1 PEP
    group = First_num + Second_num
    try:
        cores['combined_genera_' + str(group)].append(cluster)
    except KeyError:
        cores['combined_genera_>'].append(cluster)
#@profile
def calc_Second_only_core(cluster, cores, Second_num):
    try:
        cores['Second_genera_' + str(Second_num)].append(cluster)
    except KeyError:
        cores['Second_genera_>'].append(cluster)
#@profile
def calc_only_Second_only_core(cluster, cores, Second_num): # only count the true storf onlies
    try:
        cores['only_Second_genera_' + str(Second_num)].append(cluster)
    except:
        cores['only_Second_genera_>'].append(cluster)



#@profile
def cluster(options):

    if options.cluster_format == 'CD-HIT':
        genus_dict, pangenome_clusters_First, pangenome_clusters_First_genomes, pangenome_clusters_First_sequences, reps = cluster_CDHIT(options, '_')
    elif 'TSV' in options.cluster_format or 'CSV' in options.cluster_format:
        genus_dict, pangenome_clusters_First, pangenome_clusters_First_genomes, pangenome_clusters_First_sequences, reps = cluster_EdgeList(options, '_')

    ###
    cores, groups = get_cores(options)
    ###

    if options.reclustered != None:
        if options.cluster_format == 'CD-HIT':
            combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids,combined_pangenome_clusters_Second, combined_pangenome_clusters_Second_sequences = combined_clustering_CDHIT(options, genus_dict, '_')
        elif 'TSV' in options.cluster_format or 'CSV' in options.cluster_format:
            combined_pangenome_clusters_First_Second_clustered,not_Second_only_cluster_ids,combined_pangenome_clusters_Second, combined_pangenome_clusters_Second_sequences = combined_clustering_Edge_List(options, '_')
        pangenome_clusters_Type = combined_clustering_counting(options, pangenome_clusters_First, reps, combined_pangenome_clusters_First_Second_clustered, pangenome_clusters_First_genomes, '_')
    else:
        pangenome_clusters_Type = single_clustering_counting(pangenome_clusters_First, reps)



    Number_Of_Second_Extending_But_Same_Genomes = 0

    sorted_first_keys = sort_keys_by_values(pangenome_clusters_First, pangenome_clusters_First_sequences)
    pangenome_clusters_First_sorted = reorder_dict_by_keys(pangenome_clusters_First, sorted_first_keys)
    pangenome_clusters_First_sequences_sorted = reorder_dict_by_keys(pangenome_clusters_First_sequences, sorted_first_keys)
    pangenome_clusters_Type_sorted = reorder_dict_by_keys(pangenome_clusters_Type, sorted_first_keys)

    print("Calculating Groups")
    seen_groupings = []
    for cluster, numbers in pangenome_clusters_Type_sorted.items():
    ############################### Calculate First only
        cluster = str(cluster)
        for grouping in numbers[2]: #!!# Could do with a more elegant solution
            current_cluster = grouping[0].split(':')[0]
            if current_cluster not in seen_groupings:
                seen_groupings.append(current_cluster)
                current_cluster_size = grouping[0].split(':')[1]
                calc_First_only_core(current_cluster, current_cluster_size, cores)
            ############################# Calculate First and Reclustered-Second
                if numbers[0] == 1 and numbers[3] >= 1:  # If Seconds did not combine First reps
                    calc_single_First_extended_Second_only_core(cluster, numbers[1], cores, numbers[3])
                elif numbers[0] > 1 and numbers[3] >= 1:  # If unique Seconds combined multiple Firsts
                    calc_multi_First_extended_Second_only_core(cluster, numbers[1], cores, numbers[3])
                elif numbers[4] >= 1:
                    Number_Of_Second_Extending_But_Same_Genomes += 1
            else:
                if options.verbose == True:
                    print("First cluster " + current_cluster + " already processed - This is likely because it was clustered with another First representative.")

    if options.reclustered != None:
        combined_pangenome_clusters_ONLY_Second_Type = defaultdict(list)
        combined_pangenome_clusters_Second_Type = defaultdict(list)
        for cluster, genomes in combined_pangenome_clusters_Second.items():
            if cluster in not_Second_only_cluster_ids:
                combined_pangenome_clusters_Second_Type[cluster] = [cluster, len(genomes)]
            else:
                combined_pangenome_clusters_ONLY_Second_Type[cluster] = [cluster, len(genomes)]
        for cluster, data in combined_pangenome_clusters_Second_Type.items():
            if data[1] >= 1:
                calc_Second_only_core(cluster, cores, data[1])
        for cluster, data in combined_pangenome_clusters_ONLY_Second_Type.items():
            if data[1] >= 1:
                calc_only_Second_only_core(cluster,  cores, data[1])
    ###########################
    ### Output
    output_path = os.path.abspath(options.output_dir)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    stats_out = os.path.join(output_path,'summary_statistics.txt')
    key_order = list(cores.keys())
    with open(stats_out,'w') as outfile:
        print("Genus Groups:")
        outfile.write("Genus Groups\n")
        for key in key_order:
            print(key+':\t'+str(len(cores[key])))
            outfile.write(key + ':\t' + str(len(cores[key]))+'\n')
        print("Total Number of First Gene Groups (Including Singletons): " + str(len(pangenome_clusters_First_sequences_sorted)))
        outfile.write("Total Number of Gene Groups (Including Singletons): " + str(len(pangenome_clusters_First_sequences_sorted)))
        if options.reclustered!= None:
            print("Total Number of Second Gene Groups (Including Singletons): " + str(
                len(combined_pangenome_clusters_Second_sequences)))
            print("Total Number of First Gene Groups That Had Additional Second Sequences But Not New Genomes: " + str(
                Number_Of_Second_Extending_But_Same_Genomes))
            outfile.write("\nTotal Number of Second Gene Groups (Including Singletons): " + str(
                len(combined_pangenome_clusters_Second_sequences)))
            outfile.write("\nTotal Number of First Gene Groups That Had Additional Second Sequences But Not New Genomes: " + str(
                Number_Of_Second_Extending_But_Same_Genomes))

    if options.gene_presence_absence_out != False:
        gene_presence_absence_output(options,genus_dict, pangenome_clusters_First_sorted, pangenome_clusters_First_sequences_sorted)

    if options.run_mode == 'Full':
        if options.reclustered == None:
            combined_pangenome_clusters_Second_sequences = None
        if options.write_groups != None:
            print("Outputting gene group FASTA files")
            sequences = read_fasta(options.fasta)
            #output_dir = os.path.dirname(os.path.abspath(options.output_dir))
            output_dir = os.path.join(options.output_dir, 'Gene_Groups_Output')
            write_groups_func(options,output_dir, key_order, cores, sequences,
                         pangenome_clusters_First_sequences_sorted, combined_pangenome_clusters_Second_sequences)

    elif options.run_mode == 'Partial':
        if options.reclustered == None:
            combined_pangenome_clusters_Second_sequences = None
        if options.write_groups != None and options.fasta != None:
            print("Outputting gene group FASTA files")
            sequences = read_fasta(options.fasta)
            #output_dir = os.path.dirname(os.path.abspath(options.output_dir))
            output_dir = os.path.join(options.output_dir, 'Gene_Groups_Output')
            write_groups_func(options,output_dir, key_order, cores, sequences,
                         pangenome_clusters_First_sequences_sorted, combined_pangenome_clusters_Second_sequences)


    # if options.write_groups != None and options.fasta != None:
    #     sequences = read_fasta(options.fasta)
    #     output_dir = os.path.join(output_path, 'Gene_Families_Output')
    #
    #     write_groups(options,output_dir, key_order, cores, sequences,
    #                  pangenome_clusters_First_sequences_sorted, combined_pangenome_clusters_Second_sequences)


    #!!# - Currently only align in Species Mode
    #if options.align_core != None and options.fasta != None and options.write_groups != None:
    #    process_gene_families(options, os.path.join(output_path, 'Gene_Families_Output'), 'concatenated_genes_aligned.fasta')




