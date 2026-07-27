from collections import defaultdict, OrderedDict


try:
    from .constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *
    from utils import *

def run_cd_hit(options, input_file, clustering_output, clustering_mode):
    logger = logging.getLogger("PyamilySeq.Group_Splitter")
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
    logger.debug("Group-Splitter CD-HIT command: %s", " ".join(cdhit_command))
    try:
        if options.verbose:
            result = subprocess.run(cdhit_command)
        else:
            result = subprocess.run(cdhit_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if result.returncode != 0:
            raise RuntimeError(f"{clustering_mode} returned exit code {result.returncode}")
        logger.info("CD-HIT completed for %s", input_file)
    except Exception:
        logger.exception("Error running CD-HIT for %s", input_file)
        raise


def calculate_new_rep_seq(cluster_data, length_weight=1.0, identity_weight=1.0):
    total_length = sum(entry['length'] for entry in cluster_data)
    avg_length = total_length / len(cluster_data)

    total_identity = sum(entry['percent_identity'] for entry in cluster_data)
    avg_identity = total_identity / len(cluster_data)

    # Normalise length and identity
    max_length = max(entry['length'] for entry in cluster_data)
    max_identity = 100  # Assuming percent_identity is out of 100

    # Calculate a score based on both length difference and percent identity
    def score(entry):
        normalized_length_diff = abs(entry['length'] - avg_length) / max_length
        normalized_identity_diff = abs(entry['percent_identity'] - avg_identity) / max_identity
        return (length_weight * normalized_length_diff) + (identity_weight * normalized_identity_diff)

    rep_entry = min(cluster_data, key=score)
    return rep_entry


def parse_group_sequence_header(header):
    """Return family, genome, and gene tokens from legacy or canonical headers."""

    clean = header.strip().lstrip('>')
    parts = clean.split('|', 2)
    if len(parts) < 2:
        raise ValueError(f"Group FASTA header lacks a genome-qualified member ID: {header!r}")
    family = parts[0]
    genome = parts[1]
    gene = parts[2] if len(parts) == 3 else parts[1]
    return family, genome, gene



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
    # Ensure selected_group_ids is always defined to satisfy static analysis
    selected_group_ids = None

    if options.sequence_type == 'AA':
        affix = '_aa.fasta'
    else:
        affix = '_dna.fasta'

    # Ensure we look for the combined file that includes the requested group level (e.g. "99")
    # groups_to_use[1] contains the numeric group level when using ('groups', <num>)
    if groups_to_use and groups_to_use[0] == 'ids':
        group_level = str(getattr(options, 'groups', None) or 99)
    else:
        group_level = str(groups_to_use[1]) if groups_to_use and len(groups_to_use) > 1 else ''
    # combined groups are read from the original PyamilySeq Gene_Groups_Output directory
    gene_groups_input_dir = getattr(options, 'gene_groups_input', os.path.join(options.input_directory, 'Gene_Groups_Output'))
    combined_groups_fasta = os.path.join(gene_groups_input_dir, f"combined_group_sequences_{group_level}{affix}")

    # Defensive check: combined_group_sequences_* file must exist (was created by PyamilySeq with -write_groups)
    if not os.path.exists(combined_groups_fasta):
        logger = logging.getLogger("PyamilySeq.Group_Splitter")
        logger.error("Required combined group sequences file not found: %s", combined_groups_fasta)
        logger.error("This usually means the upstream PyamilySeq run did not include the -write_groups and -write_individual_groups options.")
        # Helpful debug info: list contents of Gene_Groups_Output if available
        parent_dir = os.path.dirname(combined_groups_fasta)
        if os.path.isdir(parent_dir):
            try:
                files = os.listdir(parent_dir)
                logger.debug("Files in %s: %s", parent_dir, ", ".join(sorted(files)) if files else "(none)")
            except Exception as e:
                logger.debug("Could not list %s: %s", parent_dir, e)
        # Stop further processing
        raise SystemExit(1)

    if groups_to_use[0] == 'ids':
        selected_group_ids = {g.strip().lstrip('>') for g in groups_to_use[1].split(',') if g.strip()}
    elif groups_to_use[0] == 'groups':
        # The combined file was already filtered by PyamilySeq's -write_groups
        # setting, so every family in it is eligible.
        selected_group_ids = None


    group_number = None
    with open(combined_groups_fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_group is not None and (selected_group_ids is None or group_number in selected_group_ids):
                    groups[current_group].append((current_group_header, ''.join(current_sequence)))

                current_group_header = line.strip()
                family, genome, _ = parse_group_sequence_header(current_group_header)
                current_group = '>' + family
                current_sequence = []
                genome_count[genome] += 1

                # Only process if group matches the selected_groups or if no specific groups were provided
                group_number = family
                legacy_number = family.replace('Group_', '').replace('group_', '')
                if selected_group_ids is not None and family not in selected_group_ids and legacy_number not in selected_group_ids:
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
    # Report list - each entry will be a dict with information about the split of a group
    split_report = []

    # Ensure variable exists for later checks
    all_same = False

    for group_header, sequences in groups.items():
        if options.verbose == True:
            print(f"\n###\nCurrent Group: {group_header.replace('>','')}\n")

        group_name = group_header.split('|')[0]  # Get the group part (e.g., '>Group_n')

        # Count genomes with more than one gene
        genome_to_gene_count = defaultdict(int)
        for header, _ in sequences:
            _, genome, _ = parse_group_sequence_header(header)
            genome_to_gene_count[genome] += 1

        num_genomes_with_multiple_genes = sum(1 for count in genome_to_gene_count.values() if count > 1)

        # Prepare a report entry for this group
        report_entry = {
            'group_header': group_header,
            'group_name': group_name.replace('>',''),
            'total_sequences': len(sequences),
            'num_genomes_with_multiple_genes': num_genomes_with_multiple_genes,
            'percent_genomes_with_multiple_genes': (num_genomes_with_multiple_genes / options.genome_num * 100) if options.genome_num else None,
            'split': False,
            'reason': None,
            'parameters': {
                'pident': options.pident,
                'len_diff': options.len_diff,
                'fast_mode': options.fast_mode
            },
            'num_subgroups': 0,
            'subgroup_sizes': [],
            'subgroup_representatives': []
        }

        # Check if the group meets the threshold for having paralogs
        #if options.groups == None:
        if (num_genomes_with_multiple_genes / options.genome_num) * 100 < options.group_threshold:
            # If not enough genomes have multiple genes this group is skipped, record reason and continue
            report_entry['reason'] = 'below_group_threshold'
            split_report.append(report_entry)
            continue


        group_file_name = group_name.replace('>','')

        temp_fasta = f"{options.gene_groups_output}/{group_file_name}.fasta"
        write_fasta(sequences, temp_fasta)

        # Run cd-hit on the individual group
        clustering_output = f"{options.gene_groups_output}/{group_file_name}_clustering"

        run_cd_hit(options, temp_fasta, clustering_output, clustering_mode)

        # Read the clustering results to find subgroups
        clustered_sequences = read_cd_hit_output(clustering_output + '.clstr')

        # If only a single cluster present we may still need to check if all sequences are identical
        if len(clustered_sequences) == 1:
            # Detect if all sequences are identical in length and percentage identity
            all_same = check_if_all_identical(clustered_sequences)
        else:
            all_same = False

        # **Global subgroup counter for the entire major group**
        subgroup_id = 0


        if not all_same:
            # Iterate through each cluster in clustered_sequences
            for cluster_key, cluster in clustered_sequences.items():

                remaining_sequences_tmp = sequences.copy()  # Track unprocessed sequences
                remaining_sequences = [entry for entry in remaining_sequences_tmp if entry[0] in
                                       {seq_entry['header'] for seq_entry in cluster}]
                while remaining_sequences:
                    # Track subgroups for this cluster pass
                    subgroup_sequences = []
                    sequences_to_remove = []

                    # Recalculate representative sequence dynamically for this cluster
                    rep = calculate_new_rep_seq(
                        [entry for entry in cluster if entry['header'] in (h for h, _ in remaining_sequences)]
                    )

                    # Find the sequence corresponding to rep['header'] from the list of sequences
                    rep_seq = next((seq for header, seq in sequences if header == rep['header']), None)

                    # Save previously checked seqs, so we don't have to compare them again.
                    # use local defaultdict (the module 'collections' may not be imported as a name in this file)
                    checked = defaultdict(float)

                    # Process each genome to select the best matching sequence
                    for genome in genome_to_gene_count:
                        best_sequence = None
                        best_score = None  # Initialise with a very low score, so that even negative scores can be selected

                        # Iterate over each sequence in the remaining sequences for this genome
                        for header, seq in remaining_sequences:
                            _, genome_id, _ = parse_group_sequence_header(header)

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
                            _, best_genome, best_gene = parse_group_sequence_header(best_sequence[0])
                            new_header = f">{group_file_name}_subgroup_{subgroup_id}|{best_genome}|{best_gene}"
                            subgroup_sequences.append((new_header, best_sequence[1]))
                            sequences_to_remove.append(best_sequence)

                    # Write each subgroup into a separate FASTA file
                    if subgroup_sequences:
                        sequence_suffix = '_aa.fasta' if options.sequence_type == 'AA' else '_dna.fasta'
                        subgroup_file = f"{options.sub_groups_output}/{group_file_name}_subgroup_{subgroup_id}{sequence_suffix}"
                        write_fasta(subgroup_sequences, subgroup_file)

                        # Record subgroup info for report
                        report_entry['split'] = True
                        report_entry['reason'] = (
                            'cdhit_multiple_clusters'
                            if len(clustered_sequences) > 1
                            else 'nonidentical_paralog_partition'
                        )
                        report_entry['num_subgroups'] += 1
                        report_entry['subgroup_sizes'].append(len(subgroup_sequences))
                        # representative header for this subgroup is the first sequence's header (best sequence for a genome)
                        try:
                            rep_header = subgroup_sequences[0][0]
                        except Exception:
                            rep_header = None
                        report_entry['subgroup_representatives'].append(rep_header)

                        # Remove processed sequences from the remaining list
                        remaining_sequences = [item for item in remaining_sequences if
                                               item[0] not in {h for h, _ in sequences_to_remove}]

                        # Increment subgroup ID for the next subgroup
                        subgroup_id += 1
                        paralog_groups[group_name]['count'] += 1  # Count this group as a paralog group
                        paralog_groups[group_name]['sizes'].append(len(subgroup_sequences))  # Record the size of the subgroup



        else:
            # Condition 2: If sequences are identical, distribute genes evenly into subgroups
            num_subgroups = max(genome_to_gene_count.values(), default=1)
            subgroup_sequences = defaultdict(list)  # Store sequences for each subgroup
            genome_count = defaultdict(int)  # Count how many genes have been assigned to each genome

            # Iterate over all sequences regardless of whether the genome has been seen
            for header, seq in sequences:
                _, genome, gene = parse_group_sequence_header(header)

                # Determine the next subgroup for this genome
                subgroup_id = genome_count[genome] % num_subgroups
                new_header = f">{group_file_name}_subgroup_{subgroup_id}|{genome}|{gene}"
                subgroup_sequences[subgroup_id].append((new_header, seq))

                 # Increment the count for this genome
                genome_count[genome] += 1

            # Write out each subgroup to a separate FASTA file
            for subgroup_id, seqs in subgroup_sequences.items():
                sequence_suffix = '_aa.fasta' if options.sequence_type == 'AA' else '_dna.fasta'
                subgroup_file = os.path.join(
                    options.sub_groups_output,
                    f"{group_file_name}_subgroup_{subgroup_id}{sequence_suffix}",
                )
                write_fasta(seqs, subgroup_file)

                # Increment subgroup ID globally for the next subgroup
                subgroup_id += 1
                paralog_groups[group_name]['count'] += 1  # Count this group as a paralog group
                paralog_groups[group_name]['sizes'].append(len(seqs))  # Record the size of the subgroup

            # Record report for identical case
            report_entry['split'] = True
            report_entry['reason'] = 'identical_sequences_distributed'
            report_entry['num_subgroups'] = len(subgroup_sequences)
            report_entry['subgroup_sizes'] = [len(v) for v in subgroup_sequences.values()]
            # Use first header of each subgroup as representative
            report_entry['subgroup_representatives'] = [v[0][0] if v else None for v in subgroup_sequences.values()]


        # Clean up temporary fasta file if the option is set
        if options.delete_temp_files:
            if temp_fasta and os.path.exists(temp_fasta):
                os.remove(temp_fasta)
            if os.path.exists(clustering_output + '.clstr'):
                os.remove(clustering_output + '.clstr')
            if os.path.exists(clustering_output):
                os.remove(clustering_output)

        # Append the report entry for this group
        split_report.append(report_entry)


    return paralog_groups, split_report


def write_group_report(report_entries, output_dir):
    """Write a TSV and JSON report summarizing splits for each processed group.

    report_entries: list of dicts
    output_dir: directory where report files will be written
    """
    import json
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    tsv_path = os.path.join(output_dir, 'Group-Splitter_report.tsv')

    with open(tsv_path, 'w') as tsv:
        # define columns for TSV
        cols = [
            'group_header', 'group_name', 'total_sequences', 'num_genomes_with_multiple_genes', 'percent_genomes_with_multiple_genes',
            'split', 'reason', 'num_subgroups', 'subgroup_sizes', 'subgroup_representatives', 'parameters'
        ]
        tsv.write('#' + '\t'.join(cols) + '\n')
        for entry in report_entries:
            row = []
            for c in cols:
                val = entry.get(c, '')
                if isinstance(val, list):
                    val = ';'.join(str(x) for x in val)
                elif isinstance(val, dict):
                    val = json.dumps(val)
                row.append(str(val) if val is not None else '')
            tsv.write('\t'.join(row) + '\n')


    logger = logging.getLogger('PyamilySeq.Group_Splitter')
    logger.info('Wrote Group-Splitter report TSV: %s', tsv_path)


def main():
    # Early console-only logger so parser.description is emitted via logger before argparse prints usage/help.
    early_logger = configure_logger("PyamilySeq.Group_Splitter", enable_file=False, log_dir=None, verbose=False)
    # Use LoggingArgumentParser so usage/errors are emitted via the configured logger
    parser = LoggingArgumentParser(logger_name="PyamilySeq.Group_Splitter", description='Group-Splitter - A tool to split multi-copy gene groups identified by PyamilySeq.')
    ### Required Arguments
    required = parser.add_argument_group('Required Parameters')
    required.add_argument('-input_dir', action='store', dest='input_directory',
                          help='Provide the directory of a PyamilySeq run.',
                          required=True)
    required.add_argument('-seq_type', action='store', dest='sequence_type', default='AA',choices=['AA', 'DNA'],
                          help='Default - AA: Are groups "DNA" or "AA" sequences?',
                          required=True)
    required.add_argument('-genome_num', action='store', dest='genome_num', type=int,
                          help='The total number of genomes must be provide',
                          required=True)


    ### Regrouping Arguments
    regrouping_params = parser.add_argument_group('Regrouping Parameters')
    regrouping_params.add_argument('-groups', action="store", dest='groups', type=int, default=99,
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


    # Optional file logging flags (must be added before parsing)
    parser.add_argument("--log", action="store_true", dest="log", help="Create a timestamped logfile for this run.")
    parser.add_argument("--log-dir", dest="log_dir", default=None, help="Directory for logfile (default: input_directory).")
    options = parser.parse_args()
    if options.genome_num <= 0:
        parser.error('-genome_num must be a positive integer')
    if not 0 <= options.group_threshold <= 100:
        parser.error('-group_threshold must be between 0 and 100')
    if not 0 < options.pident <= 1:
        parser.error('-c must be greater than 0 and at most 1')
    if not 0 <= options.len_diff <= 1:
        parser.error('-s must be between 0 and 1')
    options.align_aa = options.sequence_type == 'AA'
    # Compute logfile directory (default to input_directory) and only enable file logging when --log is provided.
    log_dir = options.log_dir if getattr(options, "log_dir", None) else os.path.abspath(options.input_directory)
    logger = configure_logger("PyamilySeq.Group_Splitter", enable_file=getattr(options, "log", False), log_dir=log_dir, verbose=options.verbose)
    logger.info("Running Group-Splitter %s", PyamilySeq_Version)

    ###External tool checks:
    ##MAFFT
    if options.align_core is not None:
        if is_tool_installed('mafft'):
            if options.verbose == True:
                print("mafft is installed. Proceeding with alignment.")
        else:
            raise SystemExit("mafft is not installed. Please install mafft to proceed.")

    ##CD-HIT
    clustering_mode = 'cd-hit-est' if options.sequence_type == 'DNA' else 'cd-hit'
    if is_tool_installed(clustering_mode):
        logger.info("%s is installed. Proceeding with clustering.", clustering_mode)
        if options.fast_mode == True:
            options.fast_mode = 0
            if options.verbose == True:
                print("Running CD-HIT in fast mode.")
        else:
            options.fast_mode = 1
            if options.verbose == True:
                print("Running CD-HIT in slow mode.")
    else:
        logger.error("%s is not installed. Please install CD-HIT to proceed.", clustering_mode)
        raise SystemExit(f"{clustering_mode} is not installed. Please install CD-HIT to proceed.")

    ##Alignment
    if options.align_core != None:
        try:
            options.align_core_threshold = int(str(options.align_core).split(',')[0])
        except ValueError:
            raise SystemExit('"-a" must start with an integer core threshold, e.g. "-a 99".')
        if options.groups == None and options.group_ids == None:
            raise SystemExit('Must provide "-groups" or "-group_ids" when requesting alignment with "-a".')

    ## Output directories: create a top-level Group_Splitting directory for all outputs of this tool
    base_output = os.path.join(options.input_directory, 'Group_Splitting')
    options.base_output = base_output
    # original gene groups (from PyamilySeq) live under input_directory/Gene_Groups_Output
    options.gene_groups_input = os.path.join(options.input_directory, 'Gene_Groups_Output')
    # this tool's gene groups output (temp clustering files) will be inside the new Group_Splitting dir
    options.gene_groups_output = os.path.join(base_output, 'Gene_Groups_Output')
    options.sub_groups_output = os.path.join(base_output, 'Sub_Groups_Output')
    # ensure directories exist
    if not os.path.exists(base_output):
        os.makedirs(base_output, exist_ok=True)
    if not os.path.exists(options.gene_groups_output):
        os.makedirs(options.gene_groups_output, exist_ok=True)
    if not os.path.exists(options.sub_groups_output):
        os.makedirs(options.sub_groups_output, exist_ok=True)

    logger.info("Group-Splitting base output: %s", base_output)
    logger.info("Gene groups (input): %s", options.gene_groups_input)
    logger.info("Gene groups (tool output): %s", options.gene_groups_output)
    logger.info("Sub groups output: %s", options.sub_groups_output)

    # Summary statistics file produced by upstream PyamilySeq run (kept in input directory)
    summary_file = os.path.join(options.input_directory, 'summary_statistics.txt')

    # Save arguments to a text file (inside Group_Splitting)
    params_out = os.path.join(options.base_output, 'Group-Splitter_params.txt')
    with open(params_out, "w") as outfile:
        for arg, value in vars(options).items():
            outfile.write(f"{arg}: {value}\n")
    logger.info("Saved parameters to %s", params_out)

    ## Group Selection
    if options.group_ids != None:
        groups_to_use = ('ids', options.group_ids)
    elif options.groups != None:
        groups_to_use = ('groups', options.groups)
    else:
        groups_to_use = ('groups', 99)



    paralog_groups, split_report = separate_groups(options, clustering_mode, groups_to_use)
    # Write a comprehensive TSV report for the splitting operation (no JSON)
    try:
        write_group_report(split_report, os.path.abspath(options.base_output))
    except Exception:
        logger.exception('Failed to write Group-Splitter report')
    logger.info("Identified %d paralog group(s)", len(paralog_groups))
    for group_id, data in paralog_groups.items():
        logger.debug("Group %s -> new groups: %s sizes: %s", group_id, data['count'], data['sizes'])


    # Recalculate every configured First_core_* bin rather than assuming 99/95/15/0.
    with open(summary_file, encoding='utf-8') as handle:
        summary_data = handle.read().splitlines()
    summary_info = {}
    for line in summary_data:
        if ':' not in line:
            continue
        key, value = line.split(':', 1)
        try:
            summary_info[key.strip()] = int(value.strip())
        except ValueError:
            continue
    genome_num = next(
        (value for key, value in summary_info.items() if key.lower() == 'number of genomes'),
        options.genome_num,
    )
    category_counts = OrderedDict(
        (key, value) for key, value in summary_info.items() if key.startswith('First_core_')
    )
    if not category_counts:
        raise RuntimeError(f'No First_core_* categories found in {summary_file}')
    percentage_bins = sorted({
        int(key.rsplit('_', 1)[1])
        for key in category_counts
        if key.rsplit('_', 1)[1].isdigit()
    }, reverse=True)

    def source_category_for(family_id):
        family_id = family_id.lstrip('>')
        legacy_id = family_id.replace('Group_', '').replace('group_', '')
        for filename in os.listdir(options.gene_groups_input):
            if not filename.endswith(('.fasta', '.fa')):
                continue
            for candidate in (family_id, legacy_id):
                marker = f'_{candidate}_'
                if marker in filename:
                    return filename.split(marker, 1)[0]
                if filename.endswith(f'_{candidate}.fasta'):
                    return filename[:-len(f'_{candidate}.fasta')]
        return None

    updated_counts = OrderedDict(category_counts)
    original_total = next(
        (
            value for key, value in summary_info.items()
            if key.lower() in {
                'final gene families',
                'total number of first gene groups (including singletons)',
            }
        ),
        sum(category_counts.values()),
    )
    total_gene_groups = original_total
    for family_id, data in paralog_groups.items():
        source_category = source_category_for(family_id)
        if source_category not in updated_counts:
            logger.warning(
                'Could not determine a First_core category for %s from %s',
                family_id,
                options.gene_groups_input,
            )
            continue
        updated_counts[source_category] -= 1
        total_gene_groups += max(0, len(data['sizes']) - 1)
        for size in data['sizes']:
            selected_bin = next(
                (
                    percentage for percentage in percentage_bins
                    if size >= (0 if percentage == 0 else math.ceil(percentage * genome_num / 100))
                ),
                percentage_bins[-1],
            )
            target = f'First_core_{selected_bin}'
            updated_counts.setdefault(target, 0)
            updated_counts[target] += 1

    stats_out = os.path.join(
        options.base_output,
        os.path.basename(summary_file).replace('.txt', '_recalculated.txt'),
    )
    with open(stats_out, 'w', encoding='utf-8') as outfile:
        outfile.write(f'Number of Genomes: {genome_num}\n')
        outfile.write('Recalculated Gene Groups\n')
        for category, count in updated_counts.items():
            outfile.write(f'{category}: {count}\n')
            logger.info('%s: %s', category, count)
        outfile.write(
            f'Total Number of First Gene Groups (Including Singletons): {total_gene_groups}\n'
        )

    # Alignment
    if options.align_core != None:
        logger.info("Processing gene group alignment")
        group_directory = options.gene_groups_input
        sub_group_directory = options.sub_groups_output
        combined_groups_fasta = os.path.join(options.gene_groups_input, f"combined_group_sequences_{options.groups}_dna.fasta")
        genome_list = read_genomes_from_fasta(combined_groups_fasta)
        process_gene_groups(options, group_directory, sub_group_directory, paralog_groups, genome_list, 'concatenated_genes_post_splitting_aligned_dna.fasta')







if __name__ == "__main__":

    main()
