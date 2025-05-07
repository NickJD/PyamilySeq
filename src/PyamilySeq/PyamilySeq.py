import argparse
#from config import config_params

try:
    from .PyamilySeq_Species import cluster as species_cluster
    from .PyamilySeq_Genus import cluster as genus_cluster
    from .constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from PyamilySeq_Species import cluster as species_cluster
    from PyamilySeq_Genus import cluster as genus_cluster
    from constants import *
    from utils import *




def run_cd_hit(options, input_file, clustering_output, clustering_mode):
    cdhit_command = [
        clustering_mode,
        '-i', input_file,
        '-o', clustering_output,
        '-c', f"{float(options.pident):.2f}",
        '-s', f"{float(options.len_diff):.2f}",
        '-T', str(options.threads),
        '-M', str(options.mem),
        '-d', "0",
        '-g', str(options.fast_mode),
        '-sc', "1",
        '-sf', "1"
    ]
    if options.verbose == True:
        subprocess.run(cdhit_command)
    else:
        subprocess.run(cdhit_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main():
    parser = argparse.ArgumentParser(description=f"PyamilySeq {PyamilySeq_Version}: A tool for gene clustering and analysis.")

    # Add subparsers for Full and Partial modes
    subparsers = parser.add_subparsers(dest="run_mode", required=True, help="Choose a mode: 'Full' or 'Partial'.")

    # Full Mode Subparser
    full_parser = subparsers.add_parser("Full",
                                        help="Full mode: PyamilySeq to cluster with CD-HIT and process output.")
    #full_parser.add_argument("-clustering_format", choices=['CD-HIT', 'MMseqs', 'BLAST'], required=True,
    #                         help="Clustering format to use: CD-HIT, MMseqs2, or BLAST.")
    full_parser.add_argument("-output_dir", required=True,
                             help="Directory for all output files.")
    full_parser.add_argument("-input_type", choices=['separate', 'combined', 'fasta'], required=True,
                             help="Type of input files: 'separate' for matching FASTA and GFF files, 'combined' for GFF+FASTA, or 'fasta' for a prepared FASTA file.")
    full_parser.add_argument("-input_dir", required=False,
                             help="Directory containing GFF/FASTA files - Use with -input_type separate/combined.")
    full_parser.add_argument("-input_fasta", required=False,
                             help="Input FASTA file - Use with - input_type fasta.")
    full_parser.add_argument("-name_split_gff", required=False,
                             help="Substring to split filenames and extract genome names for gff files (e.g., '_combined.gff3') - Use with -input_type separate/combined.")
    full_parser.add_argument("-name_split_fasta", required=False,
                             help="Substring to split filenames and extract genome names for fasta files if named differently to paired gff files (e.g., '_dna.fasta') - Use with -input_type separate/combined.")
    full_parser.add_argument("-sequence_type", choices=['AA', 'DNA'], default="AA", required=False,
                             help="Clustering mode: 'DNA' or 'AA'.")
    full_parser.add_argument("-gene_ident", default="CDS", required=False,
                             help="Gene identifiers to extract sequences (e.g., 'CDS, tRNA').")
    full_parser.add_argument("-c", type=str, dest="pident", default="0.90", required=False,
                             help="Sequence identity threshold for clustering (default: 0.90) - CD-HIT parameter '-c'.")
    full_parser.add_argument("-s", type=str, dest="len_diff", default="0.80", required=False,
                             help="Length difference threshold for clustering (default: 0.80) - CD-HIT parameter '-s'.")

    full_parser.add_argument("-fast_mode", action="store_true",
                             help="Enable fast mode for CD-HIT (not recommended) - CD-HIT parameter '-g'.")


    # Partial Mode Subparser
    partial_parser = subparsers.add_parser("Partial", help="Partial mode: PyamilySeq to process pre-clustered data.")
    partial_parser.add_argument("-clustering_format", choices=['CD-HIT', 'MMseqs', 'BLAST'], required=True,
                                help="Clustering format used: CD-HIT, MMseqs2, or BLAST.")
    partial_parser.add_argument("-cluster_file", required=True,
                                help="Cluster file containing pre-clustered groups from CD-HIT, MMseqs, BLAST etc.")
    partial_parser.add_argument("-original_fasta", required=True,
                                help="FASTA file used in pre-clustering (Provide sequences in DNA form).")
    partial_parser.add_argument("-output_dir", required=True,
                                help="Directory for all output files.")
    partial_parser.add_argument("-reclustered", required=False,
                                help="Clustering output file from a second round of clustering.")
    partial_parser.add_argument("-seq_tag", default="StORF", dest="sequence_tag", required=False,
                                help="Tag for distinguishing reclustered sequences.")

    # Common Grouping Arguments
    for subparser in [full_parser, partial_parser]:
        subparser.add_argument("-group_mode", choices=['Species', 'Genus'], default="Species", required=False,
                               help="Grouping mode: 'Species' or 'Genus'.")
        subparser.add_argument("-species_groups", default="99,95,15", required=False,
                               help="Gene groupings for 'Species' mode (default: '99,95,15').")
        subparser.add_argument("-genus_groups", default="1,2,3,4,5,6,7,8,9,10", required=False,
                               help="Gene groupings for 'Genus' mode (default: '1-10').")
        subparser.add_argument("-write_groups", default=None, dest="write_groups", required=False,
                               help="Output gene groups as a single FASTA file (e.g., '99,95'). Triggers writing individual groups.")
        subparser.add_argument("-write_individual_groups", action="store_true", dest="write_individual_groups",
                               help="Output individual FASTA files for each group.")
        subparser.add_argument("-align", action="store_true", dest="align_core",
                               help="Align and concatenate sequences for 'core' groups (those in 99-100%% of genomes).")
        subparser.add_argument("-align_aa", action="store_true",
                               help="Align sequences as amino acids.")
        subparser.add_argument("-no_gpa", action="store_false", dest="gene_presence_absence_out",
                               help="Skip creation of gene_presence_absence.csv.")
        subparser.add_argument("-M", type=int, default=4000, dest="mem", required=False,
                                 help="Memory allocation for clustering (MB) - CD-HIT parameter '-M'.")
        subparser.add_argument("-T", type=int, default=8, dest="threads", required=False,
                                 help="Number of threads for clustering/alignment - CD-HIT parameter '-T' | MAFFT parameter '--thread'.")

        # Miscellaneous Arguments
        subparser.add_argument("-verbose", action="store_true",
                            help="Print verbose output.")
        subparser.add_argument("-v", "--version", action="version",
                            version=f"PyamilySeq {PyamilySeq_Version}: Exiting.")

    # Parse Arguments
    options = parser.parse_args()
    ## Configuration


    if options.write_groups != None and options.write_individual_groups == False:
        options.write_individual_groups = True

    # Example of conditional logic based on selected mode
    print(f"Running PyamilySeq {PyamilySeq_Version} in {options.run_mode} mode:")
    if options.run_mode == "Full" and options.verbose == True:
        print("Processing Full mode with options:", vars(options))
    elif options.run_mode == "Partial" and options.verbose == True:
        print("Processing Partial mode with options:", vars(options))

    ### Checking all required parameters are provided by user #!!# Doesn't seem to work
    if options.run_mode == 'Full':
        options.clustering_format = 'CD-HIT'
        if getattr(options, 'reclustered', None) is not None:
            sys.exit("Currently reclustering only works on Partial Mode.")
        required_full_mode = [options.input_type, options.pident, options.len_diff]
        if options.input_type != 'fasta':
            required_full_mode.extend([options.input_dir, options.name_split_gff])
        if all(required_full_mode):
            # Proceed with the Full mode
            pass
        else:
            missing_options = [opt for opt in
                               ['input_type', 'input_dir', 'name_split_gff', 'clustering_format', 'pident', 'len_diff'] if
                               not options.__dict__.get(opt)]
            sys.exit(f"Missing required options for Full mode: {', '.join(missing_options)}")
        if options.align_core:
            options.write_individual_groups = True
            if options.write_groups == None:
                sys.exit('Must provide "-write_groups" to output gene groups before alignment "-align" can be done.')
    elif options.run_mode == 'Partial':
        required_partial_mode = [options.cluster_file, options.original_fasta]
        if all(required_partial_mode):
            # Proceed with the Partial mode
            pass
        else:
            missing_options = [opt for opt in
                               ['cluster_file','original_fasta'] if
                               not options.__dict__[opt]]
            sys.exit(f"Missing required options for Partial mode: {', '.join(missing_options)}")
        if options.align_core:
            options.write_individual_groups = True
            if options.write_groups == None or options.original_fasta == None:
                sys.exit('Must provide "-w" to output gene groups before alignment "-a" can be done.')

    if options.clustering_format == 'CD-HIT':
        clust_affix = '.clstr'
    elif options.clustering_format == 'TSV':
        clust_affix = '.tsv'
    elif options.clustering_format == 'CSV':
        clust_affix = '.csv'

    ###External tool checks:
    ##MAFFT
    if options.align_core == True:
        if is_tool_installed('mafft'):
            if options.verbose == True:
                print("mafft is installed. Proceeding with alignment.")
        else:
            exit("mafft is not installed. Please install mafft to proceed.")
    ##CD-HIT
    if options.run_mode == 'Full':
        if is_tool_installed('cd-hit'):
            if options.verbose == True:
                print("cd-hit is installed. Proceeding with clustering.")
            if options.sequence_type == 'DNA':
                clustering_mode = 'cd-hit-est'
            elif options.sequence_type == 'AA':
                clustering_mode = 'cd-hit'
            if options.fast_mode == True:
                options.fast_mode = 1
                if options.verbose == True:
                    print("Running CD-HIT in fast mode.")
            else:
                options.fast_mode = 0
                if options.verbose == True:
                    print("Running CD-HIT in accurate mode.")
        else:
            exit("cd-hit is not installed. Please install cd-hit to proceed.")


    # if options.write_groups != None and options.original_fasta == False:
    #     exit("-fasta must br provided if -w is used")

    if hasattr(options, 'cluster_file') and options.cluster_file:
        options.cluster_file = fix_path(options.cluster_file)
    if hasattr(options, 'reclustered') and options.reclustered:
        options.reclustered = fix_path(options.reclustered)
    if hasattr(options, 'input_dir') and options.input_dir:
        options.input_dir = fix_path(options.input_dir)
    if hasattr(options, 'output_dir') and options.output_dir:
        options.output_dir = fix_path(options.output_dir)

    output_path = os.path.abspath(options.output_dir)
    combined_out_file = os.path.join(output_path, "combined_sequences_dna.fasta")
    clustering_output = os.path.join(output_path, 'clustering_' + options.clustering_format)

    if options.group_mode == 'Species':
        options.species_groups = options.species_groups + ',0'
        groups_to_use = options.species_groups
    elif options.group_mode == 'Genus':
        options.genus_groups = options.genus_groups + ',>'
        groups_to_use = options.genus_groups
        if options.align_core != None:
            sys.exit("-a align_core not a valid option in Genus mode.")


    if options.run_mode == 'Full':
        if options.clustering_format != 'CD-HIT':
            sys.exit('Only CD-HIT clsutering works in Full Mode')

        if not os.path.exists(output_path):
            os.makedirs(output_path)
        if options.sequence_type == 'AA':
            clustering_mode = 'cd-hit'
            file_to_cluster = combined_out_file.replace('_dna.fasta','_aa.fasta')
            translate = True
        elif options.sequence_type == 'DNA':
            clustering_mode = 'cd-hit-est'
            translate = False
            file_to_cluster = combined_out_file
        if options.input_type == 'separate':
            read_separate_files(options.input_dir, options.name_split_gff, options.name_split_fasta, options.gene_ident, combined_out_file, translate, False)
            run_cd_hit(options, file_to_cluster, clustering_output, clustering_mode)
        elif options.input_type == 'combined':
            read_combined_files(options.input_dir, options.name_split_gff, options.gene_ident, combined_out_file, translate, False)
            run_cd_hit(options, file_to_cluster, clustering_output, clustering_mode)
        elif options.input_type == 'fasta':
            combined_out_file = options.input_fasta
            ### FIX write code to detect if DNA or AA and if sequence tpye is AA then translate
            # Detect if the input FASTA file contains DNA or AA sequences
            is_dna = detect_sequence_type(options.input_fasta)
            # If the sequence type is AA and the input is DNA, translate the DNA to AA
            if options.sequence_type == 'AA' and is_dna:
                translated_fasta = os.path.join(output_path, os.path.splitext(os.path.basename(options.input_fasta))[0] + '_aa.fasta')
                translate_dna_to_aa(options.input_fasta, translated_fasta)
                file_to_cluster = translated_fasta
            else:
                file_to_cluster = options.input_fasta
            run_cd_hit(options, file_to_cluster, clustering_output, clustering_mode)


        class clustering_options:
            def __init__(self):
                self.run_mode = options.run_mode
                self.cluster_format = options.clustering_format
                self.sequence_type = options.sequence_type
                self.reclustered = None
                self.sequence_tag = None
                self.species_groups = groups_to_use
                self.clusters = clustering_output + clust_affix
                self.output_dir = options.output_dir
                self.gene_presence_absence_out = options.gene_presence_absence_out
                self.write_groups = options.write_groups
                self.write_individual_groups = options.write_individual_groups
                self.threads = options.threads
                self.align_core = options.align_core
                self.align_aa = options.align_aa
                self.fasta = combined_out_file
                self.verbose = options.verbose

        clustering_options = clustering_options()

    elif options.run_mode == 'Partial':
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        class clustering_options:
            def __init__(self):
                self.run_mode = options.run_mode
                self.cluster_format = options.clustering_format
                self.sequence_type = None
                self.reclustered = options.reclustered
                self.sequence_tag = options.sequence_tag
                self.species_groups = groups_to_use
                self.clusters = options.cluster_file
                self.output_dir = options.output_dir
                self.gene_presence_absence_out = options.gene_presence_absence_out
                self.write_groups = options.write_groups
                self.write_individual_groups = options.write_individual_groups
                self.threads = options.threads
                self.align_core = options.align_core
                self.align_aa = options.align_aa
                self.fasta = options.original_fasta
                self.verbose = options.verbose

        clustering_options = clustering_options()


    if options.group_mode == 'Species':
        species_cluster(clustering_options)
    elif options.group_mode == 'Genus':
        genus_cluster((clustering_options))


    # Save arguments to a text file
    with open(output_path+"/PyamilySeq_params.txt", "w") as outfile:
        for arg, value in vars(options).items():
            outfile.write(f"{arg}: {value}\n")

    print("Thank you for using PyamilySeq -- A detailed user manual can be found at https://github.com/NickJD/PyamilySeq\n"
          "Please report any issues to: https://github.com/NickJD/PyamilySeq/issues\n#####")

if __name__ == "__main__":
    main()
