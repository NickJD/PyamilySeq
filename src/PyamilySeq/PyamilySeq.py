import argparse
import collections
import os
import glob
import subprocess
from PyamilySeq_Species import *


try:
    from .PyamilySeq_Species import cluster
    from .Constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from PyamilySeq_Species import cluster
    from Constants import *
    from utils import *




def run_cd_hit(input_file, clustering_output, options):
    cdhit_command = [
        'cd-hit-est',
        '-i', input_file,
        '-o', clustering_output,
        '-c', str(options.pident),
        '-s', str(options.len_diff),
        '-T', "20",
        '-d', "0",
        '-sc', "1",
        '-sf', "1"
    ]
    if options.verbose == True:
        subprocess.run(cdhit_command)
    else:
        subprocess.run(cdhit_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main():
    parser = argparse.ArgumentParser(description='PyamilySeq ' + PyamilySeq_Version + ': PyamilySeq Run Parameters.')
    ### Required Arguments
    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-run_mode', action='store', dest='run_mode', choices=['Full','Partial'],
                          help='Run Mode: Should PyamilySeq be run in "Full" or "Partial" mode?',
                          required=True)
    required.add_argument('-group_mode', action='store', dest='group_type', choices=['Species','Genus'],
                          help='Group Mode: Should PyamilySeq be run in "Species" or "Genus" mode?',
                          required=True)
    required.add_argument("-clust_tool", action="store", dest="clust_tool", choices=['CD-HIT'],
                          help="Clustering tool to use: CD-HIT, DIAMOND, BLAST or MMseqs2.",
                          required=True)
    required.add_argument("-output_dir", action="store", dest="output_dir",
                          help="Directory for all output files.",
                          required=True)
    ### Full-Mode Arguments
    full_mode_args = parser.add_argument_group('Full-Mode Arguments - Required when "-run_mode Full" is used')
    full_mode_args.add_argument("-input_type", action="store", dest="input_type", choices=['separate', 'combined'],
                          help="Type of input files: 'separate' for separate FASTA and GFF files,"
                             " 'combined' for GFF files with embedded FASTA sequences.",
                          required=False)
    full_mode_args.add_argument("-input_dir", action="store", dest="input_dir",
                          help="Directory containing GFF/FASTA files.",
                          required=False)
    full_mode_args.add_argument("-name_split", action="store", dest="name_split",
                          help="substring used to split the filename and extract the genome name ('_combined.gff3' or '.gff').",
                          required=False)
    full_mode_args.add_argument("-pid", action="store", dest="pident", type=float, default=0.95,
                          help="Default 0.95: Pident threshold for clustering.",
                          required=False)
    full_mode_args.add_argument("-len_diff", action="store", dest="len_diff", type=float, default=0.80,
                          help="Default 0.80: Minimum length difference between clustered sequences - (-s) threshold for CD-HIT clustering.",
                          required=False)


    ###Partial-Mode Arguments
    partial_mode_args = parser.add_argument_group('Partial-Mode Arguments - Required when "-run_mode Partial" is used')
    partial_mode_args.add_argument('-cluster_file', action='store', dest='cluster_file',
                        help='Clustering output file containing CD-HIT, TSV or CSV Edge List',
                        required=False)

    ###Grouping Arguments
    grouping_args = parser.add_argument_group('Grouping Arguments - Use to fine-tune grouping of genes after clustering')
    grouping_args.add_argument('-reclustered', action='store', dest='reclustered', help='Clustering output file from secondary round of clustering',
                        required=False)
    grouping_args.add_argument('-seq_tag', action='store', dest='sequence_tag', default='StORF',
                        help='Default - "StORF": Unique identifier to be used to distinguish the second of two rounds of clustered sequences',
                        required=False)
    grouping_args.add_argument('-groups', action="store", dest='core_groups', default="99,95,15",
                        help='Default - (\'99,95,15\'): Gene family groups to use',
                        required=False)

    ###Output Arguments
    output_args = parser.add_argument_group('Output Parameters')
    output_args.add_argument('-w', action="store", dest='write_families', default=None,
                          help='Default - No output: Output sequences of identified families (provide levels at which to output "-w 99,95"'
                               ' - Must provide FASTA file with -fasta',
                          required=False)
    output_args.add_argument('-con', action="store", dest='con_core', default=None,
                          help='Default - No output: Output aligned and concatinated sequences of identified families - used for MSA (provide levels at which to output "-w 99,95"'
                               ' - Must provide FASTA file with -fasta',
                          required=False)
    output_args.add_argument('-original_fasta', action='store', dest='original_fasta',
                          help='FASTA file to use in conjunction with "-w" or "-con" when running in Partial Mode.',
                          required=False)
    output_args.add_argument('-gpa', action='store', dest='gene_presence_absence_out', help='Default - False: If selected, a Roary formatted gene_presence_absence.csv will be created - Required for Coinfinder and other downstream tools',
                        required=False)

    ### Misc Arguments
    misc = parser.add_argument_group('Misc')
    misc.add_argument('-verbose', action='store', dest='verbose', default=False, type=eval, choices=[True, False],
                        help='Default - False: Print out runtime messages',
                        required = False)
    misc.add_argument('-v', action='store_true', dest='version',
                        help='Default - False: Print out version number and exit',
                        required=False)

    options = parser.parse_args()

    ### Checking all required parameters are provided by user
    if options.run_mode == 'Full':
        required_full_mode = [options.input_type, options.input_dir, options.name_split, options.clust_tool,
                              options.pident, options.len_diff]
        if all(required_full_mode):
            # Proceed with the Full mode
            pass
        else:
            missing_options = [opt for opt in
                               ['input_type', 'input_dir', 'name_split', 'clust_tool', 'pident', 'len_diff'] if
                               not options.__dict__[opt]]
            print(f"Missing required options for Full mode: {', '.join(missing_options)}")
    elif options.run_mode == 'Partial':
        required_partial_mode = [options.cluster_file, ]
        if all(required_partial_mode):
            # Proceed with the Partial mode
            pass
        else:
            missing_options = [opt for opt in
                               ['cluster_file',] if
                               not options.__dict__[opt]]
            print(f"Missing required options for Partial mode: {', '.join(missing_options)}")

    if options.clust_tool == 'CD-HIT':
        clust_affix = '.clstr'
    elif options.clust_tool == 'TSV':
        clust_affix = '.tsv'
    elif options.clust_tool == 'CSV':
        clust_affix = '.csv'



    ###External tool checks:
    ##MAFFT
    if options.con_core == True:
        if is_tool_installed('mafft'):
            if options.verbose == True:
                print("mafft is installed. Proceeding with alignment.")
        else:
            exit("mafft is not installed. Please install mafft to proceed.")
    ##CD-HIT
    if options.clust_tool == 'CD-HIT':
        if is_tool_installed('cd-hit'):
            if options.verbose == True:
                print("cd-hit is installed. Proceeding with clustering.")
        else:
            exit("cd-hit is not installed. Please install cd-hit to proceed.")

    if options.write_families != None and options.original_fasta == False:
        exit("-fasta must br provided if -w is used")

    options.core_groups = options.core_groups + ',0'


    if options.cluster_file:
        options.cluster_file = fix_path(options.cluster_file)
    if options.reclustered:
        options.reclustered = fix_path(options.reclustered)
    if options.input_dir:
        options.input_dir = fix_path(options.input_dir)
    if options.output_dir:
        options.output_dir = fix_path(options.output_dir)

    output_path = os.path.abspath(options.output_dir)
    combined_out_file = os.path.join(output_path, "combined_sequences.fasta")
    clustering_output = os.path.join(output_path, 'clustering_' + options.clust_tool)


    if options.run_mode == 'Full':



        if options.input_type == 'separate':
            read_separate_files(options.input_dir, options.name_split, combined_out_file)
        else:
            read_combined_files(options.input_dir, options.name_split, combined_out_file)

        run_cd_hit(combined_out_file, clustering_output, options)
        class clustering_options:
            def __init__(self):
                self.cluster_format = options.clust_tool
                self.reclustered = options.reclustered
                self.sequence_tag = options.sequence_tag
                self.core_groups = '99,95,15,0'
                self.clusters = clustering_output + clust_affix
                self.gene_presence_absence_out = options.gene_presence_absence_out
                self.write_families = options.write_families
                self.con_core = options.con_core
                self.fasta = combined_out_file
                self.verbose = options.verbose

        clustering_options = clustering_options()

    elif options.run_mode == 'Partial':
        class clustering_options:
            def __init__(self):
                self.cluster_format = options.clust_tool
                self.reclustered = options.reclustered
                self.sequence_tag = options.sequence_tag
                self.core_groups = '99,95,15,0'
                self.clusters = options.cluster_file
                self.gene_presence_absence_out = options.gene_presence_absence_out
                self.write_families = options.write_families
                self.con_core = options.con_core
                self.fasta = options.original_fasta
                self.verbose = options.verbose

        clustering_options = clustering_options()




    cluster(clustering_options)

    print("Thank you for using PyamilySeq -- A detailed user manual can be found at https://github.com/NickJD/PyamilySeq\n"
          "Please report any issues to: https://github.com/NickJD/PyamilySeq/issues\n#####")

if __name__ == "__main__":
    main()
