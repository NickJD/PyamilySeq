import argparse


try:
    from .constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *
    from utils import *



def main():
    parser = argparse.ArgumentParser(description='PyamilySeq ' + PyamilySeq_Version + ': Seq-Combiner - A tool to extract sequences from GFF/FASTA files and prepare them for PyamilySeq.')
    ### Required Arguments
    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-input_dir', action='store', dest='input_dir',
                          help='Directory location where the files are located.',
                          required=True)
    required.add_argument('-input_type', action='store', dest='input_type', choices=['separate', 'combined', 'fasta'],
                          help='Type of input files: "separate" for separate FASTA and GFF files,'
                             ' "combined" for GFF files with embedded FASTA sequences and "fasta" for combining multiple '
                               'FASTA files together.',
                          required=True)
    required.add_argument("-name_split_gff", action="store", dest="name_split_gff",
                          help="Substring used to split the filename and extract the genome name ('_combined.gff3' or '.gff'). - Not needed with -input_type fasta",
                          required=False)
    required.add_argument("-name_split_fasta", action="store", dest="name_split_fasta",
                          help="Substring used to split filenames and extract genome names for fasta files if named differently to paired gff files (e.g., '_dna.fasta').",
                          required=False)
    required.add_argument("-output_dir", action="store", dest="output_dir",
                          help="Directory for all output files.",
                          required=True)
    required.add_argument("-output_name", action="store", dest="output_file",
                          help="Output file name.",
                          required=True)

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-gene_ident', action='store', dest='gene_ident', default='CDS',
                          help='Default - "CDS": Identifier used for extraction of sequences such as "misc_RNA,gene,mRNA,CDS,rRNA,tRNA,tmRNA,CRISPR,ncRNA,regulatory_region,oriC,pseudo"'
                               ' - Not compatible with "fasta" input mode.',
                          required=False)
    optional.add_argument('-translate', action='store_true', dest='translate', default=None,
                          help='Default - False: Translate extracted sequences to their AA counterpart? - appends _aa.fasta to given output_name',
                          required=False)
    misc = parser.add_argument_group('Misc Arguments')
    misc.add_argument("-v", "--version", action="version",
                      version=f"PyamilySeq: Seq-Combiner version {PyamilySeq_Version} - Exiting",
                      help="Print out version number and exit")

    options = parser.parse_args()


    if options.input_type == 'separate' and options.name_split_gff is None:
        print("Please provide a substring to split the filename and extract the genome name.")
        exit(1)
    if options.input_type == 'combined' and options.name_split_gff is None:
        print("Please provide a substring to split the filename and extract the genome name.")
        exit(1)
    if options.input_type == 'fasta' and options.name_split_fasta is None:
        print("Please provide a substring to split the filename and extract the genome name.")
        exit(1)

    output_path = os.path.abspath(options.output_dir)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    #output_file = options.output_file + '.fasta'
    if os.path.exists(os.path.join(output_path, options.output_file)):
        print(f"Output file {options.output_file} already exists in the output directory. Please delete or rename the file and try again.")
        exit(1)

    combined_out_file = os.path.join(output_path, options.output_file )

    if options.input_type == 'separate':
        read_separate_files(options.input_dir, options.name_split_gff, options.name_split_fasta, options.gene_ident, combined_out_file, options.translate, True)
    elif options.input_type == 'combined':
        read_combined_files(options.input_dir, options.name_split_gff, options.gene_ident, combined_out_file, options.translate, True)
    elif options.input_type == 'fasta':
        read_fasta_files(options.input_dir, options.name_split_fasta, combined_out_file, options.translate, True)

if __name__ == "__main__":
    main()
