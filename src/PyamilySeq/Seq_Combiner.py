import argparse


try:
    from .Constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from Constants import *
    from utils import *



def main():
    parser = argparse.ArgumentParser(description='Seq-Combiner ' + PyamilySeq_Version + ': A tool to extract sequences from GFF files.')
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
    required.add_argument("-name_split", action="store", dest="name_split",
                          help="substring used to split the filename and extract the genome name ('_combined.gff3' or '.gff').",
                          required=True)
    required.add_argument("-output_dir", action="store", dest="output_dir",
                          help="Directory for all output files.",
                          required=True)
    required.add_argument("-output_name", action="store", dest="output_file",
                          help="Output file name.",
                          required=True)
    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-gene_ident', action='store', dest='gene_ident', default='CDS',
                          help='Identifier used for extraction of sequences such as "misc_RNA,gene,mRNA,CDS,rRNA,tRNA,tmRNA,CRISPR,ncRNA,regulatory_region,oriC,pseudo"'
                               ' - Not compatible with "fasta" input mode.',
                          required=False)
    misc = parser.add_argument_group('Misc Arguments')
    misc.add_argument('-v', action='store_true', dest='version',
                         help='Print out version number and exit',
                         required=False)

    options = parser.parse_args()

    if options.version:
        sys.exit(PyamilySeq_Version)

    output_path = os.path.abspath(options.output_dir)
    combined_out_file = os.path.join(output_path, options.output_file)

    if options.input_type == 'separate':
        read_separate_files(options.input_dir, options.name_split, options.gene_ident, combined_out_file)
    elif options.input_type == 'combined':
        read_combined_files(options.input_dir, options.name_split, options.gene_ident, combined_out_file)
    elif options.input_type == 'fasta':
        read_fasta_files(options.input_dir, options.name_split, combined_out_file)

if __name__ == "__main__":
    main()
