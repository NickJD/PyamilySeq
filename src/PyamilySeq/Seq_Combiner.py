import argparse


try:
    from .Constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from Constants import *
    from utils import *



def main():
    parser = argparse.ArgumentParser(description='Seq-Combiner ' + PyamilySeq_Version + ': Seq-Combiner Run Parameters.')
    ### Required Arguments
    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-input_dir', action='store', dest='input_dir',
                          help='Directory location where the files are located.',
                          required=True)
    required.add_argument("-input_type", action="store", dest="input_type", choices=['separate', 'combined'],
                          help="Type of input files: 'separate' for separate FASTA and GFF files,"
                             " 'combined' for GFF files with embedded FASTA sequences.",
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
    options = parser.parse_args()

    output_path = os.path.abspath(options.output_dir)
    combined_out_file = os.path.join(output_path, options.output_file)

    if options.input_type == 'separate':
        read_separate_files(options.input_dir, options.name_split, )
    else:
        read_combined_files(options.input_dir, options.name_split, combined_out_file)

if __name__ == "__main__":
    main()
