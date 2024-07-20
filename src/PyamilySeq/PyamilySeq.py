import argparse
import os
import glob
import subprocess
from PyamilySeq_Species import *


def read_and_rename_files(input_dir, split_char, combined_out):
    with open(combined_out, 'w') as combined_out_file:
        for file in glob.glob(os.path.join(input_dir, '*.fasta')) + glob.glob(os.path.join(input_dir, '*.gff')):
            genome_name = os.path.basename(file).split(split_char)[0]
            with open(file, 'r') as genome:
                for line in genome:
                    if line.startswith('#'):
                        continue
                    elif line.startswith('>'):
                        line = line.split(' ')[0]  # Take only the header line
                        line = line.replace('>', '>' + genome_name + '|')
                        combined_out_file.write(line + '\n')
                    else:
                        combined_out_file.write(line)


def run_cd_hit(input_file, output_file, cdhit_threshold):
    cdhit_command = [
        'cd-hit',
        '-i', input_file,
        '-o', output_file,
        '-c', str(cdhit_threshold)
    ]
    subprocess.run(cdhit_command)


def run_pyamilyseq(clusters, output_format, fasta_file, pyamilyseq_args):
    pyamilyseq_command = [
                             'pyamilyseq',  # Assuming `pyamilyseq` is executable in PATH
                             '-c', clusters,
                             '-f', output_format,
                         ] + pyamilyseq_args
    if fasta_file:
        pyamilyseq_command += ['-fasta', fasta_file]

    subprocess.run(pyamilyseq_command)


def main():
    parser = argparse.ArgumentParser(
        description='PyamilySeq ' + PyamilySeq_Version + ': PyamilySeq Run Parameters.')
    required = parser.add_argument_group('Required Arguments')
    required.add_argument("-id", action="store", dest="input_dir",
                          help="Directory containing GFF/FASTA files.")
    required.add_argument("-sc", action="store", dest="split_char",
                          help="Character used to split the filename and extract the genome name.")
    required.add_argument("-pid", action="store", dest="pident", type=float,
                          help="Pident threshold for CD-HIT clustering.")
    required.add_argument("cdhit_output", help="Output file for CD-HIT clustering.")
    required.add_argument("pyamilyseq_output_format", choices=['CD-HIT', 'CSV', 'TSV'],
                        help="Output format for PyamilySeq.")
    parser.add_argument("pyamilyseq_args", nargs=argparse.REMAINDER, help="Additional arguments for PyamilySeq.")
    args = parser.parse_args()

    combined_out_file = "combined_sequences.fasta"

    # Step 1: Read and rename sequences from files
    read_and_rename_files(args.input_dir, args.split_char, combined_out_file)

    # Step 2: Run CD-HIT on the renamed sequences
    run_cd_hit(combined_out_file, args.cdhit_output, args.cdhit_threshold)

    # Step 3: Run PyamilySeq with the CD-HIT output
    run_pyamilyseq(args.cdhit_output, args.pyamilyseq_output_format, combined_out_file, args.pyamilyseq_args)


if __name__ == "__main__":
    main()
