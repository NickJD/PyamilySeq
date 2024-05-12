
import argparse
import gzip
import glob

def combine_files(files, split, glob_location, combined_out):
    count = 0

    for file in glob.glob(glob_location + '/' + files):
        count += 1
        try:
            with gzip.open(file, 'rb') as genome:

                for line in genome:
                    if line.startswith(b'#'):
                        continue
                    elif line.startswith(b'>'):
                        genome_name = bytes(file.split(split)[0].split('/')[-1], 'utf-8')
                        line = line.split(b' ')[0]
                        line = line.replace(b'>', b'>' + genome_name + b'|')
                        combined_out.write(line.decode('utf-8')+'\n')
                    else:
                        combined_out.write(line.decode('utf-8'))
        except gzip.BadGzipFile:
            with open(file, 'r') as genome:

                for line in genome:
                    if line.startswith('#'):
                        continue
                    elif line.startswith('>'):
                        genome_name = file.split(split)[0].split('/')[-1]
                        line = line.replace('>', '>' + genome_name + '|')
                        combined_out.write(line)
                    else:
                        combined_out.write(line)

def main():
    parser = argparse.ArgumentParser(description="Combine gzipped fasta files.")
    parser.add_argument("files", help="File pattern to match within the specified directory.")
    parser.add_argument("split", help="String used to split the file path and extract the genome name.")
    parser.add_argument("glob_location", help="Directory location where the files are located.")
    parser.add_argument("combined_out", help="Output file where the combined data will be written.")
    args = parser.parse_args()

    with open(args.combined_out, 'w') as combined_out:
        combine_files(args.files, args.split, args.glob_location, combined_out)

if __name__ == "__main__":
    main()
