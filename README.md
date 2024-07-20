# PyamilySeq - !BETA!
**PyamilySeq** (Family Seek) is a Python tool for clustering gene sequences into families based on sequence similarity identified by tools such as CD-HIT, BLAST, DIAMOND or MMseqs2.
This work is an extension of the gene family / pangenome tool developed for the StORF-Reporter publication in NAR (https://doi.org/10.1093/nar/gkad814).

## Features
- **End-to-End**: PyamilySeq can take a directory of GFF+FASTA files, run CD-HIT for clustering and process the results.
- **Clustering**: Supports input from CD-HIT formatted files as well as CSV and TSV edge lists (-outfmt 6 from BLAST/DIAMOND).
- **Reclustering**: Allows for the addition of new sequences post-initial clustering.
- **Output**: Generates a gene 'Roary/Panaroo' formatted presence-absence CSV formatted file for downstream analysis.
  - Align representative sequences using MAFFT.
  - Output concatenated aligned sequences for downstream analysis.
  - Optionally output sequences of identified families.


### Installation
PyamilySeq requires Python 3.6 or higher. Install using pip:

```bash
pip install PyamilySeq
```

## Usage - Menu
```
usage: PyamilySeq.py [-h] -id INPUT_DIR -od OUTPUT_DIR -it {separate,combined} -ns NAME_SPLIT -pid PIDENT -ld LEN_DIFF -co CLUSTERING_OUT -ct {CD-HIT,BLAST,DIAMOND,MMseqs2} [-w WRITE_FAMILIES] [-con CON_CORE] [-fasta FASTA] [-rc RECLUSTERED] [-st SEQUENCE_TAG] [-groups CORE_GROUPS]
                     [-gpa GENE_PRESENCE_ABSENCE_OUT]
                     ...

PyamilySeq v0.4.0: PyamilySeq Run Parameters.

positional arguments:
  pyamilyseq_args       Additional arguments for PyamilySeq.

options:
  -h, --help            show this help message and exit

Required Arguments:
  -id INPUT_DIR         Directory containing GFF/FASTA files.
  -od OUTPUT_DIR        Directory for all output files.
  -it {separate,combined}
                        Type of input files: 'separate' for separate FASTA and GFF files, 'combined' for GFF files with embedded FASTA sequences.
  -ns NAME_SPLIT        Character used to split the filename and extract the genome name.
  -pid PIDENT           Pident threshold for CD-HIT clustering.
  -ld LEN_DIFF          Length difference (-s) threshold for CD-HIT clustering.
  -co CLUSTERING_OUT    Output file for initial clustering.
  -ct {CD-HIT,BLAST,DIAMOND,MMseqs2}
                        Clustering format for PyamilySeq.

Output Parameters:
  -w WRITE_FAMILIES     Default - No output: Output sequences of identified families (provide levels at which to output "-w 99,95" - Must provide FASTA file with -fasta
  -con CON_CORE         Default - No output: Output aligned and concatinated sequences of identified families - used for MSA (provide levels at which to output "-w 99,95" - Must provide FASTA file with -fasta
  -fasta FASTA          FASTA file to use in conjunction with "-w" or "-con"

Optional Arguments:
  -rc RECLUSTERED       Clustering output file from secondary round of clustering
  -st SEQUENCE_TAG      Default - "StORF": Unique identifier to be used to distinguish the second of two rounds of clustered sequences
  -groups CORE_GROUPS   Default - ('99,95,15'): Gene family groups to use
  -gpa GENE_PRESENCE_ABSENCE_OUT
                        Default - False: If selected, a Roary formatted gene_presence_absence.csv will be created - Required for Coinfinder and other downstream tools
```

### Example Run End-to-End - 'genomes' is a test-directory containing GFF files with ##FASTA at the bottom

```bash
PyamilySeq -id .../genomes -it combined -ns _combined.gff3 -pid 0.90 -ld 0.60 -co testing_cd-hit -ct CD-HIT -od .../testing
```

```Calculating Groups
Calculating Groups
Gene Groups:
first_core_99: 3103
first_core_95: 0
first_core_15: 3217
first_core_0: 4808
Total Number of Gene Groups (Including Singletons): 11128
```


