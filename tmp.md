# PyamilySeq - !BETA!
**PyamilySeq** (Family Seek) is a Python tool for clustering gene sequences into families based on sequence similarity identified by tools such as CD-HIT, BLAST, DIAMOND or MMseqs2.
This work is an extension of the gene family / pangenome tool developed for the StORF-Reporter publication in NAR (https://doi.org/10.1093/nar/gkad814).

## Features

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
usage: PyamilySeq_Species.py [-h] -c CLUSTERS -f {CD-HIT,CSV,TSV} [-w WRITE_FAMILIES] [-con CON_CORE] [-fasta FASTA] [-rc RECLUSTERED] [-st SEQUENCE_TAG]
                             [-groups CORE_GROUPS] [-gpa GENE_PRESENCE_ABSENCE_OUT] [-verbose {True,False}] [-v]

PyamilySeq v0.3.0: PyamilySeq Run Parameters.

Required Arguments:
  -c CLUSTERS           Clustering output file from CD-HIT, TSV or CSV Edge List
  -f {CD-HIT,CSV,TSV}   Which format to use (CD-HIT or Comma/Tab Separated Edge-List (such as MMseqs2 tsv output))

Output Parameters:
  -w WRITE_FAMILIES     Default - No output: Output sequences of identified families (provide levels at which to output "-w 99,95" - Must provide FASTA file
                        with -fasta
  -con CON_CORE         Default - No output: Output aligned and concatinated sequences of identified families - used for MSA (provide levels at which to
                        output "-w 99,95" - Must provide FASTA file with -fasta
  -fasta FASTA          FASTA file to use in conjunction with "-w" or "-con"

Optional Arguments:
  -rc RECLUSTERED       Clustering output file from secondary round of clustering
  -st SEQUENCE_TAG      Default - "StORF": Unique identifier to be used to distinguish the second of two rounds of clustered sequences
  -groups CORE_GROUPS   Default - ('99,95,15'): Gene family groups to use
  -gpa GENE_PRESENCE_ABSENCE_OUT
                        Default - False: If selected, a Roary formatted gene_presence_absence.csv will be created - Required for Coinfinder and other
                        downstream tools

Misc:
  -verbose {True,False}
                        Default - False: Print out runtime messages
  -v                    Default - False: Print out version number and exit


```
### Example Run

```bash
python pyamilyseq.py -c ..test_data/CD-HIT/combined_Ensmbl_pep_CD_90_60.clstr -f CD-HIT
```

```Calculating Groups
Gene Groups:
first_core_99: 3099
first_core_95: 0
first_core_15: 3208
first_core_0: 4988
Total Number of Gene Groups (Including Singletons): 11295
Outputting gene_presence_absence file
Thank you for using PyamilySeq -- A detailed user manual can be found at https://github.com/NickJD/PyamilySeq
Please report any issues to: https://github.com/NickJD/PyamilySeq/issues
#####
Done
```



### Command Line Arguments

#### Required Arguments
- `-c <clusters_file>`: Clustering output file from CD-HIT, TSV, or CSV Edge List.
- `-f <format>`: Format of the clustering output file. Choices are 'CD-HIT', 'CSV', or 'TSV'.

#### Output Parameters
- `-w <levels>`: Output sequences of identified families at specified levels. Must provide a FASTA file with `-fasta`. Example: `-w 99,95`
- `-con <levels>`: Output aligned and concatenated sequences of identified families for MSA at specified levels. Must provide a FASTA file with `-fasta`. Example: `-con 99,95`
- `-fasta <file>`: FASTA file to use in conjunction with `-w` or `-con`.

#### Optional Arguments
- `-rc <file>`: Clustering output file from a secondary round of clustering.
- `-st <tag>`: Unique identifier to distinguish the second of two rounds of clustered sequences. Default is 'StORF'.
- `-groups <groups>`: Gene family groups to use. Default is '99,95,15'.
- `-gpa`: If selected, a Roary formatted `gene_presence_absence.csv` will be created. Required for Coinfinder and other downstream tools.

#### Miscellaneous
- `-verbose`: Print out runtime messages. Default is `False`.
- `-v`: Print out version number and exit.

### Advanced Example

```bash
python pyamilyseq.py -c clusters.tsv -f TSV -w 99,95 -fasta sequences.fasta -verbose True
```

This command processes the clustering output file `clusters.tsv` in TSV format, outputs sequences of identified families at levels 99 and 95, uses `sequences.fasta` as the input FASTA file, and enables verbose runtime messages.
