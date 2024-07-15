# PyamilySeq
PyamilySeq (Family Seek) is a Python tool for clustering gene sequences into families based on sequence similarity identified by tools such as CD-HIT, DIAMOND or MMseqs2.
This work is an extension of the gene family / pangenome tool developed for the StORF-Reporter publication in NAR (https://doi.org/10.1093/nar/gkad814).

## Features

- **Clustering**: Supports input from CD-HIT formatted files as well as TSV and CSV Edge List formats.
- **Reclustering**: Allows for the addition of new sequences post-initial clustering.
- **Output**: Generates a gene 'Roary' presence-absence CSV formatted file for downstream analysis.

## Installation

PyamilySeq requires Python 3.6 or higher. Install dependencies using pip:

```bash
pip install PyamilySeq
```

## Usage - Menu
```
PyamilySeq_Species.py -h
usage: PyamilySeq_Species.py [-h] -c CLUSTERS -f {CD-HIT,CSV,TSV} [-w WRITE_FAMILIES] [-fasta FASTA] [-rc RECLUSTERED] [-st SEQUENCE_TAG]
                             [-groups CORE_GROUPS] [-gpa GENE_PRESENCE_ABSENCE_OUT] [-verbose {True,False}] [-v]

PyamilySeq v0.2.0: PyamilySeq Run Parameters.

Required Arguments:
  -c CLUSTERS           Clustering output file from CD-HIT, TSV or CSV Edge List
  -f {CD-HIT,CSV,TSV}   Which format to use (CD-HIT or Comma/Tab Separated Edge-List (such as MMseqs2 tsv output))

Output Parameters:
  -w WRITE_FAMILIES     Default - No output: Output sequences of identified families (provide levels at which to output "-w 99 95" - Must provide
                        FASTA file with -fasta
  -fasta FASTA          FASTA file to use in conjunction with "-w"

Optional Arguments:
  -rc RECLUSTERED       Clustering output file from secondary round of clustering
  -st SEQUENCE_TAG      Default - "StORF": Unique identifier to be used to distinguish the second of two rounds of clustered sequences
  -groups CORE_GROUPS   Default - ('99,95,90,80,15'): Gene family groups to use
  -gpa GENE_PRESENCE_ABSENCE_OUT
                        Default - False: If selected, a Roary formatted gene_presence_absence.csv will be created - Required for Coinfinder and other
                        downstream tools

Misc:
  -verbose {True,False}
                        Default - False: Print out runtime messages
  -v                    Default - False: Print out version number and exit

```

### Clustering Analysis

To perform clustering analysis:

```bash
python pyamilyseq.py -c clusters_file -f format
```

Replace `clusters_file` with the path to your clustering output file and `format` with one of: `CD-HIT`, `CSV`, or `TSV`.

### Reclustering

To add new sequences and recluster:

```bash
PyamilySeq -c clusters_file -f format --reclustered reclustered_file
```

Replace `reclustered_file` with the path to the file containing additional sequences.

## Output

PyamilySeq generates various outputs, including:

- **Gene Presence-Absence File**: This CSV file details the presence and absence of genes across genomes.
- **FASTA Files for Each Gene Family**:  

## Gene Family Groups

After analysis, PyamilySeq categorizes gene families into several groups:

- **First Core**: Gene families present in all analysed genomes initially.
- **Extended Core**: Gene families extended with additional sequences.
- **Combined Core**: Gene families combined with both initial and additional sequences.
- **Second Core**: Gene families identified only in the additional sequences.
- **Only Second Core**: Gene families exclusively found in the additional sequences.
