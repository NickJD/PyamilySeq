---
title: 'PyamilySeq (Family Seek): a Python-based tool for clustering gene sequences into families based on sequence similarity identified by tools such as CD-HIT, DIAMOND, BLAST or MMseqs2.'
tags:
  - Python
  - bioinformatics
  - gene families
  - clustering
authors:
  - name: Nicholas J. Dimonaco
    orcid: 0000-0002-3808-206X
    affiliation: 1, 2, 3, 4
affiliations:
  - name: Department of Medicine, McMaster University, Hamilton, Canada
    index: 1
  - name: Farncombe Family Digestive Health Research Institute, McMaster University, Hamilton, Canada
    index: 2
  - name: Department of Computer Science, Aberystwyth University, Wales, UK
    index: 3
  - name: School of Biological Sciences, Queen's University Belfast, Northern Ireland, UK
    index: 4
date: 12 May 2024
bibliography: paper.bib
---

# Summary
PyamilySeq (Family Seek) is a Python-based tool for grouping gene sequences into families based on sequence similarity identified by tools such as CD-HIT[@fu2012cd], DIAMOND[@buchfink2021sensitive], BLAST[@camacho2009blast+] or MMseqs2[@steinegger2017mmseqs2]. This work is an extension of the gene family / pangenome tool developed for the StORF-Reporter[@dimonaco2023storf] publication in NAR (https://doi.org/10.1093/nar/gkad814). 


# Statement of need
PyamilySeq is a user-friendly tool for detecting gene groups or 'families'. Unlike similar tools that focus on a rigid definition of a species pangenome, PyamilySeq gives users more control over the process. It works with leading clustering tools to group DNA and amino acid sequences, including non-coding elements like tRNA and ncRNA. Plus, it can also work with re-clustered sequences from various annotation sources, as seen in projects like StORF-Reporter. PyamilySeq fills a gap in existing pangenome/gene family analysis tools by catering to these specific needs.


# Method
PyamilySeq works by taking the output from a sequence clustering tool such as CD-HIT and using user-defined parameters, groups the sequences into non-overlapping user-defined sets. This is a similar approach to pangenome analysis tools such as Roary[@page2015roary] or Panaroo[@tonkin2020producing] but does not utilise gene synteny or 'care' about the types of sequences (CDS, tRNA, ncRNA etc) that it works with. 

For species-wide clustering, sequence IDs must be formatted with the name of the genome separated from the rest of the ID (gene name) by a '|'.
We provide a set of tools to allow for the correct formatting of sequences before clustering is performed with one of the compatible tools. 


Below is the example output of 10 <em>Escherichia coli (E. coli)</em> genomes from Ensembl Bacteria[@yates2022ensembl], clustered with CD-HIT and analysed by PyamilySeq.


Table 1:
This table reports an example output of using PyamilySeq on a set of 10 <em>E. coli</em> genomes with the canonical annotations from Ensembl and additional annotations from StORF-Reporter.

| **Cluster type**       | **Core** | **Soft-Core** | **Accessory** |
|------------------------|----------|---------------|---------------|
| Ensembl-Only           | 3,370    | 210           | 2,728         |
| Ensembl-StORF          | 10       | 3             | 33            |
| StORF-Combined-Ensembl | 0        | 0             | 2             |
| StORF                  | 0        | 0             | 54            |
| StORF-Only             | 348      | 267           | 3,224         |


Table 2:
This table.. Is an example of using PyamilySeq for the cross-genera analysis - Will update

| **Cluster Type**       | **1 Genus** | **2 Genera** | **3 Genera** | **4 Genera** | **5 Genera** | **6 Genera** | **>6 Genera** |
|------------------------|-------------|--------------|--------------|--------------|--------------|--------------|---------------|
| Ensembl-Only           | 2,296,340   | 203,993      | 19,566       | 4861         | 2,36         | 990          | 2,388         |
| Ensembl-StORF          | 0           | 0            | 177          | 64           | 35           | 120          | 44            |
| StORF-Combined-Ensembl | 0           | 0            | 0            | 0            | 1            | 2            | 5             |
| StORF                  | 98,671      | 2,573        | 199          | 63           | 27           | 25           | 27            |
| StORF-Only             | 1,328,805   | 70,111       | 3483         | 964          | 359          | 244          | 368           |




# Availability
PyamilySeq is distributed on PyPi under the [GNU Lesser General Public License](https://www.gnu.org/licenses/lgpl-3.0). Source code is available on [Github](https://github.com/NickJD/PyamilySeq) and feature requests/bug notices can be made [here](https://github.com/NickJD/PyamilySeq/issues).


# Acknowledgments


# Funding
N.J.D was supported by Farncombe Digestive Health Disease Institute (McMaster University) and the Weston Family Microbiome Initiative.


# References
