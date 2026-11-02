# PyamilySeq

[![PyPI Version](https://shields.io)](https://pypi.org/project/PyamilySeq/)
[![PyPI Downloads](https://shields.io)](https://pypi.org/project/PyamilySeq/)
[![Bioconda Version](https://shields.io)](https://anaconda.org/bioconda/pyamilyseq)
[![GitHub License](https://shields.io)](https://github.com)
[![GitHub Stars](https://shields.io)](https://github.com/NickJD/PyamilySeq)


PyamilySeq builds gene-family and pangenome summaries from CD-HIT, MMseqs2,
BLAST, or DIAMOND clustering results. It supports species- and genus-level
prevalence groups, staged reclustering, sequence export and alignment, visual
exploration, input validation, and genomic-location analysis.

PyamilySeq extends the gene-family workflow developed for StORF-Reporter. If
you use it in published work, please cite the
[PyamilySeq publication](https://doi.org/10.1093/nargab/lqaf198).

## Core Guarantees

PyamilySeq uses one canonical family ledger for every downstream output.

## Features
- **End-to-End**: PyamilySeq can take a directory of GFF+FASTA files, run CD-HIT for clustering and process the results.
- **Clustering input**: Supports input from CD-HIT formatted files as well as CSV and TSV node-edge lists (MMseqs2 and -outfmt 6 from BLAST/DIAMOND).
- **Reclustering**: Allows for the addition of new sequences post-initial clustering - Ensures continuity of contemporary clustering results and highlights impact of novel gene predictions.
- **'Genus Mode'**: Unlike other 'pangenome' tools, PyamilySeq can identify gene groups found across multiple genera as unique entities (see below).  
- **Output**: Generates a 'Roary/Panaroo' formatted presence-absence CSV formatted file for downstream analysis.
  - User-define species-/genus-wide gene groups - User has control over grouping parameters (core = 99/95% or min 6 genera etc).
  - Aligns representative sequences using MAFFT.
  - Output concatenated aligned sequences for tree building.
  - Optionally output sequences of each separate identified gene group.
  - Group-Splitter tool to split multi-copy gene groups.
  - Numerous additional tools to assist in the pre- and post-processing of data.


These checks address a common failure mode in pangenome pipelines: genes being
silently dropped because they do not appear in a cluster file.

## Installation

PyamilySeq requires Python 3.12 or newer.

```bash
python -m pip install PyamilySeq
```

For optional analysis libraries used by auxiliary workflows:

```bash
python -m pip install 'PyamilySeq[analysis]'
```

Full mode requires CD-HIT (`cd-hit-est` for DNA or `cd-hit` for amino-acid
input). Alignment output requires MAFFT. External command failures stop the run
and report the failing command instead of producing partial output.

## Quick Start

Partial mode consumes an existing clustering result and the authoritative FASTA:

```bash
PyamilySeq Partial \
  -clustering_format CD-HIT \
  -cluster_file clusters.clstr \
  -original_fasta genes.fasta.gz \
  -output_dir PyamilySeq_Output \
  -gene_group_table \
  -core_group_table \
  --visualize
```

Full mode extracts or combines input sequences, runs CD-HIT, and processes the
result:

```bash
PyamilySeq Full \
  -input_type combined \
  -input_dir genomes \
  -name_split_gff _combined.gff3 \
  -output_dir PyamilySeq_Output \
  --visualize
```

Run `PyamilySeq Full --help` or `PyamilySeq Partial --help` for the complete
option list. `--log` and `--log-dir` are accepted before or after the mode name.

## One-command supplementation workflow

`Supplementation` runs the full analysis on the original genomes, then clusters
the original representatives with the added genomes and writes the combined
analysis. Each input set can independently be a directory of combined GFF files,
a directory of paired GFF/FASTA files, or an extracted FASTA. Results are kept
under `original/` and `supplemented/` in the output directory.

```bash
PyamilySeq Supplementation \
  -original_input_type combined \
  -original_input_dir original_genomes \
  -original_name_split_gff .gff \
  -supplement_input_type combined \
  -supplement_input_dir added_genomes \
  -supplement_name_split_gff .gff \
  -output_dir supplementation_results \
  -sequence_type AA \
  -c 0.90 -s 0.80 -T 16 -M 8000
```

For extracted FASTA inputs, use `-original_input_type fasta
-original_input_fasta original_genes.fasta` and the corresponding
`-supplement_input_type fasta -supplement_input_fasta added_genes.fasta`.
Protein or nucleotide alphabets are detected automatically. New FASTA headers
must identify their genome using the configured sequence delimiter (by default,
`genome|gene`).

## Identifier And Taxon Rules

By default, sequence IDs use `genome|gene`, and the text before the first `|`
is treated as the genome. Change this with `--sequence-delimiter`.

For ambiguous identifiers or genus analyses, provide `--taxon-map`. The map may
be tab- or comma-separated and may contain either:

```text
identifier<TAB>taxon
identifier<TAB>genome<TAB>taxon
```

An identifier may be a complete sequence ID or a genome prefix. Exact sequence
mappings take precedence. In Genus mode, the historical fallback is the text
before the first underscore in the genome ID; an explicit map is recommended.

## Add genomes to a reclustered pangenome

Run Partial mode with the original cluster file, a FASTA containing the original
and added genes, and a second clustering result made from that same combined
FASTA. `--added-genomes` writes gene-level family links, a per-genome summary,
an SVG chart, and an HTML report. `--added-genomes` accepts comma-separated IDs
or a text file with one genome ID per line. `source_families` and `source_clusters` identify
the original families reached by each added gene; `only_Second` families have
no original family link.

```bash
PyamilySeq Partial \\
  -clustering_format CD-HIT \\
  -cluster_file original.clstr \\
  -reclustered original_plus_new.clstr \\
  -original_fasta original_plus_new_genes.fasta \\
  -output_dir extension_output \\
  --added-genomes new_isolate_1,new_isolate_2
```

## Query genomes in an existing pangenome

`PyamilySeq-Query` accepts a PyamilySeq/Roary presence-absence CSV and a text
file containing one genome column name per line. It reports family and gene
counts for core and accessory families (core defaults to presence in at least
99% of genomes), plus a family-level table with matching gene IDs.

```bash
PyamilySeq-Query --gpa PyamilySeq_Output/gene_presence_absence.csv \\
  --genomes requested_genomes.txt --output-dir genome_query --core-percent 99
```

## Clustering Inputs

### CD-HIT

Use the native `.clstr` file. Representatives and the final cluster are parsed
without relying on cluster numbering.

### MMseqs2

Use the standard two-column representative/member TSV. Repeated or
non-contiguous representative rows are combined into one family.

### BLAST Or DIAMOND

Two-column files are interpreted as explicit undirected family edges. Scored
tabular files use reciprocal-best edges by default to reduce weak-hit chaining.
For standard outfmt 6 input:

```bash
PyamilySeq Partial \
  -clustering_format BLAST \
  -cluster_file hits.tsv \
  -original_fasta genes.fasta \
  -output_dir output \
  --blast-min-identity 90 \
  --blast-max-evalue 1e-10
```

Use `--blast-columns` when the input column order differs. Query and subject
coverage filters require `qlen` and `slen` in that schema. The legacy connected
component interpretation remains available as
`--blast-graph-mode connected-components`, but can overmerge families through
weak bridge hits.

## Reclustering

Pass a second clustering result with `-reclustered`:

```bash
PyamilySeq Partial \
  -clustering_format CD-HIT \
  -cluster_file original.clstr \
  -reclustered original_and_novel.clstr \
  -original_fasta original_and_novel.fasta \
  -output_dir reclustered_output
```

Final families receive mutually exclusive provenance stages:

- `First`: an unchanged source family containing no second-round-only genes.
- `extended`: one source family plus new genes.
- `combined`: multiple source families merged by reclustering, with or without
  new genes.
- `only_Second`: genes that were absent from the first clustering, including
  singleton backfills.

The `source_families` and `source_clusters` columns preserve the reconciliation
history.

## Main Outputs

Every run writes:

- `gene_family_assignments.tsv`: the authoritative one-row-per-input-gene
  ledger, including genome, taxon, representative, stage, category, and source
  provenance.
- `family_summary.tsv`: one row per selected final family with prevalence,
  copy-number, representative, and provenance statistics.
- `gene_presence_absence.csv`: Roary/Panaroo-compatible family-by-genome table,
  unless `-no_gpa` is used.
- `pan_genome_reference.fa`: one deterministic representative per selected
  final family.
- `summary_statistics.txt`: category and assignment totals.
- `run_manifest.json`: inputs, parameters, stage counts, and the
  `exactly_once_assignment` audit result.

`-gene_group_table` and `-core_group_table` add compact TSV mappings.
`-write_groups`, `-write_individual_groups`, and `-align` produce sequence and
alignment outputs from the same final-family selection. Core alignment requires
exactly one sequence from every genome in each aligned family.

Species prevalence groups are configured with `-species_groups` (default
`99,95,15`). Genus groups use inclusive ranges such as `1-10`. Categories are
mutually exclusive and include every selected family, including the zero-percent
accessory bin.

## Single-Copy Filtering

`--single_copy_only` restricts analytical outputs to families with no
multi-copy occupied genome. `--single_copy_tolerance PERCENT` permits that
percentage of occupied genomes to be multi-copy. The complete
`gene_family_assignments.tsv` is always retained so filtering never erases the
input-to-family audit trail.

## Visual Outputs

Add `--visualize` to write `Visualisations/` with:

- category and prevalence SVG charts;
- presence/absence and copy-number heatmaps;
- family-size distribution;
- deterministic pangenome and strict-core accumulation curves with confidence
  ranges in `pangenome_accumulation.tsv`;
- matching TSV data tables;
- a searchable, filterable canvas-based pangenome explorer; and
- an HTML index linking all views.

Use `--visualize-heatmap-limit`, `--visualize-permutations`, and
`--visualize-seed` to control scale and reproducibility. Static plots show a
bounded number of rows, while TSV tables and the interactive explorer retain all
selected families.

## Input Validation

Validate clustering, FASTA, and optional GFF inputs independently:

```bash
PyamilySeq-Validate \
  --clusters clusters.clstr \
  --cluster-format CD-HIT \
  --fasta genes.fasta.gz \
  --gff-dir genomes \
  --outdir validation
```

The validator reports duplicate family membership, cluster IDs absent from the
FASTA, FASTA genes requiring singleton backfill, ambiguous GFF aliases, mapping
confidence, and mapping rate. Use `--fail-on error`, `warning`, or `never` for
automation. Reports are machine-readable TSV/JSON plus a text summary.

## Family Location Tracking

`Track-Family-Locations` maps final family members to GFF/GFF3 features and
describes their genomic context:

```bash
Track-Family-Locations \
  --clusters PyamilySeq_Output/gene_presence_absence.csv \
  --cluster-format gpa \
  --gff-dir genomes \
  --outdir family_locations \
  --minimum-mapping-rate 0.95
```

The tracker supports plain or gzip-compressed GFF/GFF3 files. Its default
`--feature-layer auto` prefers CDS loci and collapses parent/child and multipart
features to prevent duplicate biological loci. Alternative layers are `CDS`,
`gene`, and `all`.

Mapping is exact first, alias-aware second, and uniquely canonicalized last.
Each result records its method and confidence. Use
`--min-mapping-confidence`, `--strict-mapping`, and
`--minimum-mapping-rate` to enforce quality thresholds. `--id-map` resolves
known identifier differences explicitly.

Neighborhoods are transcription-relative by default, so upstream and downstream
remain biologically meaningful on the minus strand. Outputs record:

- signed family architectures and strand patterns;
- neighbor multiplicity and conserved-order metrics;
- intergenic base-pair gaps;
- contig length, topology, topology provenance, and relative-position basis;
- boundary censoring and circular wrapping; and
- per-member mapping diagnostics.

Use `--genomic-neighbors` for genomic left/right order. Circular contigs may be
declared with `--circular-contigs`, and combined GFF provenance can be supplied
with `--contig-genome-map`. `--unmapped-family-report-limit` and
`--unmapped-token-report-limit` bound terminal diagnostics while the complete
member-level report remains on disk.

Outputs include `family_locations.csv`, `family_summary.csv`,
`member_mapping.tsv`, `unmapped_members.tsv`, `location_tracking_report.txt`,
`neighbor_family_summary.tsv`, `neighborhood_architectures.tsv`, GraphML, and
HTML/SVG neighborhood views. Singleton families are included by default; use
`--skip-singletons` to omit them.

## Auxiliary Commands

- `Seq-Combiner`: extract and combine GFF/FASTA sequences.
- `Group-Splitter`: split multi-copy groups using sequence clustering.
- `Group-Extractor`: select groups from a presence/absence table.
- `Group-Summary`: summarize prevalence categories, including 0% and 100%.
- `Group-Compare`: compare family partitions independent of cluster row order.
- `Group-Sizes`: summarize family sizes.
- `Seq-Finder`: find matching group rows for sequence IDs.
- `Seq-Extractor`: extract requested full sequence identifiers.
- `Track-Genome-Order`: summarize genome contig architecture and strand-aware gene order from GFF/GFF3 files.
- `PyamilySeq-Cohort`: prepare reproducible genome cohorts, calculate Sourmash similarity groups, run a global PyamilySeq pangenome, and stratify breakpoints by similarity group.

Each command provides its full interface with `--help`.

Genome-order and cohort workflows are also available from the main command:
`PyamilySeq Analyze order ...`, `PyamilySeq Analyze location ...`, and
`PyamilySeq Analyze cohort ...`. The location workflow reports family
neighborhoods and writes its own HTML/SVG visual summaries.

`Track-Genome-Order --gff-dir genomes --outdir genome_order --no-plots` writes
tabular contig and gene-order reports. Add `--family-map` to compare shared
families and infer contig orientation against a reference genome. Install
`PyamilySeq[analysis]` to enable the figures. Orientation inference uses
single-copy family anchors and requires both `--min-orientation-anchors` and
`--min-orientation-tau`; ambiguous reference-contig matches remain unoriented.
The tool expects unique GFF feature IDs for multipart features and reports
partial family-map coverage and annotation/FASTA mismatches in `run_summary.txt`.

For reference-definition reports, `--subcohort-map` accepts a TSV with
`genome` and `subcohort` columns, and `--reference-map` accepts a TSV with
`subcohort` and `reference_genome` columns. Every genome and sub-cohort must be
assigned exactly once. Without these maps, all genomes are listed in one
`default` cohort and no reference is selected. The tool writes
`reference_candidates.tsv` and `reference_manifest.tsv`; candidate ranks are a
transparent review order, not a completeness score or an automatic selection.
Optional `--assembly-metadata` (`genome`, `assembly_level`, and/or
`completeness`) and `--similarity-matrix` inputs add candidate diagnostics.
Similarity values never determine sub-cohort boundaries. Supplying an explicit
reference map does not perform contig placement or endpoint recurrence.

### Cohort workflow for 100–1,000 genomes

The cohort workflow validates one matching GFF and FASTA per genome, writes a
manifest and managed input links, runs one global PyamilySeq pangenome, and uses
Sourmash similarity groups only to stratify downstream summaries. It does not
build separate pangenomes for each similarity group.

Sourmash is an external command-line requirement for cohort similarity
analysis. Bioconda installations should provide it as a package dependency;
PyPI installations must install the `sourmash` command separately with
`python -m pip install sourmash` before running the similarity stage.

```bash
PyamilySeq-Cohort prepare \
  --fasta-dir assemblies/fasta --gff-dir assemblies/gff \
  --metadata cohort_metadata.tsv --outdir ecoli_cohort

PyamilySeq-Cohort similarity \
  --manifest ecoli_cohort/genome_manifest.tsv --outdir ecoli_cohort \
  --min-similarity 0.55

PyamilySeq-Cohort pangenome \
  --manifest ecoli_cohort/genome_manifest.tsv --outdir ecoli_cohort \
  --threads 16 --memory-mb 16000

PyamilySeq-Cohort order \
  --manifest ecoli_cohort/genome_manifest.tsv --outdir ecoli_cohort \
  --reference-genome MG1655 --no-plots

PyamilySeq-Cohort locations \
  --manifest ecoli_cohort/genome_manifest.tsv --outdir ecoli_cohort
```

The Sourmash threshold is optional; omit it to build the signature collection
and similarity matrix before choosing a threshold. The threshold should be selected after
reviewing the similarity matrix. Groups default to deterministic greedy
complete-linkage assignment, so every pair in a group meets the threshold; use
`--grouping-method single_linkage` only when chained links are intended. The
edge table is retained to inspect these relationships. Similarity groups
stratify downstream summaries; the main pangenome remains global.
Order outputs include per-assembly endpoint evidence and a
`breakpoint_support_by_similarity_group.tsv` summary with callable-genome
denominators for each candidate region. `contig_end_gene_candidates.tsv`
lists genes close to physical contig ends or carrying partial/disruption
annotation flags; it labels review candidates and does not claim that a gene is
truncated or lost.
Pyrodigal calls are ab initio predictions, not curated annotations. For complete
chromosomes or plasmids, sequence endpoints may represent the chosen origin of
a closed replicon rather than an assembly break; use draft assemblies and
sequence-level validation before interpreting recurrent endpoints as hotspots.

For one multi-genome FASTA with BV-BRC/RefSeq-style bracketed sample labels,
split records by the text before the taxon ID. Assembly metrics can be prepared
without annotations:

```bash
PyamilySeq-Cohort prepare \\
  --combined-fasta BVBRC_genome_sequence.fasta --outdir ecoli_assemblies
```

To add gene-family and reference-projected analyses, predict CDSs with
Pyrodigal (install `PyamilySeq[annotation]`), use a local Bakta installation,
or supply matching per-genome GFF3 files with `--gff-dir`:

```bash
PyamilySeq-Cohort prepare \\
  --combined-fasta BVBRC_genome_sequence.fasta \\
  --annotate-with-pyrodigal --outdir ecoli_cohort
```

Custom headers can be grouped with `--genome-header-regex`, which must contain
a named `(?P<genome>...)` capture. The split manifest records the mapping from
header-derived genome labels to per-genome FASTAs. Contig-size reports can run
without annotations; reference-projected breakpoint analysis requires matching
GFF3 annotations and a family map.

After running preparation, Sourmash, pangenome, and order analysis, generate an
interactive, self-contained dashboard:

```bash
PyamilySeq Analyze cohort report \
  --manifest ecoli_cohort/genome_manifest.tsv --outdir ecoli_cohort
```

The report filters genomes by metadata, similarity group, N50, and sequence
record count, and can select neighbors above a chosen Jaccard similarity to a
seed genome. Selected genomes can be downloaded as a TSV, then turned into a
new managed cohort for independent downstream analysis:

```bash
PyamilySeq Analyze cohort subset \
  --manifest ecoli_cohort/genome_manifest.tsv \
  --genomes-file selected_genomes.tsv --outdir selected_cohort
```

The subset command creates its own manifest and staging directory, so rerunning
the pangenome stage on it discovers only the selected genomes. The dashboard's
family table is a descriptive audit of calls near contig ends; it does not
estimate false absences or establish that fragmentation caused family-call
differences. Endpoint analyses require FASTA or explicit GFF sequence lengths;
records whose length is inferred only from their last annotated feature are
excluded from endpoint-support counts. Reference-projected assembly ends are
screening observations, not statistical hotspots or biological breakpoints.

## Development

```bash
python -m pip install -e '.[dev]'
pytest
python -m compileall -q src scripts
```

Continuous integration tests Python 3.10 and 3.12. Bug reports and reproducible
examples are welcome at the
[issue tracker](https://github.com/NickJD/PyamilySeq/issues).
