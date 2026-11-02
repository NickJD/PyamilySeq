# PyamilySeq annotation metadata and deterministic-output update

This update applies to the uploaded PyamilySeq source modules. It makes GFF extraction deterministic, adds strict sequence validation to all supported extraction paths, writes gene annotations to the standard `gene_presence_absence.csv` metadata columns, and makes intermediate FASTA filenames identify their role in the workflow.

## 1. GFF annotation extraction

When PyamilySeq extracts selected features (normally `CDS`) from GFF/GFF3/GTF input, it now also writes a compressed sidecar:

```text
gene_annotations.tsv.gz
```

The sidecar contains:

| Column | Meaning |
|---|---|
| `sequence_id` | The exact FASTA identifier used by PyamilySeq (`genome|feature_id`) |
| `gene_name` | GFF `gene`, falling back to `gene_name`; blank when neither is present |
| `annotation` | GFF `product`, falling back to `Name`, then `Note` |
| `source_gff` | Source GFF basename |

Attribute values are URL-decoded, so values such as `%2C` are restored to commas. A protein description is not guessed to be a gene symbol: when no `gene`/`gene_name` attribute exists, **Non-unique Gene name** remains blank.

The shared Species/Genus workflow uses this sidecar to populate the existing first three columns of `gene_presence_absence.csv`:

- `Gene`: stable PyamilySeq family identifier (unchanged)
- `Non-unique Gene name`: gene symbol from the GFF
- `Annotation`: functional product annotation from the GFF

There is one GPA row per family, but annotations can differ between family members. PyamilySeq therefore chooses the most frequent non-empty value independently for gene name and annotation. If there is a tie, values are ordered lexically and the first is selected. If no family member has a value, the field is blank. This rule is deterministic and does not change clustering or family assignment.

The metadata sidecar is compressed while it is written. The GPA writer imports sidecars into a temporary, disk-backed SQLite lookup rather than building a Python dictionary containing millions of annotations in RAM. The temporary database is closed and removed after GPA writing.

### Prepared FASTA input

FASTA files do not contain GFF attributes. When a prepared FASTA is used, annotations can be supplied explicitly or auto-detected where documented:

- Full: `-gene_annotations PATH`
- Partial: `-gene_annotations PATH`
- Incremental original FASTA: `-gene_annotations PATH`
- Incremental new FASTA: `-new_gene_annotations PATH`
- Supplementation original FASTA: `-original_gene_annotations PATH`
- Supplementation supplemental FASTA: `-supplement_gene_annotations PATH`

For Full, Partial, and Incremental original FASTA inputs, the application also looks for `gene_annotations.tsv.gz` or `gene_annotations.tsv` beside the FASTA. GFF inputs generate their sidecars automatically; no additional option is needed.

## 2. Deterministic extraction and sequence validation

The combined-GFF, paired-GFF/FASTA, and per-genome FASTA extraction helpers now sort input filenames by basename before processing. This removes dependence on filesystem `glob` enumeration order and makes extraction order reproducible when the input files are unchanged.

The GFF extraction paths validate feature coordinates, contig references, strands, and extracted nucleotide sequences. Strict IUPAC DNA validation accepts `ACGTRYSWKMBDHVN` (upper- or lowercase) and rejects other characters rather than silently replacing them. Input text is decoded as strict UTF-8.

FASTA and annotation outputs are written to temporary files beside their destinations and only promoted after successful completion. A failed extraction therefore does not replace a previously complete output with a partial file.

## 3. Output filename changes

New runs use filenames that distinguish source data from intermediate reconciliation inputs:

| Previous name | New name |
|---|---|
| `original/genes_dna.fasta` | `original/original_genes_dna.fasta` |
| `original/original_sequences_aa.fasta` | `original/original_clustering_sequences_aa.fasta` |
| `supplemented/genes_dna.fasta` | `supplemented/supplemental_genes_dna.fasta` |
| `incremental_new_genes_dna.fasta` | `supplemental_new_genes_dna.fasta` |
| `incremental_new_genes_aa.fasta` | `supplemental_new_genes_aa.fasta` |
| `incremental_representatives_dna.fasta` | `reclustering_input_dna.fasta` |
| `incremental_representatives_aa.fasta` | `reclustering_input_aa.fasta` |
| `incremental_combined_sequences_dna.fasta` | `supplemented_full_sequences_dna.fasta` |
| Full-mode `combined_sequences_dna.fasta` | `input_genes_dna.fasta` |

`gene_annotations.tsv.gz` is written alongside each GFF-extracted input set. In a Supplementation run, the original and supplemental sidecars remain separate; both are supplied to the final reconciliation GPA writer.

These are names for **new runs**. Existing output directories and their manifest paths are not renamed or rewritten by this update.

## 4. Incremental and Supplementation behaviour

The underlying workflow is unchanged:

1. The original set is clustered and analysed as its own baseline.
2. Original CD-HIT representatives are combined with new genes for reclustering.
3. The full original sequence set is retained and the new genes are appended for final reconciliation.
4. Family provenance and original-core lineage reporting continue to use the existing reconciliation path.

The new filenames make these distinct FASTAs easier to identify. Python return-dictionary keys from `prepare_incremental_inputs()` remain stable for callers; only the generated filenames have changed.

## 5. Other source modules

The `order_analysis.py` module had a mis-indented reference-map preset block in the uploaded source. The block is corrected so optional `minimap2_preset` values from the reference map are validated and retained as intended. Direct-execution import fallbacks were also added to `cohort_analysis.py` and `nucleotide_placement.py`, while package-style imports remain supported.

The Species and Genus wrappers continue to call the shared `run_pangenome_workflow()`. The metadata change is therefore handled by the same GPA writer in both modes; it does not alter genus grouping or family assignment.

## 6. Validation performed

The update was tested with:

- Python compilation (`python -m compileall -q .`)
- Eight unit/regression tests covering deterministic file order, combined GFF extraction, paired Prokka-like GFF/FASTA extraction, compressed annotation sidecars, strict validation and atomic-output preservation, incremental output naming, and Species/Genus GPA metadata writing
- Import of all 17 uploaded source modules
- CLI `--help` smoke tests for PyamilySeq Full, Partial, Incremental, Supplementation and Analyze, plus the order-analysis, cohort-analysis, genome-report, supplementation-impact, location-tracking and validation subtools
- The supplied Bakta GFF (`GCF_018986875.1.gff3.gz`): 2,225 CDS records extracted, 2,225 product annotations recovered, and 750 gene symbols recovered
- A small end-to-end shared-workflow test in both Species and Genus modes using the Bakta-derived FASTA/annotations and synthetic singleton CD-HIT clusters; each produced 2,225 families with all sequences assigned exactly once, and the GPA contained gene symbols and product annotations

The Bakta workflow check isolates extraction, annotation sidecar loading and GPA writing; it does not represent a full CD-HIT reclustering run. No full 5,889-genome production run was performed as part of this code update.

## 7. Applying the update

The ZIP contains the uploaded source modules, the regression tests, and this README. Review the changes against your working tree, then replace the corresponding modules in your project. Keep your existing packaging metadata and project-specific files, which were not part of the uploaded code set.

Run the tests from the source package directory:

```bash
python -m compileall -q .
PYTHONPATH=. python -m unittest discover -s tests -v
```

For a production-scale rerun, use a fresh output directory so files from the previous naming convention are not mistaken for outputs from the updated run.
