try:
    from .PyamilySeq_Species import cluster as species_cluster
    from .PyamilySeq_Genus import cluster as genus_cluster
    from .incremental import prepare_incremental_inputs, IncrementalInputError
    from .constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from PyamilySeq_Species import cluster as species_cluster
    from PyamilySeq_Genus import cluster as genus_cluster
    from incremental import prepare_incremental_inputs, IncrementalInputError
    from constants import *
    from utils import *

import traceback
import sys
import argparse
from types import SimpleNamespace


def _normalise_clustering_format(value):
    if value is None:
        return None
    v = value.strip()
    low = v.lower()
    if low in ('cd-hit', 'cdhit', 'cd_hit', 'cd-hit-est', 'cdhit-est'):
        return 'CD-HIT'
    if low.startswith('mmseq'):
        return 'MMseqs'
    if low.startswith('blast'):
        return 'BLAST'
    raise argparse.ArgumentTypeError(f"invalid clustering format: '{value}'. valid: CD-HIT, MMseqs, BLAST")




def _gene_annotation_sidecar_for_fasta(fasta_path):
    """Return a neighbouring gene-annotation sidecar, preferring gzip."""
    if not fasta_path:
        return None
    directory = os.path.dirname(os.path.abspath(fasta_path))
    candidates = (
        os.path.join(directory, "gene_annotations.tsv.gz"),
        os.path.join(directory, "gene_annotations.tsv"),
    )
    return next((path for path in candidates if os.path.isfile(path)), None)


def run_cd_hit(options, input_file, clustering_output, clustering_mode):
    logger = logging.getLogger("PyamilySeq.PyamilySeq")
    cdhit_command = [
        clustering_mode,
        '-i', input_file,
        '-o', clustering_output,
        '-c', f"{float(options.pident):.2f}",
        '-s', f"{float(options.len_diff):.2f}",
        '-T', str(options.threads),
        '-M', str(options.mem),
        '-d', "0",
        '-g', str(options.fast_mode),
        '-sc', "1",
        '-sf', "1"
    ]
    logger.debug("CD-HIT command: %s", " ".join(cdhit_command))
    run_mode = "fast" if getattr(options, "fast_mode", 1) == 0 else "accurate"
    logger.info(
        "Starting %s in %s mode (threads=%s; memory=%s MB; input=%s).",
        clustering_mode.upper(), run_mode, options.threads, options.mem, input_file,
    )
    try:
        if options.verbose:
            ret = subprocess.run(cdhit_command)
        else:
            ret = subprocess.run(cdhit_command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if ret.returncode != 0:
            raise RuntimeError(f"{clustering_mode} returned non-zero exit code {ret.returncode}")
        logger.info("cd-hit completed successfully: %s", clustering_output)
    except Exception as e:
        logger.exception("Failed to run cd-hit: %s", e)
        raise


def _run_supplementation(options, output_path, groups_to_use, logger):
    """Run the original and supplemented pangenomes with explicit input naming."""
    os.makedirs(output_path, exist_ok=True)
    original_dir = os.path.join(output_path, "original")
    supplemented_dir = os.path.join(output_path, "supplemented")
    os.makedirs(original_dir, exist_ok=True)
    os.makedirs(supplemented_dir, exist_ok=True)

    def extract_input(prefix, input_type, input_dir, input_fasta, gff_suffix, fasta_suffix, target):
        supplied_annotations = getattr(options, f"{prefix}_gene_annotations", None)
        if input_type == "fasta":
            if not input_fasta or not os.path.isfile(input_fasta):
                raise ValueError(f"Input FASTA does not exist: {input_fasta}")
            if supplied_annotations is None:
                supplied_annotations = _gene_annotation_sidecar_for_fasta(input_fasta)
            if supplied_annotations and not os.path.isfile(supplied_annotations):
                raise ValueError(f"Gene-annotation file does not exist: {supplied_annotations}")
            return input_fasta, supplied_annotations

        dna_name = "original_genes_dna.fasta" if prefix == "original" else "supplemental_genes_dna.fasta"
        dna_path = os.path.join(target, dna_name)
        annotation_path = os.path.join(target, "gene_annotations.tsv.gz")
        if input_type == "combined":
            read_combined_files(
                input_dir, gff_suffix, options.gene_ident, dna_path, False, False,
                progress_logger=logger, annotation_out=annotation_path,
            )
        else:
            read_separate_files(
                input_dir, gff_suffix, fasta_suffix, options.gene_ident,
                dna_path, False, False, annotation_out=annotation_path,
            )
        return dna_path, annotation_path

    def clustering_input(fasta_path, target, prefix):
        is_dna = detect_sequence_type(fasta_path)
        if options.sequence_type == "DNA":
            if not is_dna:
                raise ValueError(f"DNA clustering requested but {fasta_path} contains protein sequences")
            return fasta_path, "cd-hit-est", True
        if is_dna:
            aa_path = os.path.join(target, prefix + "_aa.fasta")
            translate_dna_to_aa(fasta_path, aa_path)
            return aa_path, "cd-hit", True
        return fasta_path, "cd-hit", False

    original_fasta, original_annotations = extract_input(
        "original", options.original_input_type, options.original_input_dir,
        options.original_input_fasta, options.original_name_split_gff,
        options.original_name_split_fasta, original_dir,
    )
    supplement_fasta, supplement_annotations = extract_input(
        "supplement", options.supplement_input_type, options.supplement_input_dir,
        options.supplement_input_fasta, options.supplement_name_split_gff,
        options.supplement_name_split_fasta, supplemented_dir,
    )
    original_cluster_input, clustering_mode, original_is_dna = clustering_input(
        original_fasta, original_dir, "original_clustering_sequences"
    )
    if options.sequence_type == "AA" and not original_is_dna and detect_sequence_type(supplement_fasta):
        # The incremental preparer translates DNA supplements when the source
        # collection is proteins, preserving one alphabet for reconciliation.
        pass
    elif options.sequence_type == "DNA" and not detect_sequence_type(supplement_fasta):
        raise ValueError(f"DNA clustering requested but {supplement_fasta} contains protein sequences")

    old_prefix = os.path.join(original_dir, "clustering_CD-HIT_original")
    run_cd_hit(options, original_cluster_input, old_prefix, clustering_mode)
    original_clusters = old_prefix + ".clstr"
    if not os.path.isfile(original_clusters):
        raise RuntimeError(f"CD-HIT did not create the original cluster file: {original_clusters}")

    def annotation_list(*paths):
        return [os.path.abspath(path) for path in paths if path]

    def analysis_options(run_mode, clusters, fasta, target, reclustered=None,
                         sequence_db_out=None, gene_annotation_files=None):
        return SimpleNamespace(
            run_mode=run_mode, cluster_format="CD-HIT", sequence_type=options.sequence_type,
            reclustered=reclustered, sequence_tag=options.sequence_tag,
            species_groups=groups_to_use, core_groups=groups_to_use, clusters=clusters,
            output_dir=target, gene_presence_absence_out=options.gene_presence_absence_out,
            write_groups=options.write_groups, write_individual_groups=options.write_individual_groups,
            threads=options.threads, mem=options.mem,
            sequence_buffer_cap_mb=options.sequence_buffer_cap_mb, tmp_dir=options.tmp_dir,
            align_core=options.align_core, align_aa=options.align_aa, fasta=fasta,
            sequence_db=None, sequence_db_out=sequence_db_out, verbose=options.verbose,
            single_copy_only=options.single_copy_only, single_copy_tolerance=options.single_copy_tolerance,
            gene_group_table=options.gene_group_table, core_group_table=options.core_group_table,
            visualise=options.visualise, visualise_heatmap_limit=options.visualise_heatmap_limit,
            visualise_permutations=options.visualise_permutations, visualise_seed=options.visualise_seed,
            taxon_map=options.taxon_map, sequence_delimiter=options.sequence_delimiter,
            blast_graph_mode=options.blast_graph_mode, blast_min_identity=options.blast_min_identity,
            blast_max_evalue=options.blast_max_evalue,
            blast_min_query_coverage=options.blast_min_query_coverage,
            blast_min_subject_coverage=options.blast_min_subject_coverage,
            blast_columns=options.blast_columns, added_genomes=None,
            gene_annotation_files=annotation_list(*(gene_annotation_files or [])),
        )

    original_options = analysis_options(
        "Full", original_clusters, original_fasta, original_dir,
        gene_annotation_files=[original_annotations],
    )
    logger.info("Analysing original genome set")
    species_cluster(original_options)

    prepared = prepare_incremental_inputs(
        cluster_file=original_clusters,
        original_fasta=original_fasta,
        new_genes_fasta=supplement_fasta,
        output_dir=supplemented_dir,
        sequence_delimiter=options.sequence_delimiter,
        logger=logger,
        generated_ids_use_pipe=False,
    )
    if prepared["original_is_dna"]:
        recluster_input = prepared["representatives_dna_fasta"]
        recluster_mode = "cd-hit-est" if options.sequence_type == "DNA" else "cd-hit"
        if options.sequence_type == "AA":
            recluster_input = os.path.join(supplemented_dir, "reclustering_input_aa.fasta")
            translate_dna_to_aa(prepared["representatives_dna_fasta"], recluster_input)
    else:
        recluster_input = prepared["representatives_fasta"]
        recluster_mode = "cd-hit"

    new_prefix = os.path.join(supplemented_dir, "clustering_CD-HIT_supplemented")
    run_cd_hit(options, recluster_input, new_prefix, recluster_mode)
    new_clusters = new_prefix + ".clstr"
    if not os.path.isfile(new_clusters):
        raise RuntimeError(f"CD-HIT did not create the supplemented cluster file: {new_clusters}")

    options.added_genomes = ",".join(prepared["new_genomes"])
    combined_annotations = [path for path in (original_annotations, supplement_annotations) if path]
    final_options = analysis_options(
        "Partial", original_clusters, prepared["combined_fasta"], supplemented_dir,
        reclustered=new_clusters,
        sequence_db_out=getattr(options, "sequence_db_out", None),
        gene_annotation_files=combined_annotations,
    )
    final_options.added_genomes = options.added_genomes
    final_options.supplemental_sequence_ids = prepared["new_gene_ids"]
    final_options.baseline_output_dir = original_dir

    logger.info(
        "Tracking %s supplementary sequence IDs for provenance",
        len(final_options.supplemental_sequence_ids),
    )
    logger.info("Analysing original representatives plus supplemental genomes")
    species_cluster(final_options)
    logger.info(
        "Supplementation complete: original outputs in %s; supplemented outputs in %s",
        original_dir, supplemented_dir,
    )


def main():
    # Initial console-only logger so welcome and parser.description are logged before argparse outputs.
    early_logger = configure_logger("PyamilySeq.PyamilySeq", enable_file=False, log_dir=None, verbose=False)
    # Use LoggingArgumentParser so usage/errors are emitted via logger
    parser = LoggingArgumentParser(logger_name="PyamilySeq.PyamilySeq")  # , description="PyamilySeq entrypoint")

    # Add subparsers for Full, Partial, Incremental, and Supplementation modes.
    subparsers = parser.add_subparsers(
        dest="run_mode", required=True,
        help="Choose Full, Partial, Incremental, or Supplementation mode."
    )

    # Full Mode Subparser
    full_parser = subparsers.add_parser("Full",
                                        help="Full mode: PyamilySeq to cluster with CD-HIT and process output.")
    # Full mode currently supports CD-HIT only. Keep the default for backward
    # compatibility while making the runtime check below reachable.
    full_parser.add_argument("-clustering_format", choices=['CD-HIT'], default='CD-HIT',
                             help="Clustering format used in Full mode (CD-HIT).")
    full_parser.add_argument("-output_dir", required=True,
                             help="Directory for all output files.")
    full_parser.add_argument("-input_type", choices=['separate', 'combined', 'fasta'], required=True,
                             help="Type of input files: 'separate' for matching FASTA and GFF files, 'combined' for GFF+FASTA, or 'fasta' for a prepared FASTA file.")
    full_parser.add_argument("-input_dir", required=False,
                             help="Directory containing GFF/FASTA files - Use with -input_type separate/combined.")
    full_parser.add_argument("-input_fasta", required=False,
                             help="Input FASTA file - Use with - input_type fasta.")
    full_parser.add_argument("-name_split_gff", required=False,
                             help="Substring to split filenames and extract genome names for gff files (e.g., '_combined.gff3') - Use with -input_type separate/combined.")
    full_parser.add_argument("-name_split_fasta", required=False,
                             help="Substring to split filenames and extract genome names for fasta files if named differently to paired gff files (e.g., '_dna.fasta') - Use with -input_type separate/combined.")
    full_parser.add_argument("-sequence_type", choices=['AA', 'DNA'], default="AA", required=False,
                             help="Clustering mode: 'DNA' or 'AA'.")
    full_parser.add_argument("-gene_ident", default="CDS", required=False,
                             help="Gene identifiers to extract sequences (e.g., 'CDS, tRNA').")
    full_parser.add_argument("-gene_annotations", dest="gene_annotation_file", default=None,
                             help="Optional gene-annotation TSV for prepared FASTA input. GFF input generates annotations automatically.")
    full_parser.add_argument("-c", type=str, dest="pident", default="0.90", required=False,
                             help="Sequence identity threshold for clustering (default: 0.90) - CD-HIT parameter '-c'.")
    full_parser.add_argument("-s", type=str, dest="len_diff", default="0.80", required=False,
                             help="Length difference threshold for clustering (default: 0.80) - CD-HIT parameter '-s'.")

    full_parser.add_argument("-fast_mode", action="store_true",
                             help="Enable fast mode for CD-HIT (not recommended) - CD-HIT parameter '-g'.")

    # Partial Mode Subparser
    partial_parser = subparsers.add_parser("Partial", help="Partial mode: PyamilySeq to process pre-clustered data.")
    partial_parser.add_argument("-clustering_format", type=_normalise_clustering_format,
                                choices=['CD-HIT', 'MMseqs', 'BLAST'], required=True,
                                help="Clustering format used: CD-HIT, MMseqs2, or BLAST. ")
    partial_parser.add_argument("-cluster_file", required=True,
                                help="Cluster file containing pre-clustered groups from CD-HIT, MMseqs, BLAST etc.")
    partial_parser.add_argument("-original_fasta", required=False,
                                help="FASTA file used in pre-clustering (Provide sequences in DNA form). "
                                     "Not required if -sequence_db is given.")
    partial_parser.add_argument("-sequence_db", required=False,
                                help="Checkpoint: path to a pyamilyseq_sequences_*.sqlite3 file already built "
                                     "created by a previous run using --sequence-db-out. If given, the "
                                     "original FASTA is NOT re-parsed - the existing database is reused directly.")
    partial_parser.add_argument("-output_dir", required=True,
                                help="Directory for all output files.")
    partial_parser.add_argument("-reclustered", required=False,
                                help="Clustering output file from a second round of clustering.")
    partial_parser.add_argument("-seq_tag", default="StORF", dest="sequence_tag", required=False,
                                help="Legacy compatibility label; reclustering provenance is inferred from first-round membership.")
    partial_parser.add_argument("-gene_annotations", dest="gene_annotation_file", default=None,
                                help="Optional TSV sidecar generated during GFF extraction, with sequence_id, gene_name and annotation columns.")

    # Incremental mode: reuse an existing CD-HIT partition and recluster only
    # the original representatives plus genes from the new genomes.
    incremental_parser = subparsers.add_parser(
        "Incremental",
        help="Incremental mode: recluster prior CD-HIT representatives with newly added genes, then reconcile all genes."
    )
    incremental_parser.add_argument("-clustering_format", choices=['CD-HIT'], default='CD-HIT',
                                    help="The original cluster file must be native CD-HIT .clstr format.")
    incremental_parser.add_argument("-cluster_file", required=True,
                                    help="Original full CD-HIT .clstr file from the existing pangenome run.")
    incremental_parser.add_argument("-original_fasta", required=True,
                                    help="Full nucleotide or protein FASTA used for the original clustering; its sequence alphabet is auto-detected and it must contain every old gene.")
    incremental_parser.add_argument("-output_dir", required=True,
                                    help="Directory for incremental intermediate files and final PyamilySeq outputs.")
    incremental_parser.add_argument("-new_input_type", choices=['combined', 'separate', 'fasta'], required=True,
                                    help="New gene input: combined GFF+FASTA, paired GFF/FASTA files, or a prepared nucleotide genes FASTA.")
    incremental_parser.add_argument("-new_input_dir", required=False,
                                    help="Directory containing only the new genomes' GFF/FASTA inputs.")
    incremental_parser.add_argument("-new_input_fasta", required=False,
                                    help="Prepared nucleotide genes FASTA for only the new genomes (IDs must be genome|gene).")
    incremental_parser.add_argument("-name_split_gff", required=False,
                                    help="Filename suffix for new GFF files, e.g. '_combined.gff3'.")
    incremental_parser.add_argument("-name_split_fasta", required=False,
                                    help="Optional FASTA filename suffix for paired inputs in separate mode.")
    incremental_parser.add_argument("-sequence_type", choices=['AA', 'DNA'], default='AA',
                                    help="Clustering sequence type used by the original run; keep it and -c/-s consistent with that run.")
    incremental_parser.add_argument("-gene_ident", default='CDS',
                                    help="GFF feature type(s) to extract from new genomes, e.g. 'CDS' or 'CDS,tRNA'.")
    incremental_parser.add_argument("-c", type=str, dest="pident", default='0.90',
                                    help="CD-HIT identity threshold; use the same value as the original run (default: 0.90).")
    incremental_parser.add_argument("-s", type=str, dest="len_diff", default='0.80',
                                    help="CD-HIT length threshold; use the same value as the original run (default: 0.80).")
    incremental_parser.add_argument("-fast_mode", action='store_true',
                                    help="Enable CD-HIT fast mode (not recommended; match the original run settings).")
    incremental_parser.add_argument("-seq_tag", default='StORF', dest='sequence_tag', required=False,
                                    help="Legacy compatibility label; family provenance is inferred from source membership.")
    incremental_parser.add_argument("-gene_annotations", dest="gene_annotation_file", default=None,
                                    help="Optional annotation TSV/TSV.GZ for the original FASTA; auto-detected beside original_fasta as gene_annotations.tsv.gz or gene_annotations.tsv when present.")
    incremental_parser.add_argument("-new_gene_annotations", dest="new_gene_annotation_file", default=None,
                                    help="Optional annotation TSV when the new input is a prepared FASTA. GFF inputs generate it automatically.")

    supplementation_parser = subparsers.add_parser(
        "Supplementation",
        help="Cluster and analyse an original genome set, then recluster its representatives with supplemental genomes."
    )
    supplementation_parser.add_argument(
        "-clustering_format", choices=["CD-HIT"], default="CD-HIT",
        help="Clustering format used by the supplementation workflow (CD-HIT)."
    )
    supplementation_parser.add_argument("-output_dir", required=True, help="Root directory for both analyses.")
    for prefix, description in (("original", "original"), ("supplement", "supplemental")):
        supplementation_parser.add_argument(
            f"-{prefix}_input_type", choices=["combined", "separate", "fasta"], required=True,
            help=f"{description.capitalize()} input type: combined GFF+FASTA, paired GFF/FASTA, or extracted FASTA."
        )
        supplementation_parser.add_argument(
            f"-{prefix}_input_dir", help=f"Directory containing {description} GFF or paired GFF/FASTA files."
        )
        supplementation_parser.add_argument(
            f"-{prefix}_input_fasta", help=f"Extracted FASTA for {description} genomes."
        )
        supplementation_parser.add_argument(
            f"-{prefix}_name_split_gff", help=f"GFF suffix for {description} input directory."
        )
        supplementation_parser.add_argument(
            f"-{prefix}_name_split_fasta", help=f"Optional FASTA suffix for paired {description} inputs."
        )
        supplementation_parser.add_argument(
            f"-{prefix}_gene_annotations",
            help=f"Optional existing gene-annotation TSV when {description} input is a prepared FASTA. "
                 "GFF inputs generate this metadata automatically.",
        )
    supplementation_parser.add_argument("-sequence_type", choices=["AA", "DNA"], default="AA")
    supplementation_parser.add_argument("-gene_ident", default="CDS")
    supplementation_parser.add_argument("-c", type=str, dest="pident", default="0.90")
    supplementation_parser.add_argument("-s", type=str, dest="len_diff", default="0.80")
    supplementation_parser.add_argument("-fast_mode", action="store_true")
    supplementation_parser.add_argument("-seq_tag", default="StORF", dest="sequence_tag")

    # Common Grouping Arguments
    for subparser in [full_parser, partial_parser, incremental_parser, supplementation_parser]:
        subparser.add_argument("-group_mode", choices=['Species', 'Genus'], default="Species", required=False,
                               help="Grouping mode: 'Species' or 'Genus'.")
        subparser.add_argument("-species_groups", default="99,95,15", required=False,
                               help="Gene groupings for 'Species' mode (default: '99,95,15').")
        subparser.add_argument("-genus_groups", default="1,2,3,4,5,6,7,8,9,10", required=False,
                               help="Gene groupings for 'Genus' mode (default: '1-10').")
        subparser.add_argument("-write_groups", default=None, dest="write_groups", required=False,
                               help="Output gene groups as a single FASTA file (e.g., '99,95'). Triggers writing individual groups.")
        subparser.add_argument("-write_individual_groups", action="store_true", dest="write_individual_groups",
                               help="Output individual FASTA files for each group.")
        subparser.add_argument("-align", action="store_true", dest="align_core",
                               help="Align and concatenate sequences for 'core' groups (those in 99-100%% of genomes).")
        subparser.add_argument("-align_aa", action="store_true",
                               help="Align sequences as amino acids.")
        subparser.add_argument("-no_gpa", action="store_false", dest="gene_presence_absence_out",
                               help="Skip creation of gene_presence_absence.csv.")
        subparser.add_argument("-M", type=int, default=4000, dest="mem", required=False,
                               help="Memory budget in MB. Used for FASTA/cluster processing; in Full, Incremental, and Supplementation modes it is passed to CD-HIT.")
        subparser.add_argument("--sequence-buffer-cap-mb", type=int, default=256, dest="sequence_buffer_cap_mb",
                               help="Maximum chunk size in MB for the FASTA-backed sequence store used during loading.")
        subparser.add_argument("--tmp-dir", dest="tmp_dir", default=None,
                               help="Directory for temporary files such as chunked sequence stores and MAFFT scratch files.")
        subparser.add_argument(
            "--sequence-db-out",
            dest="sequence_db_out",
            default=None,
            help=("Create and retain a validated SQLite sequence-store checkpoint at this path. "
                  "In Supplementation mode this retains the final combined database only. "
                  "The path must not already exist."),
        )
        subparser.add_argument("-T", type=int, default=8, dest="threads", required=False,
                               help="Number of threads for clustering/alignment - CD-HIT parameter '-T' | MAFFT parameter '--thread'.")
        # Options for selecting single-copy families
        subparser.add_argument("--single_copy_only", action="store_true", dest="single_copy_only",
                               help="If set, only report gene families that are strictly single-copy across genomes.")
        subparser.add_argument("--single_copy_tolerance", type=float, dest="single_copy_tolerance", default=0.0,
                               help=("Allow up to this percentage (0-100) of genomes in a family to have multi-copies. "
                                     "If >0, families with multi-copy genomes less than or equal to this percentage are kept."))
        subparser.add_argument("-gene_group_table", action="store_true", dest="gene_group_table",
                               help="Output a tab-separated gene-to-family map using stable hashed family IDs.")
        subparser.add_argument("-core_group_table", action="store_true", dest="core_group_table",
                               help="Output a tab-separated gene-to-prevalence-category map.")
        subparser.add_argument("--added-genomes", dest="added_genomes", default=None,
                               help="Comma-separated genome IDs added in this reclustering run; write extension impact tables and chart.")
        subparser.add_argument("-visualise", "--visualize", action="store_true", dest="visualise",
                               help="Write SVG/HTML visual summaries of the pangenome outputs.")
        subparser.add_argument("-visualise_heatmap_limit", "--visualize-heatmap-limit", type=int, default=200,
                               dest="visualise_heatmap_limit",
                               help="Maximum number of gene families to show in the presence/absence heatmap (default: 200).")
        subparser.add_argument("--visualize-permutations", type=int, default=100, dest="visualise_permutations",
                               help="Permutation count for deterministic pangenome accumulation curves (default: 100).")
        subparser.add_argument("--visualize-seed", type=int, default=1, dest="visualise_seed",
                               help="Random seed for pangenome accumulation curves (default: 1).")
        subparser.add_argument("--taxon-map",
                               help="Optional two/three-column identifier-to-taxon mapping, recommended for Genus mode.")
        subparser.add_argument("--sequence-delimiter", default="|",
                               help="Delimiter between genome and gene in sequence IDs (default: '|').")
        subparser.add_argument("--blast-graph-mode", choices=["reciprocal-best", "connected-components"],
                               default="reciprocal-best",
                               help="How scored BLAST hits are converted into family edges (default: reciprocal-best).")
        subparser.add_argument("--blast-min-identity", type=float, default=0.0,
                               help="Minimum BLAST percent identity (default: 0).")
        subparser.add_argument("--blast-max-evalue", type=float, default=10.0,
                               help="Maximum BLAST e-value (default: 10).")
        subparser.add_argument("--blast-min-query-coverage", type=float, default=0.0,
                               help="Minimum query coverage percentage; requires qlen in --blast-columns.")
        subparser.add_argument("--blast-min-subject-coverage", type=float, default=0.0,
                               help="Minimum subject coverage percentage; requires slen in --blast-columns.")
        subparser.add_argument(
            "--blast-columns",
            default="qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore",
            help="Comma-separated BLAST column names matching the input file.",
        )
        # Mirror global options so argparse accepts them before or after the subcommand.
        subparser.add_argument("--log", action="store_true", dest="log", default=argparse.SUPPRESS,
                               help="Create a timestamped logfile for this run.")
        subparser.add_argument("--log-dir", dest="log_dir", default=argparse.SUPPRESS,
                               help="Directory for logfile (default: output dir or cwd).")
        subparser.add_argument("-verbose", action="store_true", default=argparse.SUPPRESS,
                               help="Print verbose output.")

    analysis_parser = subparsers.add_parser("Analyze",
                                            help="Run genome-order, gene-neighborhood, or cohort-scale analyses.")
    analysis_parser.add_argument("analysis_tool", choices=("order", "location", "cohort"),
                                 help="Analysis workflow to run.")
    analysis_parser.add_argument("analysis_args", nargs=argparse.REMAINDER,
                                 help="Arguments forwarded to the selected analysis workflow.")

    # Miscellaneous Arguments
    # Global logging options (user controls logfile creation)
    parser.add_argument("--log", action="store_true", dest="log", help="Create a timestamped logfile for this run.")
    parser.add_argument("--log-dir", dest="log_dir", default=None,
                        help="Directory for logfile (default: output dir or cwd).")
    parser.add_argument("-verbose", action="store_true",
                        help="Print verbose output.")
    parser.add_argument("-v", "--version", action="version",
                        version=f"PyamilySeq {PyamilySeq_Version}: Exiting.")

    # Parse Arguments
    options = parser.parse_args()

    if options.run_mode == "Analyze":
        if options.analysis_tool == "order":
            try:
                from .order_analysis import main as analysis_main
            except ImportError:
                from order_analysis import main as analysis_main
        elif options.analysis_tool == "location":
            try:
                from .location_tracking import run_from_cli as analysis_main
            except ImportError:
                from location_tracking import run_from_cli as analysis_main
        else:
            try:
                from .cohort_analysis import main as analysis_main
            except ImportError:
                from cohort_analysis import main as analysis_main
        return analysis_main(options.analysis_args)

    # Setup logger once we know output paths/options
    # after we resolve output_path / options.output_dir:
    resolved_log_dir = options.log_dir if getattr(options, "log_dir", None) else (
        os.path.abspath(options.output_dir) if getattr(options, "output_dir", None) else os.getcwd())
    logger = configure_logger("PyamilySeq.PyamilySeq", enable_file=getattr(options, "log", False),
                              log_dir=resolved_log_dir, verbose=options.verbose)
    logger.info("Running PyamilySeq %s in %s mode", PyamilySeq_Version, getattr(options, "run_mode", "N/A"))
    if options.verbose:
        logger.debug("Options: %s", vars(options))

    if options.run_mode == 'Full':
        options.clustering_format = 'CD-HIT'
        if getattr(options, 'reclustered', None) is not None:
            sys.exit("Currently reclustering only works on Partial Mode.")
        required_full_mode = [options.input_type, options.pident, options.len_diff]
        if options.input_type != 'fasta':
            required_full_mode.extend([options.input_dir, options.name_split_gff])
        else:
            required_full_mode.append(options.input_fasta)
        if all(required_full_mode):
            # Proceed with the Full mode
            pass
        else:
            required_names = ['input_type', 'pident', 'len_diff']
            if options.input_type == 'fasta':
                required_names.append('input_fasta')
            else:
                required_names.extend(['input_dir', 'name_split_gff'])
            missing_options = [opt for opt in required_names if not options.__dict__.get(opt)]
            logger.error("Missing required options for Full mode: %s", ', '.join(missing_options))
            sys.exit(f"Missing required options for Full mode: {', '.join(missing_options)}")
        if options.align_core:
            options.write_individual_groups = True
            if options.write_groups == None:
                sys.exit('Must provide "-write_groups" to output gene groups before alignment "-align" can be done.')
    elif options.run_mode == 'Partial':
        required_partial_mode = [options.cluster_file, options.original_fasta]
        if all(required_partial_mode):
            pass
        else:
            missing_options = [opt for opt in ['cluster_file', 'original_fasta']
                               if not options.__dict__.get(opt)]
            sys.exit(f"Missing required options for Partial mode: {', '.join(missing_options)}")
        if options.align_core:
            options.write_individual_groups = True
            if options.write_groups is None or options.original_fasta is None:
                sys.exit('Must provide "-write_groups" to output gene groups before alignment "-align" can be done.')
    elif options.run_mode == 'Incremental':
        missing_options = [
            name for name in ('cluster_file', 'original_fasta', 'output_dir', 'new_input_type')
            if not options.__dict__.get(name)
        ]
        if options.new_input_type == 'fasta':
            if not options.new_input_fasta:
                missing_options.append('new_input_fasta')
        else:
            if not options.new_input_dir:
                missing_options.append('new_input_dir')
            if not options.name_split_gff:
                missing_options.append('name_split_gff')
        if missing_options:
            sys.exit(f"Missing required options for Incremental mode: {', '.join(missing_options)}")
        if options.align_core:
            options.write_individual_groups = True
            if options.write_groups is None:
                sys.exit('Must provide "-write_groups" to output gene groups before alignment "-align" can be done.')
    elif options.run_mode == 'Supplementation':
        missing = []
        for prefix in ('original', 'supplement'):
            input_type = getattr(options, f"{prefix}_input_type")
            if input_type == 'fasta':
                if not getattr(options, f"{prefix}_input_fasta"):
                    missing.append(f"{prefix}_input_fasta")
            else:
                for suffix in ('input_dir', 'name_split_gff'):
                    if not getattr(options, f"{prefix}_{suffix}"):
                        missing.append(f"{prefix}_{suffix}")
        if missing:
            sys.exit(f"Missing required options for Supplementation mode: {', '.join(missing)}")
        if options.align_core:
            options.write_individual_groups = True
            if options.write_groups is None:
                sys.exit('Must provide "-write_groups" to output gene groups before alignment "-align" can be done.')

    if getattr(options, "clustering_format", "").upper() == "CD-HIT":
        clust_affix = '.clstr'
    if getattr(options, "clustering_format", "").upper() == "TSV":
        clust_affix = '.tsv'
    if getattr(options, "clustering_format", "").upper() == "CSV":
        clust_affix = '.csv'

    ###External tool checks:
    ##MAFFT
    if options.align_core == True:
        if is_tool_installed('mafft'):
            logger.info("mafft is installed. Proceeding with alignment.")
        else:
            logger.error("mafft is not installed. Please install mafft to proceed.")
            sys.exit("mafft is not installed. Please install mafft to proceed.")
    ##CD-HIT
    if options.run_mode in ('Full', 'Incremental', 'Supplementation'):
        clustering_mode = 'cd-hit-est' if options.sequence_type == 'DNA' else 'cd-hit'
        if is_tool_installed(clustering_mode):
            logger.info("%s executable found.", clustering_mode)
            if options.fast_mode is True:
                options.fast_mode = 0
                logger.info("CD-HIT configured to run in fast mode.")
            else:
                options.fast_mode = 1
                logger.info("CD-HIT configured to run in accurate mode.")
        else:
            logger.error("%s is not installed. Please install CD-HIT to proceed.", clustering_mode)
            sys.exit(f"{clustering_mode} is not installed. Please install CD-HIT to proceed.")

    if hasattr(options, 'cluster_file') and options.cluster_file:
        options.cluster_file = fix_path(options.cluster_file)
    if hasattr(options, 'reclustered') and options.reclustered:
        options.reclustered = fix_path(options.reclustered)
    if getattr(options, 'input_dir', None):
        options.input_dir = fix_path(options.input_dir)
    if getattr(options, 'input_fasta', None):
        options.input_fasta = fix_path(options.input_fasta)
    if getattr(options, 'new_input_dir', None):
        options.new_input_dir = fix_path(options.new_input_dir)
    if getattr(options, 'new_input_fasta', None):
        options.new_input_fasta = fix_path(options.new_input_fasta)
    if getattr(options, 'original_fasta', None):
        options.original_fasta = fix_path(options.original_fasta)
    for prefix in ('original', 'supplement'):
        for suffix in ('input_dir', 'input_fasta'):
            attr = f"{prefix}_{suffix}"
            if getattr(options, attr, None):
                setattr(options, attr, fix_path(getattr(options, attr)))
    if hasattr(options, 'output_dir') and options.output_dir:
        options.output_dir = fix_path(options.output_dir)
    if getattr(options, 'taxon_map', None):
        options.taxon_map = fix_path(options.taxon_map)
    for annotation_attribute in (
        "gene_annotation_file", "new_gene_annotation_file",
        "original_gene_annotations", "supplement_gene_annotations",
    ):
        annotation_path = getattr(options, annotation_attribute, None)
        if annotation_path:
            setattr(options, annotation_attribute, fix_path(annotation_path))
    if getattr(options, 'sequence_buffer_cap_mb', None) is not None and options.sequence_buffer_cap_mb <= 0:
        sys.exit("--sequence-buffer-cap-mb must be greater than 0")
    if getattr(options, 'tmp_dir', None):
        options.tmp_dir = fix_path(options.tmp_dir)
        os.makedirs(options.tmp_dir, exist_ok=True)

    if getattr(options, 'sequence_db_out', None):
        options.sequence_db_out = os.path.abspath(fix_path(options.sequence_db_out))
        db_parent = os.path.dirname(options.sequence_db_out) or os.getcwd()
        os.makedirs(db_parent, exist_ok=True)
        if any(os.path.exists(candidate) for candidate in (
            options.sequence_db_out,
            options.sequence_db_out + "-wal",
            options.sequence_db_out + "-shm",
            options.sequence_db_out + "-journal",
        )):
            sys.exit(
                f"Sequence database output already exists: {options.sequence_db_out}. "
                "Choose a new path; existing checkpoints are never overwritten."
            )

    if options.run_mode == 'Partial':
        if getattr(options, 'sequence_db', None):
            options.sequence_db = fix_path(options.sequence_db)
            if getattr(options, 'sequence_db_out', None):
                sys.exit("Partial mode cannot combine -sequence_db with --sequence-db-out; reuse the existing database or build a new one from -original_fasta.")
        elif not getattr(options, 'original_fasta', None):
            sys.exit("Partial mode requires either -original_fasta or -sequence_db.")
        if getattr(options, 'original_fasta', None):
            options.original_fasta = fix_path(options.original_fasta)

    output_path = os.path.abspath(options.output_dir)
    combined_out_file = os.path.join(output_path, "input_genes_dna.fasta")
    clustering_output = os.path.join(output_path, 'clustering_' + options.clustering_format)

    if options.group_mode == 'Species':
        options.species_groups = options.species_groups + ',0'
        groups_to_use = options.species_groups
    elif options.group_mode == 'Genus':
        options.genus_groups = options.genus_groups + ',>'
        groups_to_use = options.genus_groups
        if options.align_core:
            sys.exit("-a align_core not a valid option in Genus mode.")

    if options.run_mode == 'Full':
        if getattr(options, "clustering_format", "").upper() != "CD-HIT":
            sys.exit('Only CD-HIT clustering works in Full Mode')

        if not os.path.exists(output_path):
            os.makedirs(output_path)
        if options.sequence_type == 'AA':
            clustering_mode = 'cd-hit'
            file_to_cluster = combined_out_file.replace('_dna.fasta', '_aa.fasta')
            translate = True
        elif options.sequence_type == 'DNA':
            clustering_mode = 'cd-hit-est'
            translate = False
            file_to_cluster = combined_out_file
        gene_annotation_files = []
        if options.input_type == 'separate':
            annotation_path = os.path.join(output_path, "gene_annotations.tsv.gz")
            read_separate_files(
                options.input_dir, options.name_split_gff, options.name_split_fasta,
                options.gene_ident, combined_out_file, translate, False,
                annotation_out=annotation_path,
            )
            gene_annotation_files = [annotation_path]
            run_cd_hit(options, file_to_cluster, clustering_output, clustering_mode)
        elif options.input_type == 'combined':
            annotation_path = os.path.join(output_path, "gene_annotations.tsv.gz")
            read_combined_files(
                options.input_dir, options.name_split_gff, options.gene_ident,
                combined_out_file, translate, False,
                progress_logger=logger, annotation_out=annotation_path,
            )
            gene_annotation_files = [annotation_path]
            run_cd_hit(options, file_to_cluster, clustering_output, clustering_mode)
        elif options.input_type == 'fasta':
            combined_out_file = options.input_fasta
            # Detect if the input FASTA file contains DNA or AA sequences
            is_dna = detect_sequence_type(options.input_fasta)
            # If the sequence type is AA and the input is DNA, translate the DNA to AA
            if options.sequence_type == 'AA' and is_dna:
                translated_fasta = os.path.join(output_path, os.path.splitext(os.path.basename(options.input_fasta))[
                    0] + '_aa.fasta')
                translate_dna_to_aa(options.input_fasta, translated_fasta)
                file_to_cluster = translated_fasta
            else:
                file_to_cluster = options.input_fasta
            if options.input_type == 'fasta':
                annotation_path = getattr(options, "gene_annotation_file", None)
                if annotation_path is None:
                    annotation_path = _gene_annotation_sidecar_for_fasta(options.input_fasta)
                if annotation_path:
                    gene_annotation_files = [annotation_path]
            run_cd_hit(options, file_to_cluster, clustering_output, clustering_mode)

        class clustering_options:
            def __init__(self):
                self.run_mode = options.run_mode
                self.cluster_format = options.clustering_format
                self.sequence_type = options.sequence_type
                self.reclustered = None
                self.sequence_tag = None
                self.species_groups = groups_to_use
                self.core_groups = groups_to_use
                self.clusters = clustering_output + clust_affix
                self.output_dir = options.output_dir
                self.gene_presence_absence_out = options.gene_presence_absence_out
                self.write_groups = options.write_groups
                self.write_individual_groups = options.write_individual_groups
                self.threads = options.threads
                self.mem = options.mem
                self.sequence_buffer_cap_mb = options.sequence_buffer_cap_mb
                self.tmp_dir = options.tmp_dir
                self.align_core = options.align_core
                self.align_aa = options.align_aa
                self.fasta = combined_out_file
                self.sequence_db = None
                self.sequence_db_out = getattr(options, "sequence_db_out", None)
                self.verbose = options.verbose
                # single-copy reporting options
                self.single_copy_only = options.single_copy_only
                self.single_copy_tolerance = options.single_copy_tolerance
                self.gene_group_table = options.gene_group_table
                self.core_group_table = options.core_group_table
                self.visualise = options.visualise
                self.visualise_heatmap_limit = options.visualise_heatmap_limit
                self.visualise_permutations = options.visualise_permutations
                self.visualise_seed = options.visualise_seed
                self.taxon_map = options.taxon_map
                self.sequence_delimiter = options.sequence_delimiter
                self.blast_graph_mode = options.blast_graph_mode
                self.blast_min_identity = options.blast_min_identity
                self.blast_max_evalue = options.blast_max_evalue
                self.blast_min_query_coverage = options.blast_min_query_coverage
                self.blast_min_subject_coverage = options.blast_min_subject_coverage
                self.blast_columns = options.blast_columns
                self.added_genomes = options.added_genomes
                self.gene_annotation_files = gene_annotation_files

        clustering_options = clustering_options()

    elif options.run_mode == 'Incremental':
        os.makedirs(output_path, exist_ok=True)
        if getattr(options, "clustering_format", "").upper() != "CD-HIT":
            sys.exit("Incremental mode currently requires a native CD-HIT .clstr file.")

        # Extract only the newly added gene sequences. Existing genes are copied
        # from the authoritative full DNA FASTA, not re-extracted from old GFFs.
        new_genes_dna = os.path.join(output_path, "supplemental_genes_dna.fasta")
        new_annotation_path = None
        if options.new_input_type == 'combined':
            new_annotation_path = os.path.join(output_path, "supplemental_gene_annotations.tsv.gz")
            read_combined_files(
                options.new_input_dir, options.name_split_gff, options.gene_ident,
                new_genes_dna, False, False, annotation_out=new_annotation_path,
            )
            new_genes_source = new_genes_dna
        elif options.new_input_type == 'separate':
            new_annotation_path = os.path.join(output_path, "supplemental_gene_annotations.tsv.gz")
            read_separate_files(
                options.new_input_dir, options.name_split_gff, options.name_split_fasta,
                options.gene_ident, new_genes_dna, False, False,
                annotation_out=new_annotation_path,
            )
            new_genes_source = new_genes_dna
        else:
            new_genes_source = options.new_input_fasta
            new_annotation_path = getattr(options, "new_gene_annotation_file", None)
            if new_annotation_path is None:
                new_annotation_path = _gene_annotation_sidecar_for_fasta(new_genes_source)

        original_annotation_path = getattr(options, "gene_annotation_file", None)
        if original_annotation_path is None:
            original_annotation_path = _gene_annotation_sidecar_for_fasta(options.original_fasta)
        gene_annotation_files = [
            path for path in (original_annotation_path, new_annotation_path) if path
        ]

        try:
            prepared = prepare_incremental_inputs(
                cluster_file=options.cluster_file,
                original_fasta=options.original_fasta,
                new_genes_fasta=new_genes_source,
                output_dir=output_path,
                sequence_delimiter=options.sequence_delimiter,
                logger=logger,
                generated_ids_use_pipe=(options.new_input_type in ('combined', 'separate')),
            )
        except IncrementalInputError as error:
            logger.error("Cannot prepare incremental clustering inputs: %s", error)
            sys.exit(2)

        # Match the source alphabet used by the original clustering. When the
        # original FASTA is protein, preparation has translated only the new
        # nucleotide genes, so the representative input is ready for CD-HIT.
        if options.sequence_type == 'AA':
            if prepared["original_is_dna"]:
                file_to_cluster = os.path.join(output_path, "reclustering_input_aa.fasta")
                translate_dna_to_aa(prepared["representatives_dna_fasta"], file_to_cluster)
            else:
                file_to_cluster = prepared["representatives_fasta"]
            clustering_mode = 'cd-hit'
        else:
            file_to_cluster = prepared["representatives_dna_fasta"]
            clustering_mode = 'cd-hit-est'

        recluster_prefix = os.path.join(output_path, "clustering_CD-HIT_incremental")
        run_cd_hit(options, file_to_cluster, recluster_prefix, clustering_mode)
        reclustered_file = recluster_prefix + ".clstr"
        if not os.path.isfile(reclustered_file):
            logger.error("CD-HIT finished but expected reclustering output was not created: %s", reclustered_file)
            sys.exit(1)

        # Pass the original cluster file plus all old/new sequences and the new
        # representative-level clustering into the existing reconciliation path.
        # This is what restores original non-representative genes to their families.
        if options.added_genomes is None:
            options.added_genomes = ",".join(prepared["new_genomes"])

        class clustering_options:
            def __init__(self):
                self.run_mode = 'Partial'
                self.cluster_format = 'CD-HIT'
                self.sequence_type = None
                self.reclustered = reclustered_file
                self.sequence_tag = options.sequence_tag
                self.species_groups = groups_to_use
                self.core_groups = groups_to_use
                self.clusters = options.cluster_file
                self.output_dir = options.output_dir
                self.gene_presence_absence_out = options.gene_presence_absence_out
                self.write_groups = options.write_groups
                self.write_individual_groups = options.write_individual_groups
                self.threads = options.threads
                self.mem = options.mem
                self.sequence_buffer_cap_mb = options.sequence_buffer_cap_mb
                self.tmp_dir = options.tmp_dir
                self.align_core = options.align_core
                self.align_aa = options.align_aa
                self.fasta = prepared["combined_fasta"]
                self.sequence_db = None
                self.sequence_db_out = getattr(options, "sequence_db_out", None)
                self.verbose = options.verbose
                self.single_copy_only = options.single_copy_only
                self.single_copy_tolerance = options.single_copy_tolerance
                self.gene_group_table = options.gene_group_table
                self.core_group_table = options.core_group_table
                self.visualise = options.visualise
                self.visualise_heatmap_limit = options.visualise_heatmap_limit
                self.visualise_permutations = options.visualise_permutations
                self.visualise_seed = options.visualise_seed
                self.taxon_map = options.taxon_map
                self.sequence_delimiter = options.sequence_delimiter
                self.blast_graph_mode = options.blast_graph_mode
                self.blast_min_identity = options.blast_min_identity
                self.blast_max_evalue = options.blast_max_evalue
                self.blast_min_query_coverage = options.blast_min_query_coverage
                self.blast_min_subject_coverage = options.blast_min_subject_coverage
                self.blast_columns = options.blast_columns
                self.added_genomes = options.added_genomes
                self.gene_annotation_files = gene_annotation_files

        clustering_options = clustering_options()
        logger.info(
            "Incremental clustering prepared: %s original families, %s new genes across %s genomes; output %s",
            f"{prepared['original_cluster_count']:,}", f"{prepared['new_gene_count']:,}",
            f"{len(prepared['new_genomes']):,}", reclustered_file,
        )

    elif options.run_mode == 'Supplementation':
        if getattr(options, "clustering_format", None) not in (None, "CD-HIT"):
            sys.exit("Supplementation mode currently requires CD-HIT.")
        try:
            _run_supplementation(options, output_path, groups_to_use, logger)
        except UnicodeDecodeError as error:
            logger.exception(
                "Supplementation encountered a text decoding error. "
                "Encoding=%s, invalid bytes=%s, position=%s.",
                error.encoding,
                error.object[error.start:error.end].hex(" "),
                error.start,
            )
            sys.exit(2)

        except (ValueError, IncrementalInputError) as error:
            logger.error("Supplementation could not continue: %s", error)
            sys.exit(2)
        except Exception:
            logger.exception("Supplementation workflow failed")
            sys.exit(1)
        return

    elif options.run_mode == 'Partial':
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        class clustering_options:
            def __init__(self):
                self.run_mode = options.run_mode
                self.cluster_format = options.clustering_format
                self.sequence_type = None
                self.reclustered = options.reclustered
                self.sequence_tag = options.sequence_tag
                self.species_groups = groups_to_use
                self.core_groups = groups_to_use
                self.clusters = options.cluster_file
                self.output_dir = options.output_dir
                self.gene_presence_absence_out = options.gene_presence_absence_out
                self.write_groups = options.write_groups
                self.write_individual_groups = options.write_individual_groups
                self.threads = options.threads
                self.mem = options.mem
                self.sequence_buffer_cap_mb = options.sequence_buffer_cap_mb
                self.tmp_dir = options.tmp_dir
                self.align_core = options.align_core
                self.align_aa = options.align_aa
                self.fasta = options.original_fasta
                self.sequence_db = getattr(options, 'sequence_db', None)
                self.sequence_db_out = getattr(options, "sequence_db_out", None)
                self.verbose = options.verbose
                # single-copy reporting options
                self.single_copy_only = options.single_copy_only
                self.single_copy_tolerance = options.single_copy_tolerance
                self.gene_group_table = options.gene_group_table
                self.core_group_table = options.core_group_table
                self.visualise = options.visualise
                self.visualise_heatmap_limit = options.visualise_heatmap_limit
                self.visualise_permutations = options.visualise_permutations
                self.visualise_seed = options.visualise_seed
                self.taxon_map = options.taxon_map
                self.sequence_delimiter = options.sequence_delimiter
                self.blast_graph_mode = options.blast_graph_mode
                self.blast_min_identity = options.blast_min_identity
                self.blast_max_evalue = options.blast_max_evalue
                self.blast_min_query_coverage = options.blast_min_query_coverage
                self.blast_min_subject_coverage = options.blast_min_subject_coverage
                self.blast_columns = options.blast_columns
                self.added_genomes = options.added_genomes
                annotation_path = getattr(options, "gene_annotation_file", None)
                if annotation_path is None and getattr(options, "original_fasta", None):
                    annotation_path = _gene_annotation_sidecar_for_fasta(options.original_fasta)
                self.gene_annotation_files = [annotation_path] if annotation_path else []

        clustering_options = clustering_options()

    if options.group_mode == 'Species':
        try:
            logger.info("Dispatching species clustering")
            species_cluster(clustering_options)
            logger.info("Invoked species clustering.")
        except FileNotFoundError as e:
            logger.error("File not found during species clustering: %s", e)
            logger.debug("Traceback:\n%s", traceback.format_exc())
            sys.exit(1)
        except MemoryError:
            logger.exception(
                "Out of memory during pangenome processing. The reconciliation workflow may load substantial cluster metadata into RAM."
            )
            sys.exit(1)
        except Exception:
            logger.exception("Unexpected error during species clustering")
            sys.exit(1)
    elif options.group_mode == 'Genus':
        try:
            logger.info("Dispatching genus clustering")
            genus_cluster(clustering_options)
            logger.info("Invoked genus clustering.")
        except FileNotFoundError as e:
            logger.error("File not found during genus clustering: %s", e)
            logger.debug("Traceback:\n%s", traceback.format_exc())
            sys.exit(1)
        except MemoryError:
            logger.exception(
                "Out of memory during pangenome processing. The reconciliation workflow may load substantial cluster metadata into RAM."
            )
            sys.exit(1)
        except Exception:
            logger.exception("Unexpected error during genus clustering")
            sys.exit(1)

    # Save arguments to a text file
    from datetime import datetime
    with open(output_path + "/PyamilySeq_params.txt", "w") as outfile:
        outfile.write(f"Timestamp: {datetime.now().isoformat()}\n")
        for arg, value in vars(options).items():
            outfile.write(f"{arg}: {value}\n")
    logger.info("Saved parameters to %s", os.path.join(output_path, "PyamilySeq_params.txt"))


if __name__ == "__main__":
    main()