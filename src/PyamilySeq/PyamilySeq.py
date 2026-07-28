try:
    from .PyamilySeq_Species import cluster as species_cluster
    from .PyamilySeq_Genus import cluster as genus_cluster
    from .constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from PyamilySeq_Species import cluster as species_cluster
    from PyamilySeq_Genus import cluster as genus_cluster
    from constants import *
    from utils import *

import traceback
import sys
import argparse

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


def main():
    # Initial console-only logger so welcome and parser.description are logged before argparse outputs.
    early_logger = configure_logger("PyamilySeq.PyamilySeq", enable_file=False, log_dir=None, verbose=False)
    # Use LoggingArgumentParser so usage/errors are emitted via logger
    parser = LoggingArgumentParser(logger_name="PyamilySeq.PyamilySeq")#, description="PyamilySeq entrypoint")

    # Add subparsers for Full and Partial modes
    subparsers = parser.add_subparsers(dest="run_mode", required=True, help="Choose a mode: 'Full' or 'Partial'.")

    # Full Mode Subparser
    full_parser = subparsers.add_parser("Full",
                                        help="Full mode: PyamilySeq to cluster with CD-HIT and process output.")
    #full_parser.add_argument("-clustering_format", choices=['CD-HIT', 'MMseqs', 'BLAST'], required=True,
    #                         help="Clustering format to use: CD-HIT, MMseqs2, or BLAST.")
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
    full_parser.add_argument("-c", type=str, dest="pident", default="0.90", required=False,
                             help="Sequence identity threshold for clustering (default: 0.90) - CD-HIT parameter '-c'.")
    full_parser.add_argument("-s", type=str, dest="len_diff", default="0.80", required=False,
                             help="Length difference threshold for clustering (default: 0.80) - CD-HIT parameter '-s'.")

    full_parser.add_argument("-fast_mode", action="store_true",
                             help="Enable fast mode for CD-HIT (not recommended) - CD-HIT parameter '-g'.")


    # Partial Mode Subparser
    partial_parser = subparsers.add_parser("Partial", help="Partial mode: PyamilySeq to process pre-clustered data.")
    partial_parser.add_argument("-clustering_format", type=_normalise_clustering_format, choices=['CD-HIT', 'MMseqs', 'BLAST'], required=True,
                                help="Clustering format used: CD-HIT, MMseqs2, or BLAST. ")
    partial_parser.add_argument("-cluster_file", required=True,
                                help="Cluster file containing pre-clustered groups from CD-HIT, MMseqs, BLAST etc.")
    partial_parser.add_argument("-original_fasta", required=True,
                                help="FASTA file used in pre-clustering (Provide sequences in DNA form).")
    partial_parser.add_argument("-output_dir", required=True,
                                help="Directory for all output files.")
    partial_parser.add_argument("-reclustered", required=False,
                                help="Clustering output file from a second round of clustering.")
    partial_parser.add_argument("-seq_tag", default="StORF", dest="sequence_tag", required=False,
                                help="Legacy compatibility label; reclustering provenance is inferred from first-round membership.")

    # Common Grouping Arguments
    for subparser in [full_parser, partial_parser]:
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
                                 help="Memory budget in MB. In Partial mode this controls chunked loading of large inputs; in Full mode it is passed to CD-HIT.")
        subparser.add_argument("--sequence-buffer-cap-mb", type=int, default=256, dest="sequence_buffer_cap_mb",
                               help="Maximum chunk size in MB for the FASTA-backed sequence store used during loading.")
        subparser.add_argument("--tmp-dir", dest="tmp_dir", default=None,
                               help="Directory for temporary files such as chunked sequence stores and MAFFT scratch files.")
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



    # Setup logger once we know output paths/options
    # after we resolve output_path / options.output_dir:
    resolved_log_dir = options.log_dir if getattr(options, "log_dir", None) else (os.path.abspath(options.output_dir) if getattr(options, "output_dir", None) else os.getcwd())
    logger = configure_logger("PyamilySeq.PyamilySeq", enable_file=getattr(options, "log", False), log_dir=resolved_log_dir, verbose=options.verbose)
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
            # Proceed with the Partial mode
            pass
        else:
            missing_options = [opt for opt in
                               ['cluster_file','original_fasta'] if
                               not options.__dict__[opt]]
            sys.exit(f"Missing required options for Partial mode: {', '.join(missing_options)}")
        if options.align_core:
            options.write_individual_groups = True
            if options.write_groups == None or options.original_fasta == None:
                sys.exit('Must provide "-w" to output gene groups before alignment "-a" can be done.')


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
    if options.run_mode == 'Full':
        clustering_mode = 'cd-hit-est' if options.sequence_type == 'DNA' else 'cd-hit'
        if is_tool_installed(clustering_mode):
            logger.info("%s is installed. Proceeding with clustering.", clustering_mode)
            if options.fast_mode == True:
                options.fast_mode = 0
                logger.info("Running CD-HIT in fast mode.")
            else:
                options.fast_mode = 1
                logger.info("Running CD-HIT in accurate mode.")
        else:
            logger.error("%s is not installed. Please install CD-HIT to proceed.", clustering_mode)
            sys.exit(f"{clustering_mode} is not installed. Please install CD-HIT to proceed.")

    if hasattr(options, 'cluster_file') and options.cluster_file:
        options.cluster_file = fix_path(options.cluster_file)
    if hasattr(options, 'reclustered') and options.reclustered:
        options.reclustered = fix_path(options.reclustered)
    if hasattr(options, 'input_dir') and options.input_dir:
        options.input_dir = fix_path(options.input_dir)
    if hasattr(options, 'output_dir') and options.output_dir:
        options.output_dir = fix_path(options.output_dir)
    if getattr(options, 'taxon_map', None):
        options.taxon_map = fix_path(options.taxon_map)
    if getattr(options, 'sequence_buffer_cap_mb', None) is not None and options.sequence_buffer_cap_mb <= 0:
        sys.exit("--sequence-buffer-cap-mb must be greater than 0")
    if getattr(options, 'tmp_dir', None):
        options.tmp_dir = fix_path(options.tmp_dir)
        os.makedirs(options.tmp_dir, exist_ok=True)

    output_path = os.path.abspath(options.output_dir)
    combined_out_file = os.path.join(output_path, "combined_sequences_dna.fasta")
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
            file_to_cluster = combined_out_file.replace('_dna.fasta','_aa.fasta')
            translate = True
        elif options.sequence_type == 'DNA':
            clustering_mode = 'cd-hit-est'
            translate = False
            file_to_cluster = combined_out_file
        if options.input_type == 'separate':
            read_separate_files(options.input_dir, options.name_split_gff, options.name_split_fasta, options.gene_ident, combined_out_file, translate, False)
            run_cd_hit(options, file_to_cluster, clustering_output, clustering_mode)
        elif options.input_type == 'combined':
            read_combined_files(options.input_dir, options.name_split_gff, options.gene_ident, combined_out_file, translate, False)
            run_cd_hit(options, file_to_cluster, clustering_output, clustering_mode)
        elif options.input_type == 'fasta':
            combined_out_file = options.input_fasta
            # Detect if the input FASTA file contains DNA or AA sequences
            is_dna = detect_sequence_type(options.input_fasta)
            # If the sequence type is AA and the input is DNA, translate the DNA to AA
            if options.sequence_type == 'AA' and is_dna:
                translated_fasta = os.path.join(output_path, os.path.splitext(os.path.basename(options.input_fasta))[0] + '_aa.fasta')
                translate_dna_to_aa(options.input_fasta, translated_fasta)
                file_to_cluster = translated_fasta
            else:
                file_to_cluster = options.input_fasta
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

        clustering_options = clustering_options()

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
                "Out of memory during species clustering. Partial mode currently loads the full FASTA and "
                "cluster structures into RAM, so extremely large inputs can exhaust available memory."
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
                "Out of memory during genus clustering. Partial mode currently loads the full FASTA and "
                "cluster structures into RAM, so extremely large inputs can exhaust available memory."
            )
            sys.exit(1)
        except Exception:
            logger.exception("Unexpected error during genus clustering")
            sys.exit(1)

    # Save arguments to a text file
    from datetime import datetime
    with open(output_path+"/PyamilySeq_params.txt", "w") as outfile:
        outfile.write(f"Timestamp: {datetime.now().isoformat()}\n")
        for arg, value in vars(options).items():
            outfile.write(f"{arg}: {value}\n")
    logger.info("Saved parameters to %s", os.path.join(output_path, "PyamilySeq_params.txt"))



if __name__ == "__main__":
    main()
