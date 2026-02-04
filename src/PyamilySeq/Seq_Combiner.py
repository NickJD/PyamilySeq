try:
    from .constants import *
    from .utils import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from constants import *
    from utils import *

import threading
import time
import os
from typing import Optional
import re

def count_matching_files(input_dir: str, name_split: Optional[str], extensions):
    """
    Count input files in input_dir that match the provided extensions and, if name_split supplied,
    contain the name_split substring in the filename. This is used to compute total work units (files).
    """
    if not input_dir or not os.path.isdir(input_dir):
        return 0
    total = 0
    for fname in os.listdir(input_dir):
        low = fname.lower()
        if any(low.endswith(ext) for ext in extensions):
            if name_split:
                if name_split in fname:
                    total += 1
            else:
                total += 1
    return total

def count_files_present_in_combined(combined_file: str, name_split: Optional[str]) -> int:
    """
    Heuristic: count number of distinct input files (genomes) already present in the combined output.
    Primary approach: parse headers and take the second '|' field (header.split('|')[1]) as genome/file id.
    If that parsing fails, look for tokens containing name_split inside the header.
    """
    if not combined_file or not os.path.exists(combined_file):
        return 0
    seen = set()
    try:
        with open(combined_file, 'r') as fh:
            for line in fh:
                if not line.startswith('>'):
                    continue
                header = line[1:].strip()
                # 1) Prefer headers like ">id|genome|rest" -> take genome (second field)
                if '|' in header:
                    parts = header.split('|')
                    if len(parts) > 1 and parts[1]:
                        seen.add(parts[1])
                        continue
                # 2) If name_split provided, look for a filename-like token that includes it
                if name_split:
                    match = re.search(r'([^\s/\\]*' + re.escape(name_split) + r'[^\s/\\]*)', header)
                    if match:
                        token = os.path.basename(match.group(1))
                        seen.add(token)
                        continue
                # 3) If nothing matched, skip this header (avoids per-sequence overcounting)
    except Exception:
        return 0
    return len(seen)

# Helpers for progress reporting

def progress_reporter(stop_event, logger, total_files, combined_file, name_split=None, interval=10):
    """
    Periodically log progress. Preference: count headers in combined_file.
    Falls back to simple heartbeat if combined_file isn't yet created.
    """
    start = time.time()
    while not stop_event.is_set():
        # Use number of distinct input files represented in the combined output for "processed"
        processed = count_files_present_in_combined(combined_file, name_split) if combined_file else 0
        # Cap processed to total_files (prevents >100%)
        if total_files > 0 and processed > total_files:
            processed = total_files
        pct = (processed / total_files * 100) if total_files > 0 else 0.0
        elapsed = time.time() - start
        logger.info("Progress: %d/%d processed (%.1f%%). Elapsed: %.0fs", processed, total_files, pct, elapsed)
        # Wait with early exit support
        stop_event.wait(interval)
    # Final log when exiting
    processed = count_files_present_in_combined(combined_file, name_split) if combined_file else 0
    if total_files > 0 and processed > total_files:
        processed = total_files
    pct = (processed / total_files * 100) if total_files > 0 else 0.0
    elapsed = time.time() - start
    logger.info("Final progress: %d/%d processed (%.1f%%). Total elapsed: %.0fs", processed, total_files, pct, elapsed)

def main():
    # Early console-only logger so parser.description is logged before help/usage.
    early_logger = configure_logger("PyamilySeq.Seq_Combiner", enable_file=False, log_dir=None, verbose=False)
    parser = LoggingArgumentParser(logger_name="PyamilySeq.Seq_Combiner", description='Running Seq-Combiner - A tool to extract sequences from GFF/FASTA files and prepare them for PyamilySeq.')
    ### Required Arguments
    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-input_dir', action='store', dest='input_dir',
                          help='Directory location where the files are located.',
                          required=True)
    required.add_argument('-input_type', action='store', dest='input_type', choices=['separate', 'combined', 'fasta'],
                          help='Type of input files: "separate" for separate FASTA and GFF files,'
                             ' "combined" for GFF files with embedded FASTA sequences and "fasta" for combining multiple '
                               'FASTA files together.',
                          required=True)
    required.add_argument("-name_split_gff", action="store", dest="name_split_gff",
                          help="Substring used to split the filename and extract the genome name ('_combined.gff3' or '.gff'). - Not needed with -input_type fasta",
                          required=False)
    required.add_argument("-name_split_fasta", action="store", dest="name_split_fasta",
                          help="Substring used to split filenames and extract genome names for fasta files if named differently to paired gff files (e.g., '_dna.fasta').",
                          required=False)
    required.add_argument("-output_dir", action="store", dest="output_dir",
                          help="Directory for all output files.",
                          required=True)
    required.add_argument("-output_name", action="store", dest="output_file",
                          help="Output file name.",
                          required=True)

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-gene_ident', action='store', dest='gene_ident', default='CDS',
                          help='Default - "CDS": Identifier used for extraction of sequences such as "misc_RNA,gene,mRNA,CDS,rRNA,tRNA,tmRNA,CRISPR,ncRNA,regulatory_region,oriC,pseudo"'
                               ' - Not compatible with "fasta" input mode.',
                          required=False)
    optional.add_argument('-translate', action='store_true', dest='translate', default=None,
                          help='Default - False: Translate extracted sequences to their AA counterpart? - appends _aa.fasta to given output_name',
                          required=False)
    misc = parser.add_argument_group('Misc Arguments')
    misc.add_argument("-v", "--version", action="version",
                      version=f"PyamilySeq: Seq-Combiner version {PyamilySeq_Version} - Exiting",
                      help="Print out version number and exit")
    parser.add_argument("--log", action="store_true", dest="log", help="Create a timestamped logfile for this run.")
    parser.add_argument("--log-dir", dest="log_dir", default=None, help="Directory for logfile (default: output_dir).")

    options = parser.parse_args()

    # Setup logger for Seq-Combiner
    output_path = os.path.abspath(options.output_dir)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    log_dir = options.log_dir if getattr(options, "log_dir", None) else output_path
    logger = configure_logger("PyamilySeq.Seq_Combiner", enable_file=getattr(options, "log", False), log_dir=log_dir, verbose=False)

    # --- Progress reporting setup ------------------------------------------------
    combined_out_file = os.path.join(output_path, options.output_file)
    # Determine name_split and extensions per mode and count matching input files as total work units
    if options.input_type == 'fasta':
        name_split = options.name_split_fasta
        exts = ('.fasta', '.fa', '.fna')
    else:  # 'separate' or 'combined'
        name_split = options.name_split_gff
        exts = ('.gff', '.gff3', '.gff.gz', '.gff3.gz')

    total_work = count_matching_files(options.input_dir, name_split, exts)
    logger.info("Found %d input files (matching pattern) to process in %s", total_work, options.input_dir)

    stop_event = threading.Event()
    reporter_thread = threading.Thread(target=progress_reporter, args=(stop_event, logger, total_work, combined_out_file, name_split, 10), daemon=True)
    reporter_thread.start()
    # ---------------------------------------------------------------------------

    if options.input_type == 'separate' and options.name_split_gff is None:
        logger.error("Please provide a substring to split the filename and extract the genome name.")
        print("Please provide a substring to split the filename and extract the genome name.")
        exit(1)
    if options.input_type == 'combined' and options.name_split_gff is None:
        logger.error("Please provide a substring to split the filename and extract the genome name.")
        print("Please provide a substring to split the filename and extract the genome name.")
        exit(1)
    if options.input_type == 'fasta' and options.name_split_fasta is None:
        logger.error("Please provide a substring to split the filename and extract the genome name.")
        print("Please provide a substring to split the filename and extract the genome name.")
        exit(1)

    if options.input_type == 'separate':
        logger.info("Processing 'separate' input_type from %s", options.input_dir)
        read_separate_files(options.input_dir, options.name_split_gff, options.name_split_fasta, options.gene_ident, combined_out_file, options.translate, True)
    elif options.input_type == 'combined':
        logger.info("Processing 'combined' input_type from %s", options.input_dir)
        read_combined_files(options.input_dir, options.name_split_gff, options.gene_ident, combined_out_file, options.translate, True)
    elif options.input_type == 'fasta':
        logger.info("Processing 'fasta' input_type from %s", options.input_dir)
        read_fasta_files(options.input_dir, options.name_split_fasta, combined_out_file, options.translate, True)
    logger.info("Seq-Combiner completed.")

    # Stop reporter and wait for final log
    stop_event.set()
    reporter_thread.join(timeout=5)
    # Final summary: count number of input files represented (heuristic)
    final_files = count_files_present_in_combined(combined_out_file, name_split)
    logger.info("Completed combining. Final combined file: %s (input files represented: %d)", combined_out_file, final_files)

if __name__ == "__main__":
    main()
