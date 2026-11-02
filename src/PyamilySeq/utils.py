import subprocess
import shutil
import os
import glob
import collections
from tempfile import NamedTemporaryFile
import sys
import re
import math
import logging
import argparse
import gzip
import time
import secrets
import csv
import io
from urllib.parse import unquote
from contextlib import contextmanager

# Strict IUPAC DNA alphabet.
# Lowercase characters are accepted by converting to uppercase for validation.
DNA_IUPAC_BASES = frozenset(
    "ACGTRYSWKMBDHVN"
)

logger = logging.getLogger("PyamilySeq")  # Use the shared top-level PyamilySeq logger so all utils logs propagate to the same handlers

_startup_messages_pending = []


def emit_pending_startup_messages():
    global _startup_messages_pending
    if _startup_messages_pending:
        try:
            for msg in _startup_messages_pending:
                try:
                    logger.info("%s", msg)
                except Exception:
                    # swallow any logging errors to avoid breaking flow
                    pass
        finally:
            _startup_messages_pending.clear()

####
# Placeholder for the distance function
levenshtein_distance_cal = None
# Check for Levenshtein library once
try:
    import Levenshtein as LV
    # Assign the optimised function
    def levenshtein_distance_calc(seq1, seq2):
        return LV.distance(seq1, seq2)
except (ModuleNotFoundError, ImportError):
    # Save the notice for later emission (after logger handlers are configured)
    _startup_messages_pending.append("Levenshtein package not installed - Will fallback to slower Python implementation.")
    # Fallback implementation
    def levenshtein_distance_calc(seq1, seq2):
        # Slower Python implementation of Levenshtein distance
        len1, len2 = len(seq1), len(seq2)
        dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]

        for i in range(len1 + 1):
            dp[i][0] = i
        for j in range(len2 + 1):
            dp[0][j] = j

        for i in range(1, len1 + 1):
            for j in range(1, len2 + 1):
                if seq1[i - 1] == seq2[j - 1]:
                    cost = 0
                else:
                    cost = 1
                dp[i][j] = min(dp[i - 1][j] + 1,  # Deletion
                               dp[i][j - 1] + 1,  # Insertion
                               dp[i - 1][j - 1] + cost)  # Substitution

        return dp[len1][len2]
#####

################### We are currently fixed using Table 11
codon_table = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
      'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

def translate_frame(sequence):
    sequence = ''.join(str(sequence).split()).upper().replace('U', 'T')
    translate = ''.join([codon_table.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate

def _open_text(path, mode='rt'):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    return opener(
        path,
        mode,
        encoding="utf-8",
        errors="strict",
    )

@contextmanager
def _atomic_text_outputs(final_paths):
    """Write UTF-8 outputs to temporary files, then atomically replace targets.

    A target ending in ``.gz`` is gzip-compressed while streaming, avoiding a
    second full-size uncompressed annotation sidecar on disk. Temporary files
    are created beside each destination so ``os.replace`` stays on one filesystem.
    """
    handles = {}
    raw_handles = {}
    temporary_paths = {}
    committed = set()

    try:
        for final_path in final_paths:
            if final_path is None:
                if None in handles:
                    raise ValueError("Only one null output handle is supported.")
                handle = open(os.devnull, "w", encoding="utf-8", errors="strict", newline="\n")
                handles[None] = handle
                raw_handles[None] = handle
                continue

            final_path = os.fspath(final_path)
            if final_path in handles:
                raise ValueError(f"Duplicate output path supplied: {final_path}")

            destination_dir = os.path.dirname(os.path.abspath(final_path))
            prefix = f".{os.path.basename(final_path)}."
            for _attempt in range(100):
                temporary_path = os.path.join(
                    destination_dir, prefix + secrets.token_hex(8) + ".tmp"
                )
                try:
                    fd = os.open(
                        temporary_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666
                    )
                    break
                except FileExistsError:
                    continue
            else:
                raise FileExistsError(
                    f"Could not create a unique temporary output beside {final_path}"
                )

            temporary_paths[final_path] = temporary_path
            if str(final_path).lower().endswith(".gz"):
                raw_handle = os.fdopen(fd, "wb")
                gzip_handle = gzip.GzipFile(fileobj=raw_handle, mode="wb", mtime=0)
                text_handle = io.TextIOWrapper(
                    gzip_handle, encoding="utf-8", errors="strict", newline="\n"
                )
                handles[final_path] = text_handle
                raw_handles[final_path] = raw_handle
            else:
                handle = os.fdopen(
                    fd, "w", encoding="utf-8", errors="strict", newline="\n"
                )
                handles[final_path] = handle
                raw_handles[final_path] = handle

        yield handles

        # Flush and finalise compression streams before syncing their underlying
        # file descriptors. GzipFile writes its footer when closed.
        for final_path, handle in handles.items():
            if final_path is None:
                handle.flush()
                handle.close()
                continue
            handle.flush()
            if handle is not raw_handles[final_path]:
                handle.close()
            raw_handle = raw_handles[final_path]
            raw_handle.flush()
            os.fsync(raw_handle.fileno())
            if not raw_handle.closed:
                raw_handle.close()

        for final_path, temporary_path in temporary_paths.items():
            os.replace(temporary_path, final_path)
            committed.add(final_path)

    finally:
        for handle in handles.values():
            if not handle.closed:
                handle.close()
        for final_path, handle in raw_handles.items():
            if handle not in handles.values() and not handle.closed:
                handle.close()
        for final_path, temporary_path in temporary_paths.items():
            if final_path not in committed and os.path.exists(temporary_path):
                try:
                    os.unlink(temporary_path)
                except OSError:
                    logger.warning(
                        "Could not remove incomplete temporary output: %s", temporary_path
                    )


def _sorted_input_files(input_dir, suffix, *, gzip_fallback=False):
    """Discover files in a deterministic, basename-sorted order."""
    pattern = os.path.join(os.fspath(input_dir), "*" + str(suffix))
    files = [item for item in glob.glob(pattern) if os.path.isfile(item)]

    if not files and gzip_fallback and not str(suffix).lower().endswith(".gz"):
        compressed_pattern = pattern + ".gz"
        files = [
            item for item in glob.glob(compressed_pattern)
            if os.path.isfile(item)
        ]

    return sorted(files, key=lambda item: os.path.basename(item))


def _iter_fasta_records(fasta_path):
    """Yield ``(full_header, sequence)`` records with basic FASTA validation."""
    current_header = None
    current_id = None
    sequence_chunks = []
    seen_ids = set()
    record_count = 0

    with _open_text(fasta_path, "rt") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    sequence = "".join(sequence_chunks)
                    if not sequence:
                        raise ValueError(
                            f"Empty FASTA sequence for {current_id!r} in "
                            f"{fasta_path} before line {line_number}."
                        )
                    yield current_header, sequence
                    record_count += 1

                header_text = line[1:].strip()
                if not header_text:
                    raise ValueError(
                        f"Empty FASTA header at line {line_number} in {fasta_path}."
                    )
                current_id = header_text.split(None, 1)[0]
                if current_id in seen_ids:
                    raise ValueError(
                        f"Duplicate FASTA identifier {current_id!r} at line "
                        f"{line_number} in {fasta_path}."
                    )
                seen_ids.add(current_id)
                current_header = ">" + header_text
                sequence_chunks = []
            else:
                if current_header is None:
                    raise ValueError(
                        f"Sequence data before the first FASTA header at line "
                        f"{line_number} in {fasta_path}."
                    )
                sequence_chunks.append("".join(line.split()))

    if current_header is not None:
        sequence = "".join(sequence_chunks)
        if not sequence:
            raise ValueError(
                f"Empty FASTA sequence for {current_id!r} at end of {fasta_path}."
            )
        yield current_header, sequence
        record_count += 1

    if record_count == 0:
        raise ValueError(f"No FASTA records found in {fasta_path}.")


def _read_fasta_sequences_by_id(fasta_path):
    """Read a per-genome FASTA into an insertion-ordered contig dictionary."""
    sequences = {}
    for header, sequence in _iter_fasta_records(fasta_path):
        sequence_id = header[1:].split(None, 1)[0]
        sequences[sequence_id] = sequence
    return sequences


def _extract_gff_feature_sequence(
    contig_sequences,
    *,
    contig,
    start,
    end,
    strand,
    context,
    circular_contigs=None,
):
    """Extract and validate one 1-based inclusive GFF feature.

    GFF3 permits an origin-spanning feature on a circular landmark to have an
    end coordinate greater than the underlying sequence length. Such
    coordinates are virtual coordinates: the sequence is extracted from the
    end of the contig and continued from its beginning. This exception is
    enabled only when the GFF explicitly marks the contig ``Is_circular=true``.
    """
    if contig not in contig_sequences:
        raise ValueError(
            f"GFF feature references contig {contig!r}, which is absent from "
            f"the associated FASTA. {context}"
        )

    if start < 1 or end < start:
        raise ValueError(
            f"Invalid GFF coordinates {start}-{end}; expected 1 <= start <= end. "
            f"{context}"
        )

    contig_sequence = contig_sequences[contig]
    contig_length = len(contig_sequence)
    if contig_length == 0:
        raise ValueError(f"Contig {contig!r} has an empty FASTA sequence. {context}")

    feature_length = end - start + 1
    is_circular = circular_contigs is not None and contig in circular_contigs
    coordinates_exceed_sequence = start > contig_length or end > contig_length

    if coordinates_exceed_sequence and not is_circular:
        raise ValueError(
            f"GFF feature coordinates {start}-{end} exceed the length of linear "
            f"contig {contig!r} ({contig_length:,} bp). Origin-spanning coordinates "
            f"are accepted only when the GFF landmark explicitly declares "
            f"Is_circular=true. {context}"
        )

    if coordinates_exceed_sequence and is_circular:
        if feature_length > contig_length:
            raise ValueError(
                f"Circular GFF feature spans {feature_length:,} bp, longer than "
                f"contig {contig!r} ({contig_length:,} bp); extraction would require "
                f"more than one traversal of the circular sequence. {context}"
            )
        # Convert the virtual 1-based coordinate to a zero-based position on
        # the circular sequence, then concatenate the tail and head as needed.
        start_offset = (start - 1) % contig_length
        stop_offset = start_offset + feature_length
        if stop_offset <= contig_length:
            sequence = contig_sequence[start_offset:stop_offset]
        else:
            sequence = (
                contig_sequence[start_offset:]
                + contig_sequence[:stop_offset - contig_length]
            )
    else:
        sequence = contig_sequence[start - 1:end]

    if len(sequence) != feature_length:
        raise ValueError(
            f"Extracted sequence length ({len(sequence):,}) does not match GFF "
            f"coordinate length ({feature_length:,}). {context}"
        )

    if strand == "-":
        sequence = reverse_complement(sequence)
    elif strand != "+":
        raise ValueError(
            f"Unsupported CDS strand {strand!r}; expected '+' or '-'. {context}"
        )

    return validate_nucleotide_sequence(
        sequence,
        context=context,
        start=start,
        end=end,
        strand=strand,
    )


def validate_nucleotide_sequence(
    sequence,
    *,
    context,
    start=None,
    end=None,
    strand=None,
):
    """
    Validate an extracted DNA sequence against the IUPAC DNA alphabet.

    Parameters
    ----------
    sequence : str
        Extracted nucleotide sequence.
    context : str
        Descriptive source information, e.g. genome, gene, GFF file.
    start, end : int, optional
        Original 1-based inclusive GFF feature coordinates.
    strand : str, optional
        Feature strand ('+' or '-').

    Returns
    -------
    str
        Original sequence, unchanged.

    Raises
    ------
    ValueError
        If any character is outside the accepted DNA alphabet.
    """

    if not sequence:
        raise ValueError(f"Empty nucleotide sequence encountered: {context}")

    # Fast path: set construction and upper-casing are implemented in C. Avoid
    # a Python-level per-base loop for the many millions of valid CDS records.
    invalid_symbols = set(sequence.upper()).difference(DNA_IUPAC_BASES)
    if not invalid_symbols:
        return sequence

    invalid = []
    for position, base in enumerate(sequence, start=1):
        if base.upper() not in DNA_IUPAC_BASES:
            genomic_position = None
            if start is not None and end is not None:
                if strand == "+":
                    genomic_position = start + position - 1
                elif strand == "-":
                    genomic_position = end - position + 1
            invalid.append((position, base, genomic_position))

    if invalid:

        examples = []

        for position, base, genomic_position in invalid[:10]:

            detail = (
                f"{base!r} (U+{ord(base):04X}) "
                f"at CDS position {position}"
            )

            if genomic_position is not None:
                detail += (
                    f", genomic coordinate {genomic_position}"
                )

            examples.append(detail)

        details = "; ".join(examples)

        raise ValueError(
            "\n"
            "STRICT NUCLEOTIDE VALIDATION FAILED\n"
            f"Source: {context}\n"
            f"Sequence length: {len(sequence):,} bp\n"
            f"Invalid character count: {len(invalid):,}\n"
            f"Invalid characters: {details}\n"
            "Accepted DNA alphabet: ACGTRYSWKMBDHVN\n"
            "No sequences were automatically modified."
        )

    return sequence

def normalise_gene_identifiers(gene_ident):
    if gene_ident is None:
        return []
    if isinstance(gene_ident, str):
        return [item.strip() for item in gene_ident.split(',') if item.strip()]
    return [str(item).strip() for item in gene_ident if str(item).strip()]


def feature_matches_gene_ident(feature_type, gene_ident):
    idents = normalise_gene_identifiers(gene_ident)
    if not idents:
        return False
    return str(feature_type).lower() in {ident.lower() for ident in idents}


def _extract_gff_identifier(attributes):
    parsed = {}
    for token in (part.strip() for part in attributes.strip().split(';') if part.strip()):
        if '=' in token:
            key, value = token.split('=', 1)
        elif ' ' in token:
            key, value = token.split(None, 1)
            value = value.strip('"')
        else:
            continue
        parsed[key.strip()] = value.strip().strip('"')
    for key in ('ID', 'locus_tag', 'protein_id', 'gene', 'Name', 'Parent'):
        if parsed.get(key):
            return parsed[key].split(',')[0]
    raise ValueError(f"GFF feature has no supported identifier attribute: {attributes.strip()}")


def _parse_gff_attributes(attributes):
    """Parse GFF3/GTF-like attributes and URL-decode their values."""
    parsed = {}
    for token in (part.strip() for part in str(attributes).strip().split(';') if part.strip()):
        if '=' in token:
            key, value = token.split('=', 1)
        elif ' ' in token:
            key, value = token.split(None, 1)
            value = value.strip().strip('"')
        else:
            continue
        key = unquote(key.strip())
        value = unquote(value.strip().strip('"'))
        if key in parsed and value:
            parsed[key] = parsed[key] + ',' + value
        else:
            parsed[key] = value
    return parsed


def _gff_declares_circular(attributes):
    """Return True when a GFF landmark explicitly declares ``Is_circular=true``.

    This lightweight parser is intentionally cheaper than decoding all GFF
    attributes for every feature in a multi-million-CDS input.
    """
    for token in str(attributes).split(";"):
        key, separator, value = token.partition("=")
        if separator and key.strip().casefold() == "is_circular":
            return value.strip().strip('"').casefold() == "true"
    return False


def _gff_gene_metadata(attributes):
    """Return (non-unique gene name, annotation) without inferring missing values.

    ``gene`` is the conventional Prokka/Bakta gene-symbol attribute. ``product``
    is the functional annotation; ``Name`` is the common fallback when product
    is absent. Protein descriptions are never used as gene symbols.
    """
    parsed = _parse_gff_attributes(attributes)
    gene_name = parsed.get('gene') or parsed.get('gene_name') or ''
    annotation = parsed.get('product') or parsed.get('Name') or parsed.get('Note') or ''
    if gene_name:
        gene_name = gene_name.split(',', 1)[0].strip()
    return gene_name, annotation.strip()


def _write_gene_annotation(writer, genome_name, sequence_id, gene_name, annotation, source_gff):
    if writer is None:
        return
    writer.writerow([
        f"{genome_name}|{sequence_id}",
        gene_name or "",
        annotation or "",
        os.path.basename(source_gff),
    ])


def iter_gene_annotations(path):
    """Yield validated gene-annotation sidecar rows as dictionaries."""
    with _open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sequence_id", "gene_name", "annotation"}
        if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                f"Gene-annotation file {path} must contain columns: "
                "sequence_id, gene_name, annotation."
            )
        for line_number, row in enumerate(reader, start=2):
            identifier = (row.get("sequence_id") or "").strip()
            if not identifier:
                raise ValueError(
                    f"Empty sequence_id in gene-annotation file {path} at line {line_number}."
                )
            yield {
                "sequence_id": identifier,
                "gene_name": (row.get("gene_name") or "").strip(),
                "annotation": (row.get("annotation") or "").strip(),
            }


def translate_dna_to_aa(dna_fasta, aa_fasta):
    """Translate a DNA FASTA after strict IUPAC validation.

    The amino-acid output is committed only if every input record is valid.
    """
    aa_path = os.fspath(aa_fasta)
    record_count = 0

    with _atomic_text_outputs([aa_path]) as outputs:
        outfile = outputs[aa_path]
        for header, dna_sequence in _iter_fasta_records(dna_fasta):
            validate_nucleotide_sequence(
                dna_sequence,
                context=f"FASTA={dna_fasta}; record={header}",
            )
            aa_sequence = translate_frame(dna_sequence)
            outfile.write(f"{header}\n{wrap_sequence(aa_sequence, 60)}\n")
            record_count += 1

    if record_count == 0:
        # _iter_fasta_records already raises for an empty file; retained as a
        # defensive assertion if the parser's behaviour changes in future.
        raise ValueError(f"No DNA FASTA records were translated from {dna_fasta}.")


def detect_sequence_type(fasta_file):
    observed = []
    with _open_text(fasta_file, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                continue
            observed.extend(character for character in line.strip().upper() if character.isalpha())
            if len(observed) >= 10000:
                break
    if not observed:
        raise ValueError(f"No sequence characters found in {fasta_file}")
    amino_acid_only = set(observed) - set('ACGTUNRYKMSWBDHV')
    if amino_acid_only:
        return False
    canonical_fraction = sum(character in 'ACGTUN' for character in observed) / len(observed)
    return canonical_fraction >= 0.85


def is_tool_installed(tool_name):
    """Check if a tool is installed and available in PATH."""
    if shutil.which(tool_name) is None:
        return False

    # Try running the tool to ensure it's executable
    try:
        subprocess.run([tool_name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except subprocess.CalledProcessError:
        return True  # The tool is installed and ran, even if it returns an error code
    except FileNotFoundError:
        return False  # This shouldn't happen due to the earlier check

_comp = str.maketrans(
    "ATGCRYSWKMVBHDNatgcryswkmvbhdn-",
    "TACGYRSWMKVBHDNtacgyrswmkvbhdn-"
)

def reverse_complement(seq):
    return seq.translate(_comp)[::-1]

# def reverse_complement(seq):
#     complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N','R': 'Y',
#                   'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'V': 'B',
#                   'B': 'V', 'H': 'D', 'D': 'H'}
#     return ''.join(complement[base] for base in reversed(seq))


def fix_path(path):
    fixed_path = os.path.normpath(path)
    fixed_path = os.path.realpath(fixed_path)
    return fixed_path


def extract_identity(clustered_info):
    # Use regex to capture percentage, including optional '-' or '+' before it
    match = re.search(r'at [+-/]*(\d+\.\d+)%', clustered_info)

    if match:
        percent_identity = float(match.group(1))  # Extract the percentage value
        return percent_identity
    else:
        raise ValueError("Percent identity not found in the string.")


def wrap_sequence(sequence, width=60):
    wrapped_sequence = []
    for i in range(0, len(sequence), width):
        wrapped_sequence.append(sequence[i:i + width])
    return "\n".join(wrapped_sequence)


def read_genomes_from_fasta(fasta_file):
    """Return the genome identifiers in deterministic sorted order."""
    genomes = set()
    for header, _sequence in _iter_fasta_records(fasta_file):
        sequence_id = header[1:].split(None, 1)[0]
        parts = sequence_id.split("|")
        if len(parts) < 2 or not parts[1]:
            raise ValueError(
                f"FASTA identifier {sequence_id!r} in {fasta_file} does not "
                "contain the expected 'genome|sequence' structure."
            )
        genomes.add(parts[1])
    return sorted(genomes)


def read_fasta(fasta_file):
    sequences = {}
    current_sequence = None
    with _open_text(fasta_file, 'rt') as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            if line.startswith('>'):
                parts = line[1:].split(maxsplit=1)
                if not parts:
                    raise ValueError(f"Empty FASTA header at line {line_number} in {fasta_file}")
                current_sequence = parts[0]
                if current_sequence in sequences:
                    raise ValueError(
                        f"Duplicate FASTA identifier {current_sequence!r} at line {line_number} in {fasta_file}"
                    )
                sequences[current_sequence] = ''
            else:
                if current_sequence is None:
                    raise ValueError(
                        f"Sequence data before first FASTA header at line {line_number} in {fasta_file}"
                    )
                sequences[current_sequence] += line
    return sequences


def reorder_dict_by_keys(original_dict, sorted_keys):
    return {k: original_dict[k] for k in sorted_keys}
def custom_sort_key(k, dict1, dict2):
    return (len(dict1[k]), len(dict2[k]))

def sort_keys_by_values(dict1, dict2):
    sorted_keys = sorted(dict1.keys(), key=lambda k: custom_sort_key(k, dict1, dict2), reverse=True)
    return sorted_keys

def select_longest_gene(sequences, subgrouped):
    """Select the longest sequence for each genome."""
    longest_sequences = {}
    for seq_id, sequence in sequences.items():
        if subgrouped == False:
            genome = seq_id.split('|')[0]  # Assuming genome name can be derived from the sequence ID
        elif subgrouped == True:
            genome = seq_id.split('|')[1]
        if genome not in longest_sequences or len(sequence) > len(longest_sequences[genome][1]):
            longest_sequences[genome] = (seq_id, sequence)
    return longest_sequences


def run_mafft_on_sequences(options, sequences, output_file):
    """Run mafft on the given sequences and write to output file."""
    emit_pending_startup_messages()
    # Create a temporary input file for mafft
    temp_dir = getattr(options, "tmp_dir", None)
    with NamedTemporaryFile('w', delete=False, dir=temp_dir) as temp_input_file:
        for header, sequence in sequences.items():
            temp_input_file.write(f">{header}\n{sequence}\n")
        temp_input_file_path = temp_input_file.name

    # Run mafft
    try:
        with open(output_file, 'w') as output_f:
            if getattr(options, "verbose", False):
                subprocess.run(
                    ['mafft', '--auto', '--thread', str(options.threads), temp_input_file_path],
                    stdout=output_f,
                    stderr=sys.stderr,
                    check=True
                )
            else:
                subprocess.run(
                    ['mafft', '--auto', '--thread', str(options.threads), temp_input_file_path],
                    stdout=output_f,
                    stderr=subprocess.DEVNULL,  # Suppress stderr
                    check=True
                )
    finally:
        try:
            os.remove(temp_input_file_path)  # Clean up the temporary file
        except Exception:
            pass



def read_separate_files(
    input_dir,
    name_split_gff,
    name_split_fasta,
    gene_ident,
    combined_out,
    translate,
    run_as_combiner,
    annotation_out=None,
):
    """Extract CDS records from paired GFF and FASTA files deterministically.

    GFF files are processed in basename-sorted order. Every extracted CDS is
    checked against the IUPAC DNA alphabet and its source coordinates. Outputs
    are committed only after all input files have been processed successfully.
    """
    emit_pending_startup_messages()
    dna_output_path = os.fspath(combined_out)
    requested_aa_path = None
    if not run_as_combiner:
        requested_aa_path = dna_output_path.replace("_dna.fasta", "_aa.fasta")
        if requested_aa_path == dna_output_path:
            raise ValueError(
                "The combined output filename must contain '_dna.fasta' when "
                "run_as_combiner=False, so a separate '_aa.fasta' path can be made."
            )

    aa_output_path = requested_aa_path if translate else None
    annotation_output_path = os.fspath(annotation_out) if annotation_out else None
    output_paths = [dna_output_path, aa_output_path]
    if annotation_output_path:
        output_paths.append(annotation_output_path)
    gff_files = _sorted_input_files(
        input_dir, name_split_gff, gzip_fallback=True
    )
    if not gff_files:
        raise FileNotFoundError(
            f"No GFF files found in {input_dir} (pattern: *{name_split_gff})."
        )

    logger.info(
        "Extracting paired GFF/FASTA inputs in deterministic basename order: %s GFF files.",
        f"{len(gff_files):,}",
    )

    paired_files_found = False
    missing_pairs = []
    with _atomic_text_outputs(output_paths) as outputs:
        dna_handle = outputs[dna_output_path]
        aa_handle = outputs[aa_output_path]
        annotation_handle = outputs[annotation_output_path] if annotation_output_path else None
        annotation_writer = csv.writer(annotation_handle, delimiter="\t", lineterminator="\n") if annotation_handle else None
        if annotation_writer:
            annotation_writer.writerow(["sequence_id", "gene_name", "annotation", "source_gff"])

        for gff_file in gff_files:
            gff_basename = os.path.basename(gff_file)
            if gff_basename.endswith(name_split_gff):
                genome_name = gff_basename[:-len(name_split_gff)] if name_split_gff else gff_basename
            else:
                # This is the .gz fallback when name_split_gff excludes '.gz'.
                suffix = str(name_split_gff) + ".gz"
                genome_name = gff_basename[:-len(suffix)] if gff_basename.endswith(suffix) else os.path.splitext(gff_basename)[0]

            if name_split_fasta is None:
                # Strip both compression and annotation extensions before trying
                # the conventional per-genome FASTA suffixes.
                fasta_stem = gff_basename
                if fasta_stem.lower().endswith(".gz"):
                    fasta_stem = fasta_stem[:-3]
                fasta_stem = os.path.splitext(fasta_stem)[0]
                possible_extensions = (
                    ".fa", ".fasta", ".fna",
                    ".fa.gz", ".fasta.gz", ".fna.gz",
                )
                corresponding_fasta_file = next(
                    (
                        os.path.join(input_dir, fasta_stem + extension)
                        for extension in possible_extensions
                        if os.path.isfile(os.path.join(input_dir, fasta_stem + extension))
                    ),
                    None,
                )
                if corresponding_fasta_file is None:
                    logger.warning(
                        "Corresponding FASTA for GFF '%s' was not found; this pair cannot be processed.",
                        gff_file,
                    )
                    missing_pairs.append(gff_file)
                    continue
            else:
                corresponding_fasta_file = os.path.join(
                    input_dir, genome_name + str(name_split_fasta)
                )
                if not os.path.isfile(corresponding_fasta_file):
                    logger.warning(
                        "Corresponding FASTA for GFF '%s' was not found: expected '%s'.",
                        gff_file, corresponding_fasta_file,
                    )
                    missing_pairs.append(gff_file)
                    continue

            gff_features = []
            circular_contigs = set()
            with _open_text(gff_file, "rt") as gff_handle:
                seen_sequence_ids = collections.defaultdict(int)
                for line_number, raw_line in enumerate(gff_handle, start=1):
                    line = raw_line.rstrip("\r\n")
                    if not line or line.startswith("#"):
                        continue
                    fields = line.split("\t")
                    if len(fields) != 9:
                        raise ValueError(
                            f"Malformed GFF3 record: expected 9 tab-delimited fields, "
                            f"found {len(fields)}. GFF={gff_file}; line={line_number}"
                        )
                    if _gff_declares_circular(fields[8]):
                        circular_contigs.add(fields[0])
                    if not feature_matches_gene_ident(fields[2], gene_ident):
                        continue
                    try:
                        start = int(fields[3])
                        end = int(fields[4])
                    except ValueError as error:
                        raise ValueError(
                            f"Non-integer GFF coordinates. GFF={gff_file}; line={line_number}"
                        ) from error

                    try:
                        sequence_id = _extract_gff_identifier(fields[8])
                    except ValueError as error:
                        raise ValueError(
                            f"Could not identify GFF feature. GFF={gff_file}; line={line_number}: {error}"
                        ) from error

                    original_id = sequence_id
                    if original_id in seen_sequence_ids:
                        sequence_id = f"{original_id}_{seen_sequence_ids[original_id]}"
                        seen_sequence_ids[original_id] += 1
                    else:
                        seen_sequence_ids[original_id] = 1

                    gene_name, annotation = _gff_gene_metadata(fields[8])
                    gff_features.append(
                        (fields[0], start, end, fields[6], sequence_id, line_number,
                         gene_name, annotation)
                    )

            contig_sequences = _read_fasta_sequences_by_id(corresponding_fasta_file)
            paired_files_found = True

            for (contig, start, end, strand, sequence_id, line_number,
                 gene_name, annotation) in gff_features:
                context = (
                    f"GFF={gff_file}; FASTA={corresponding_fasta_file}; "
                    f"genome={genome_name}; gene={sequence_id}; contig={contig}; "
                    f"coordinates={start}-{end}; strand={strand}; GFF_line={line_number}"
                )
                sequence = _extract_gff_feature_sequence(
                    contig_sequences,
                    contig=contig,
                    start=start,
                    end=end,
                    strand=strand,
                    context=context,
                    circular_contigs=circular_contigs,
                )

                _write_gene_annotation(
                    annotation_writer, genome_name, sequence_id,
                    gene_name, annotation, gff_file,
                )

                if run_as_combiner:
                    output_sequence = translate_frame(sequence) if translate else sequence
                    dna_handle.write(
                        f">{genome_name}|{sequence_id}\n{wrap_sequence(output_sequence, 60)}\n"
                    )
                else:
                    if translate:
                        aa_sequence = translate_frame(sequence)
                        aa_handle.write(
                            f">{genome_name}|{sequence_id}\n{wrap_sequence(aa_sequence, 60)}\n"
                        )
                    dna_handle.write(
                        f">{genome_name}|{sequence_id}\n{wrap_sequence(sequence, 60)}\n"
                    )

        if missing_pairs:
            examples = ", ".join(os.path.basename(item) for item in missing_pairs[:10])
            suffix = " ..." if len(missing_pairs) > 10 else ""
            raise FileNotFoundError(
                f"Missing paired FASTA files for {len(missing_pairs):,} GFF inputs: "
                f"{examples}{suffix}. No extraction outputs were committed."
            )
        if not paired_files_found:
            raise FileNotFoundError(
                "No matching GFF/FASTA pairs were found. Check -name_split_gff and "
                "-name_split_fasta, or provide the corresponding per-genome FASTA files."
            )

    # Remove an obsolete AA output only after the DNA extraction has succeeded.
    if requested_aa_path and not translate and os.path.exists(requested_aa_path):
        try:
            os.remove(requested_aa_path)
        except OSError as error:
            logger.warning("Could not remove obsolete amino-acid output %s: %s", requested_aa_path, error)



class _ExtractionProgress:
    """Lightweight standard-library progress display for multi-file extraction.

    Progress is measured in input GFF files, with an additional running count
    of successfully extracted gene sequences.  It writes to stderr so normal
    PyamilySeq log messages remain readable on stdout.  When stderr is not a
    terminal, it emits periodic 5-percent progress log messages instead.
    """

    def __init__(self, total, description, progress_logger):
        self.total = max(0, int(total))
        self.description = str(description)
        self.logger = progress_logger
        self.completed = 0
        self.records = 0
        self.current = ""
        self.started = time.monotonic()
        self.is_tty = bool(getattr(sys.stderr, "isatty", lambda: False)())
        self.next_log_percent = 5
        self.last_rendered_records = 0
        self.closed = False

        self.logger.info(
            "%s: processing %s input GFF files.",
            self.description, f"{self.total:,}"
        )
        if self.is_tty:
            self._render()

    @staticmethod
    def _elapsed(seconds):
        seconds = max(0, int(seconds))
        hours, remainder = divmod(seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        if hours:
            return f"{hours:d}:{minutes:02d}:{seconds:02d}"
        return f"{minutes:02d}:{seconds:02d}"

    def set_current(self, name):
        self.current = os.path.basename(str(name))
        if self.is_tty:
            self._render()

    def add_records(self, count=1):
        self.records += int(count)
        # Refresh the live terminal display every 1,000 genes, avoiding a costly
        # terminal write for every single extracted sequence.
        if self.is_tty and self.records - self.last_rendered_records >= 1000:
            self._render()

    def file_complete(self):
        self.completed = min(self.total, self.completed + 1)
        self.current = ""
        self._render()

    def _render(self, force=False):
        if self.total:
            fraction = min(1.0, self.completed / self.total)
        else:
            fraction = 1.0
        percent = fraction * 100.0
        elapsed = self._elapsed(time.monotonic() - self.started)

        if self.is_tty:
            width = 24
            filled = int(width * fraction)
            bar = "=" * filled + (">" if filled < width else "") + " " * max(0, width - filled - 1)
            current = self.current
            if len(current) > 28:
                current = current[:12] + "..." + current[-13:]
            suffix = f" | {current}" if current else ""
            message = (
                f"\r{self.description[:18]:18s} [{bar}] {percent:5.1f}% "
                f"{self.completed:,}/{self.total:,} GFFs | "
                f"{self.records:,} genes | {elapsed}{suffix}"
            )
            sys.stderr.write(message)
            sys.stderr.flush()
            self.last_rendered_records = self.records
            return

        if force or (percent >= self.next_log_percent and self.completed < self.total):
            self.logger.info(
                "%s: %.1f%% (%s/%s GFF files); %s genes extracted; elapsed %s.",
                self.description, percent, f"{self.completed:,}",
                f"{self.total:,}", f"{self.records:,}", elapsed,
            )
            while self.next_log_percent <= percent:
                self.next_log_percent += 5

    def close(self):
        if self.closed:
            return
        self.closed = True
        if self.is_tty:
            # file_complete() has already rendered the final 100% state.
            # If extraction was interrupted, render the partial state here.
            if self.completed < self.total:
                self._render()
            sys.stderr.write("\n")
            sys.stderr.flush()
        else:
            self._render(force=True)

        elapsed = self._elapsed(time.monotonic() - self.started)
        if self.completed == self.total:
            self.logger.info(
                "Finished %s: extracted %s genes from %s GFF files in %s.",
                self.description.lower(), f"{self.records:,}",
                f"{self.completed:,}", elapsed,
            )
        else:
            self.logger.warning(
                "%s stopped after %s of %s GFF files; %s genes had been extracted in %s.",
                self.description, f"{self.completed:,}", f"{self.total:,}",
                f"{self.records:,}", elapsed,
            )


def read_combined_files(
    input_dir,
    name_split,
    gene_ident,
    combined_out,
    translate,
    run_as_combiner,
    progress_logger=None,
    annotation_out=None,
):
    """Extract CDS records from combined GFF3 files deterministically.

    GFF inputs are processed in basename-sorted order. Embedded contig FASTA
    records and CDS coordinates are checked, and every extracted CDS is
    validated against the strict IUPAC DNA alphabet before it is written.
    Outputs are atomically replaced only after successful completion.
    """
    emit_pending_startup_messages()
    log = progress_logger or logging.getLogger("PyamilySeq.PyamilySeq")
    dna_output_path = os.fspath(combined_out)
    requested_aa_path = None
    if not run_as_combiner:
        requested_aa_path = dna_output_path.replace("_dna.fasta", "_aa.fasta")
        if requested_aa_path == dna_output_path:
            raise ValueError(
                "The combined output filename must contain '_dna.fasta' when "
                "run_as_combiner=False, so a separate '_aa.fasta' path can be made."
            )

    aa_output_path = requested_aa_path if translate else None
    annotation_output_path = os.fspath(annotation_out) if annotation_out else None
    output_paths = [dna_output_path, aa_output_path]
    if annotation_output_path:
        output_paths.append(annotation_output_path)
    gff_files = _sorted_input_files(input_dir, name_split, gzip_fallback=True)
    if not gff_files:
        raise FileNotFoundError(
            f"No GFF files found in {input_dir} (pattern: *{name_split})."
        )

    log.info(
        "Input GFF files will be processed in deterministic basename order (%s files).",
        f"{len(gff_files):,}",
    )

    if run_as_combiner:
        output_description = "Extracting CDS (translated AA)" if translate else "Extracting CDS (DNA)"
    elif translate:
        output_description = "Extracting CDS (DNA + AA)"
    else:
        output_description = "Extracting CDS (DNA)"

    progress = _ExtractionProgress(
        total=len(gff_files),
        description=output_description,
        progress_logger=log,
    )

    try:
        with _atomic_text_outputs(output_paths) as outputs:
            dna_handle = outputs[dna_output_path]
            aa_handle = outputs[aa_output_path]
            annotation_handle = outputs[annotation_output_path] if annotation_output_path else None
            annotation_writer = csv.writer(annotation_handle, delimiter="\t", lineterminator="\n") if annotation_handle else None
            if annotation_writer:
                annotation_writer.writerow(["sequence_id", "gene_name", "annotation", "source_gff"])

            for gff_file in gff_files:
                genome_name = os.path.basename(gff_file)
                if genome_name.endswith(name_split):
                    genome_name = genome_name[:-len(name_split)] if name_split else genome_name
                else:
                    compressed_suffix = str(name_split) + ".gz"
                    if genome_name.endswith(compressed_suffix):
                        genome_name = genome_name[:-len(compressed_suffix)]

                progress.set_current(genome_name)
                contig_chunks = collections.defaultdict(list)
                gff_features = []
                circular_contigs = set()
                fasta_section = False
                current_contig = None

                with _open_text(gff_file, "rt") as gff_handle:
                    seen_sequence_ids = collections.defaultdict(int)
                    for line_number, raw_line in enumerate(gff_handle, start=1):
                        line = raw_line.rstrip("\r\n")
                        if not line:
                            continue

                        if line == "##FASTA":
                            fasta_section = True
                            continue

                        if fasta_section:
                            if line.startswith(">"):
                                header = line[1:].strip()
                                if not header:
                                    raise ValueError(
                                        f"Empty embedded FASTA header. GFF={gff_file}; line={line_number}"
                                    )
                                current_contig = header.split(None, 1)[0]
                                if current_contig in contig_chunks:
                                    raise ValueError(
                                        f"Duplicate embedded FASTA contig {current_contig!r}. "
                                        f"GFF={gff_file}; line={line_number}"
                                    )
                                contig_chunks[current_contig] = []
                            else:
                                if current_contig is None:
                                    raise ValueError(
                                        f"Embedded FASTA sequence occurs before a header. "
                                        f"GFF={gff_file}; line={line_number}"
                                    )
                                contig_chunks[current_contig].append("".join(line.split()))
                            continue

                        if line.startswith("#"):
                            continue

                        fields = line.split("\t")
                        if len(fields) != 9:
                            raise ValueError(
                                f"Malformed GFF3 record: expected 9 tab-delimited fields, "
                                f"found {len(fields)}. GFF={gff_file}; line={line_number}"
                            )
                        if _gff_declares_circular(fields[8]):
                            circular_contigs.add(fields[0])
                        if not feature_matches_gene_ident(fields[2], gene_ident):
                            continue
                        try:
                            start = int(fields[3])
                            end = int(fields[4])
                        except ValueError as error:
                            raise ValueError(
                                f"Non-integer GFF coordinates. GFF={gff_file}; line={line_number}"
                            ) from error

                        try:
                            sequence_id = _extract_gff_identifier(fields[8])
                        except ValueError as error:
                            raise ValueError(
                                f"Could not identify GFF feature. GFF={gff_file}; "
                                f"line={line_number}: {error}"
                            ) from error

                        original_id = sequence_id
                        if original_id in seen_sequence_ids:
                            sequence_id = f"{original_id}_{seen_sequence_ids[original_id]}"
                            seen_sequence_ids[original_id] += 1
                        else:
                            seen_sequence_ids[original_id] = 1

                        gene_name, annotation = _gff_gene_metadata(fields[8])
                        gff_features.append(
                            (fields[0], start, end, fields[6], sequence_id, line_number,
                             gene_name, annotation)
                        )

                contig_sequences = {
                    contig: "".join(chunks)
                    for contig, chunks in contig_chunks.items()
                }

                if gff_features and not contig_sequences:
                    raise ValueError(
                        f"GFF contains selected gene features but no embedded FASTA section. "
                        f"GFF={gff_file}"
                    )

                for (contig, start, end, strand, sequence_id, line_number,
                     gene_name, annotation) in gff_features:
                    context = (
                        f"GFF={gff_file}; genome={genome_name}; gene={sequence_id}; "
                        f"contig={contig}; coordinates={start}-{end}; "
                        f"strand={strand}; GFF_line={line_number}"
                    )
                    sequence = _extract_gff_feature_sequence(
                        contig_sequences,
                        contig=contig,
                        start=start,
                        end=end,
                        strand=strand,
                        context=context,
                        circular_contigs=circular_contigs,
                    )

                    _write_gene_annotation(
                        annotation_writer, genome_name, sequence_id,
                        gene_name, annotation, gff_file,
                    )

                    if run_as_combiner:
                        output_sequence = translate_frame(sequence) if translate else sequence
                        dna_handle.write(
                            f">{genome_name}|{sequence_id}\n{wrap_sequence(output_sequence, 60)}\n"
                        )
                    else:
                        if translate:
                            aa_sequence = translate_frame(sequence)
                            aa_handle.write(
                                f">{genome_name}|{sequence_id}\n{wrap_sequence(aa_sequence, 60)}\n"
                            )
                        dna_handle.write(
                            f">{genome_name}|{sequence_id}\n{wrap_sequence(sequence, 60)}\n"
                        )
                    progress.add_records()

                progress.file_complete()
    finally:
        progress.close()

    if requested_aa_path and not translate and os.path.exists(requested_aa_path):
        try:
            os.remove(requested_aa_path)
        except OSError as error:
            logger.warning("Could not remove obsolete amino-acid output %s: %s", requested_aa_path, error)


def read_fasta_files(input_dir, name_split_fasta, combined_out, translate, run_as_combiner):
    """Combine per-genome FASTA files in deterministic file and record order.

    When ``translate=True`` the input is treated as nucleotide sequence and is
    strictly IUPAC-validated before translation. When ``translate=False`` the
    input alphabet is left unchanged because this legacy helper is also used for
    non-translated sequence input.
    """
    emit_pending_startup_messages()
    dna_output_path = os.fspath(combined_out)
    requested_aa_path = None
    if not run_as_combiner:
        requested_aa_path = dna_output_path.replace("_dna.fasta", "_aa.fasta")
        if requested_aa_path == dna_output_path:
            raise ValueError(
                "The combined output filename must contain '_dna.fasta' when "
                "run_as_combiner=False, so a separate '_aa.fasta' path can be made."
            )

    aa_output_path = requested_aa_path if translate else None
    fasta_files = _sorted_input_files(input_dir, name_split_fasta)
    if not fasta_files:
        raise FileNotFoundError(
            f"No FASTA files found in {input_dir} (pattern: *{name_split_fasta})."
        )

    logger.info(
        "Input FASTA files will be processed in deterministic basename order (%s files).",
        f"{len(fasta_files):,}",
    )

    with _atomic_text_outputs([dna_output_path, aa_output_path]) as outputs:
        dna_handle = outputs[dna_output_path]
        aa_handle = outputs[aa_output_path]

        for fasta_file in fasta_files:
            filename = os.path.basename(fasta_file)
            genome_name = filename[:-len(name_split_fasta)] if name_split_fasta else filename

            for header, sequence in _iter_fasta_records(fasta_file):
                sequence_id = header[1:].split(None, 1)[0]
                if translate:
                    validate_nucleotide_sequence(
                        sequence,
                        context=(
                            f"FASTA={fasta_file}; genome={genome_name}; "
                            f"record={header}"
                        ),
                    )
                    aa_sequence = translate_frame(sequence)
                else:
                    aa_sequence = None

                if run_as_combiner:
                    output_sequence = aa_sequence if translate else sequence
                    dna_handle.write(
                        f">{genome_name}|{sequence_id}\n{wrap_sequence(output_sequence, 60)}\n"
                    )
                else:
                    if translate:
                        aa_handle.write(
                            f">{genome_name}|{sequence_id}\n{wrap_sequence(aa_sequence, 60)}\n"
                        )
                    dna_handle.write(
                        f">{genome_name}|{sequence_id}\n{wrap_sequence(sequence, 60)}\n"
                    )

    if requested_aa_path and not translate and os.path.exists(requested_aa_path):
        try:
            os.remove(requested_aa_path)
        except OSError as error:
            logger.warning("Could not remove obsolete amino-acid output %s: %s", requested_aa_path, error)


def write_individual_groups(options, output_dir, key_order, cores, sequences,
                            pangenome_clusters_First_sequences_sorted,
                            combined_pangenome_clusters_Second_sequences):
    if not getattr(options, "write_individual_groups", False):
        return

    os.makedirs(output_dir, exist_ok=True)


    for key_prefix in key_order:
        for key, values in cores.items():
            if not key.startswith(key_prefix):
                continue

            for value in values:
                sequences_to_write = (pangenome_clusters_First_sequences_sorted[value]
                                      if 'First' in key_prefix
                                      else combined_pangenome_clusters_Second_sequences[value])

                dna_path = os.path.join(output_dir, f"{key}_{value}_dna.fasta")
                aa_path = dna_path.replace('_dna.fasta', '_aa.fasta')

                if getattr(options, "sequence_type", None) == 'AA':
                    with open(dna_path, 'w') as dna_f, open(aa_path, 'w') as aa_f:
                        for header in sequences_to_write:
                            if header not in sequences:
                                if getattr(options, "verbose", False):
                                    print(f"Sequence {header} not found in original_fasta file.")
                                continue
                            seq = sequences[header]
                            dna_f.write(f">{header}\n{wrap_sequence(seq)}\n")
                            aa_f.write(f">{header}\n{wrap_sequence(translate_frame(seq))}\n")
                else:
                    with open(dna_path, 'w') as dna_f:
                        for header in sequences_to_write:
                            if header not in sequences:
                                if getattr(options, "verbose", False):
                                    print(f"Sequence {header} not found in original_fasta file.")
                                continue
                            seq = sequences[header]
                            dna_f.write(f">{header}\n{wrap_sequence(seq)}\n")

def write_groups_func(options, output_dir, key_order, cores, sequences,
                 pangenome_clusters_First_sequences_sorted, combined_pangenome_clusters_Second_sequences):

    emit_pending_startup_messages()
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
     os.makedirs(output_dir)

    for group in options.write_groups.split(','):

        combined_fasta_filename = os.path.join(output_dir, "combined_group_sequences_" + group + "_dna.fasta")

        # Open combined FASTA file for writing all sequences
        with open(combined_fasta_filename, 'w') as combined_fasta, open(combined_fasta_filename.replace('_dna.fasta','_aa.fasta'), 'w') as combined_fasta_aa:
            for key_prefix in key_order:
                for key, values in cores.items():
                    if any(part in group for part in key.split('_')):
                        if key.startswith(key_prefix):
                            for value in values:
                                output_filename = f"{key}_{value}_dna.fasta"
                                if 'First' in key_prefix:
                                    sequences_to_write = pangenome_clusters_First_sequences_sorted[value]
                                else:
                                    sequences_to_write = combined_pangenome_clusters_Second_sequences[value]

                                # Write individual FASTA file
                                with open(os.path.join(output_dir,output_filename), 'w') as outfile, open(os.path.join(output_dir, output_filename.replace('_dna.fasta','_aa.fasta')), 'w') as outfile_aa:
                                    for header in sequences_to_write:
                                        if header in sequences:
                                            sequence = sequences[header]
                                            wrapped_sequence = wrap_sequence(sequence)
                                            # Handle Amino Acid Sequences (AA)
                                            if options.sequence_type == 'AA':
                                                seq_aa = translate_frame(sequence)
                                                wrapped_sequence_aa = wrap_sequence(seq_aa)
                                                # Write individual group file for AA, if option is enabled
                                                if options.write_individual_groups:
                                                    outfile_aa.write(f">{header}\n")
                                                    outfile_aa.write(f"{wrapped_sequence_aa}\n")
                                                else:
                                                    try:
                                                        os.remove(outfile_aa.name)  # Delete individual file if option is disabled
                                                    except FileNotFoundError:
                                                        pass
                                                # Always write to the combined AA file
                                                combined_fasta_aa.write(f">Group_{value}|{header}\n")
                                                combined_fasta_aa.write(f"{wrapped_sequence_aa}\n")
                                            # Handle Nucleotide Sequences
                                            else:
                                                # If the option is disabled, delete individual AA file (if created)
                                                try:
                                                    os.remove(outfile_aa.name)  # Ensure outfile_aa is removed when sequence_type isn't 'AA'
                                                except FileNotFoundError:
                                                    pass
                                            # Write individual group file for nucleotide sequence, if option is enabled
                                            if options.write_individual_groups:
                                                outfile.write(f">{header}\n")
                                                outfile.write(f"{wrapped_sequence}\n")
                                            else:
                                                try:
                                                    os.remove(outfile.name)  # Delete individual file if option is disabled
                                                except FileNotFoundError:
                                                    pass
                                            # Always write to the combined nucleotide file
                                            combined_fasta.write(f">Group_{value}|{header}\n")
                                            combined_fasta.write(f"{wrapped_sequence}\n")

                                        else:
                                            if options.verbose == True:
                                                logger.info("Sequence " + header + " not found in original_fasta file.")
        if options.sequence_type != 'AA':
            #Clean up unused file
            try:
                os.remove(combined_fasta_aa.name)
            except FileNotFoundError:
                pass
    logger.info("Combined FASTA file saved to: " + combined_fasta_filename)


# def process_gene_groups(options, group_directory, sub_group_directory, paralog_groups, output_file):
#     """Process each gene family file to select the longest sequence per genome and concatenate aligned sequences."""
#     concatenated_sequences = {}
#     output_file = group_directory.replace('Gene_Groups_Output',output_file)
#
#     # Iterate over each gene family file
#     for gene_file in os.listdir(group_directory):
#         if gene_file.endswith('.fasta') and not gene_file.endswith('combined_group_sequences.fasta') :
#             gene_path = os.path.join(group_directory, gene_file)
#
#             # Read sequences from the gene family file
#             sequences = read_fasta(gene_path)
#
#             # Select the longest sequence for each genome
#             longest_sequences = select_longest_gene(sequences)
#
#             # Run mafft on the longest sequences
#             aligned_file = f"{group_directory}/{gene_file}_aligned.fasta.tmp"
#             run_mafft_on_sequences(options, {seq_id: seq for seq_id, seq in longest_sequences.values()}, aligned_file)
#
#             # Read aligned sequences and concatenate them
#             aligned_sequences = read_fasta(aligned_file)
#             for genome, aligned_seq in aligned_sequences.items():
#                 genome_name = genome.split('|')[0]
#                 if genome_name not in concatenated_sequences:
#                     concatenated_sequences[genome_name] = ""
#                 concatenated_sequences[genome_name] += aligned_seq
#
#             # Clean up aligned file
#             os.remove(aligned_file)
#
#     # Write the concatenated sequences to the output file
#     with open(output_file, 'w') as out:
#         for genome, sequence in concatenated_sequences.items():
#             out.write(f">{genome}\n")
#             wrapped_sequence = wrap_sequence(sequence, 60)
#             out.write(f"{wrapped_sequence}\n")

def perform_alignment(gene_path,group_directory, gene_file, options, concatenated_sequences, subgrouped):
    # Read sequences from the gene family file
    sequences = read_fasta(gene_path)
    # Select the longest sequence for each genome
    longest_sequences = select_longest_gene(sequences, subgrouped)

    aligned_file = f"{group_directory}/{gene_file}_aligned.fasta.tmp"
    if len(longest_sequences) > 1:
        run_mafft_on_sequences(
            options,
            {seq_id: seq for seq_id, seq in longest_sequences.values()},
            aligned_file,
        )
        aligned_sequences = read_fasta(aligned_file)
    else:
        aligned_sequences = {
            seq_id: sequence for seq_id, sequence in longest_sequences.values()
        }
    # Find the length of the longest sequence in aligned_sequences
    max_length = max(len(seq) for seq in aligned_sequences.values())

    aligned_by_genome = {}
    for sequence_id, sequence in aligned_sequences.items():
        parts = sequence_id.split('|')
        aligned_genome = parts[1] if subgrouped and len(parts) > 1 else parts[0]
        if aligned_genome in aligned_by_genome:
            raise ValueError(f"Alignment contains multiple representatives for genome {aligned_genome}")
        aligned_by_genome[aligned_genome] = sequence
    for genome in concatenated_sequences:
        concatenated_sequences[genome] += aligned_by_genome.get(genome, "-" * max_length)

    # Clean up aligned file
    if os.path.exists(aligned_file):
        os.remove(aligned_file)

    return concatenated_sequences

def process_gene_groups(options, group_directory, sub_group_directory, paralog_groups, genome_list, output_file):
    emit_pending_startup_messages()
    """Process each gene family file to select the longest sequence per genome and concatenate aligned sequences."""
    concatenated_sequences = {genome: "" for genome in genome_list}
    output_file = group_directory.replace('Gene_Groups_Output', output_file)
    threshold_percent = int(
        getattr(
            options,
            'align_core_threshold',
            99 if getattr(options, 'align_core', False) is True else options.align_core,
        )
    )
    threshold_size = math.ceil(threshold_percent * int(options.genome_num) / 100)

    if getattr(options, 'align_aa', False) == True:
     affix = '_aa.fasta'
    else:
     affix = '_dna.fasta'

    if getattr(options, 'align_core', False):
     for gene_file in sorted(os.listdir(group_directory)):
         if not gene_file.endswith(affix) or gene_file.startswith('combined_group_sequences'):
             continue
         match = re.match(
             r'^(?P<category>.+_core_(?P<percent>\d+))_(?P<family>group_[0-9a-f]+|\d+)_(?:dna|aa)\.fasta$',
             gene_file,
         )
         if not match or int(match.group('percent')) < threshold_percent:
             continue
         family_id = match.group('family')
         gene_path = os.path.join(group_directory, gene_file)
         paralog_key = next(
             (
                 key for key in (f'>{family_id}', f'>Group_{family_id}')
                 if paralog_groups and key in paralog_groups
             ),
             None,
         )
         if sub_group_directory and paralog_key:
             for subgroup, size in enumerate(paralog_groups[paralog_key]['sizes']):
                 if size < threshold_size:
                     continue
                 subgroup_path = os.path.join(
                     sub_group_directory,
                     f"{family_id}_subgroup_{subgroup}{affix}",
                 )
                 concatenated_sequences = perform_alignment(
                     subgroup_path,
                     group_directory,
                     gene_file,
                     options,
                     concatenated_sequences,
                     True,
                 )
         else:
             concatenated_sequences = perform_alignment(
                 gene_path,
                 group_directory,
                 gene_file,
                 options,
                 concatenated_sequences,
                 False,
             )

     # Write the concatenated sequences to the output file
     with open(output_file, 'w') as out:
         for genome, sequence in concatenated_sequences.items():
             out.write(f">{genome}\n")
             wrapped_sequence = wrap_sequence(sequence, 60)
             out.write(f"{wrapped_sequence}\n")