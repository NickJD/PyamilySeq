"""Preparation helpers for representative-based incremental PyamilySeq runs.

This module deliberately streams FASTA input one record at a time. The large
original gene FASTA is not materialised as a dictionary in RAM.
"""
from __future__ import annotations

import gzip
import logging
import os
import re
from collections import OrderedDict
from typing import Iterator

try:
    from .utils import detect_sequence_type, translate_frame
except ImportError:
    from utils import detect_sequence_type, translate_frame


class IncrementalInputError(ValueError):
    """Raised when an incremental run's FASTA/cluster inputs are inconsistent."""


def _open_text(path: str):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    return opener(path, "rt", encoding="utf-8")


def _iter_fasta_records(path: str) -> Iterator[tuple[str, str]]:
    """Yield (first-token identifier, sequence) from a FASTA, with gzip support."""
    header = None
    sequence_parts: list[str] = []
    with _open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)
                header = line[1:].strip().split(maxsplit=1)[0] if line[1:].strip() else ""
                if not header:
                    raise IncrementalInputError(
                        f"Empty FASTA identifier at line {line_number} in {path}"
                    )
                sequence_parts = []
            else:
                if header is None:
                    raise IncrementalInputError(
                        f"Sequence data before first FASTA header at line {line_number} in {path}"
                    )
                sequence_parts.append("".join(line.split()))
    if header is not None:
        yield header, "".join(sequence_parts)


def _write_fasta_record(handle, identifier: str, sequence: str, width: int = 60):
    handle.write(f">{identifier}\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start:start + width] + "\n")
    if not sequence:
        handle.write("\n")


def _read_cdhit_representatives(cluster_file: str) -> tuple[list[str], int]:
    """Read one representative identifier per CD-HIT cluster.

    This parser mirrors the native CD-HIT ``.clstr`` format and works with
    compressed ``.clstr.gz`` files too.
    """
    member_re = re.compile(r"^\s*\d+\s+[^,]+,\s+>(.*?)\.\.\.\s*(.*)$")
    representatives: list[str] = []
    seen = set()
    current_cluster = None
    current_rep = None
    cluster_count = 0

    def finish_cluster(label, representative):
        if label is None:
            return
        if representative is None:
            raise IncrementalInputError(
                f"CD-HIT cluster {label!r} in {cluster_file} has no representative ('*') entry"
            )
        if representative in seen:
            raise IncrementalInputError(
                f"Representative {representative!r} occurs in more than one cluster in {cluster_file}"
            )
        seen.add(representative)
        representatives.append(representative)

    with _open_text(cluster_file) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                finish_cluster(current_cluster, current_rep)
                match = re.match(r"^>Cluster\s+(.+?)\s*$", line, flags=re.IGNORECASE)
                if not match:
                    raise IncrementalInputError(
                        f"Expected a CD-HIT '>Cluster N' header at line {line_number} in {cluster_file}; got {line!r}"
                    )
                current_cluster = match.group(1)
                current_rep = None
                cluster_count += 1
                continue
            if current_cluster is None:
                raise IncrementalInputError(
                    f"CD-HIT member before a cluster header at line {line_number} in {cluster_file}"
                )
            match = member_re.match(line)
            if not match:
                raise IncrementalInputError(
                    f"Malformed CD-HIT member line {line_number} in {cluster_file}: {line!r}"
                )
            identifier, member_detail = match.group(1), match.group(2)
            if "*" in member_detail:
                if current_rep is not None and current_rep != identifier:
                    raise IncrementalInputError(
                        f"CD-HIT cluster {current_cluster!r} has multiple representatives"
                    )
                current_rep = identifier
    finish_cluster(current_cluster, current_rep)
    if not representatives:
        raise IncrementalInputError(f"No CD-HIT clusters/representatives found in {cluster_file}")
    if len(representatives) != cluster_count:
        raise IncrementalInputError(
            f"Parsed {len(representatives)} representatives for {cluster_count} clusters in {cluster_file}"
        )
    return representatives, cluster_count


def _normalise_new_fasta(
    source_fasta: str,
    output_fasta: str,
    sequence_delimiter: str,
    generated_ids_use_pipe: bool = False,
):
    """Write a canonical, uncompressed new-gene FASTA and collect new IDs.

    The existing GFF extractors build identifiers with ``|``. If a caller has
    configured a different delimiter, translate only the first (genome) split
    in those generated IDs so they stay consistent with the original FASTA.
    """
    identifiers = set()
    genomes = set()
    count = 0
    os.makedirs(os.path.dirname(os.path.abspath(output_fasta)), exist_ok=True)
    with open(output_fasta, "w", encoding="utf-8") as out:
        for identifier, sequence in _iter_fasta_records(source_fasta):
            if generated_ids_use_pipe and sequence_delimiter and sequence_delimiter != "|":
                if "|" not in identifier:
                    raise IncrementalInputError(
                        f"Generated GFF identifier {identifier!r} lacks the expected '|' delimiter"
                    )
                genome, gene = identifier.split("|", 1)
                identifier = genome + sequence_delimiter + gene
            if identifier in identifiers:
                raise IncrementalInputError(
                    f"Duplicate sequence identifier {identifier!r} in new-gene FASTA {source_fasta}"
                )
            if sequence_delimiter and sequence_delimiter not in identifier:
                raise IncrementalInputError(
                    f"New-gene identifier {identifier!r} does not contain the configured sequence delimiter "
                    f"{sequence_delimiter!r}. IDs must identify the genome, e.g. genome|gene."
                )
            if not sequence:
                raise IncrementalInputError(f"Empty sequence for {identifier!r} in {source_fasta}")
            identifiers.add(identifier)
            genome = identifier.split(sequence_delimiter, 1)[0] if sequence_delimiter else identifier
            genomes.add(genome)
            _write_fasta_record(out, identifier, sequence)
            count += 1
    if count == 0:
        raise IncrementalInputError(f"No FASTA records found in new-gene input {source_fasta}")
    return identifiers, sorted(genomes), count


def prepare_incremental_inputs(
    cluster_file: str,
    original_fasta: str,
    new_genes_fasta: str,
    output_dir: str,
    sequence_delimiter: str = "|",
    logger=None,
    generated_ids_use_pipe: bool = False,
) -> dict:
    """Create full sequence FASTA and representative-plus-new FASTA for CD-HIT.

    ``original_fasta`` must be the full DNA or protein FASTA corresponding to
    ``cluster_file``. Its alphabet is detected automatically. ``new_genes_fasta``
    must contain DNA sequences whose IDs already follow the configured
    genome/gene identifier convention.
    """
    logger = logger or logging.getLogger("PyamilySeq.PyamilySeq")
    for path, label in ((cluster_file, "original CD-HIT cluster file"),
                        (original_fasta, "original full DNA/protein FASTA"),
                        (new_genes_fasta, "new-gene DNA FASTA")):
        if not path or not os.path.isfile(path):
            raise IncrementalInputError(f"The {label} does not exist: {path}")

    os.makedirs(output_dir, exist_ok=True)
    representatives, original_cluster_count = _read_cdhit_representatives(cluster_file)
    # The original run may have clustered either nucleotide or protein
    # sequences. Keep the reconciliation FASTA and reclustering FASTA in that
    # same alphabet; only translate the newly extracted nucleotide genes when
    # the original input was protein.
    original_is_dna = detect_sequence_type(original_fasta)
    new_is_dna = detect_sequence_type(new_genes_fasta)
    if not original_is_dna and not new_is_dna:
        translate_new_genes = False
    elif not original_is_dna and new_is_dna:
        translate_new_genes = True
    elif original_is_dna and not new_is_dna:
        raise IncrementalInputError(
            "The original FASTA contains DNA but the new-gene FASTA contains protein sequences; "
            "incremental clustering requires matching sequence alphabets."
        )
    else:
        translate_new_genes = False
    representative_set = set(representatives)
    logger.info(
        "Loaded %s original CD-HIT representatives from %s clusters",
        f"{len(representatives):,}", f"{original_cluster_count:,}",
    )

    normalised_new_fasta = os.path.join(output_dir, "incremental_new_genes_dna.fasta")
    new_ids, new_genomes, new_count = _normalise_new_fasta(
        new_genes_fasta,
        normalised_new_fasta,
        sequence_delimiter,
        generated_ids_use_pipe=generated_ids_use_pipe,
    )
    if translate_new_genes:
        translated_new_fasta = os.path.join(output_dir, "incremental_new_genes_aa.fasta")
        with open(translated_new_fasta, "w", encoding="utf-8") as out:
            for identifier, sequence in _iter_fasta_records(normalised_new_fasta):
                _write_fasta_record(out, identifier, translate_frame(sequence))
        normalised_new_fasta = translated_new_fasta
    collisions = representative_set.intersection(new_ids)
    if collisions:
        examples = ", ".join(sorted(collisions)[:10])
        raise IncrementalInputError(
            f"New-gene FASTA reuses {len(collisions):,} identifier(s) already used by original representatives: {examples}"
        )

    alphabet = "dna" if original_is_dna else "aa"
    combined_path = os.path.join(output_dir, f"incremental_combined_sequences_{alphabet}.fasta")
    reps_path = os.path.join(output_dir, f"incremental_representatives_{alphabet}.fasta")
    protected_inputs = {os.path.realpath(cluster_file), os.path.realpath(original_fasta),
                        os.path.realpath(new_genes_fasta)}
    for output_path in (combined_path, reps_path, normalised_new_fasta):
        if os.path.realpath(output_path) in protected_inputs:
            raise IncrementalInputError(
                f"Output path would overwrite one of the supplied input files: {output_path}. "
                "Choose a different output directory or input filename."
            )
    combined_tmp = combined_path + ".tmp"
    reps_tmp = reps_path + ".tmp"
    found_representatives = set()
    original_sequence_count = 0
    try:
        with open(combined_tmp, "w", encoding="utf-8") as combined_out, \
             open(reps_tmp, "w", encoding="utf-8") as reps_out:
            for identifier, sequence in _iter_fasta_records(original_fasta):
                if identifier in new_ids:
                    raise IncrementalInputError(
                        f"New sequence identifier {identifier!r} already exists in original FASTA {original_fasta}"
                    )
                if not sequence:
                    raise IncrementalInputError(
                        f"Empty original sequence for {identifier!r} in {original_fasta}"
                    )
                _write_fasta_record(combined_out, identifier, sequence)
                original_sequence_count += 1
                if identifier in representative_set:
                    if identifier in found_representatives:
                        raise IncrementalInputError(
                            f"Representative identifier {identifier!r} occurs more than once in {original_fasta}"
                        )
                    _write_fasta_record(reps_out, identifier, sequence)
                    found_representatives.add(identifier)
            if original_sequence_count == 0:
                raise IncrementalInputError(f"No FASTA records found in original FASTA {original_fasta}")
            missing_representatives = representative_set - found_representatives
            if missing_representatives:
                examples = ", ".join(sorted(missing_representatives)[:10])
                raise IncrementalInputError(
                    f"{len(missing_representatives):,} representative(s) from {cluster_file} were absent from "
                    f"{original_fasta}: {examples}. The cluster file and original FASTA must correspond."
                )
            # Append each new gene to both: the full all-gene FASTA and the smaller
            # FASTA that CD-HIT will recluster.
            with open(normalised_new_fasta, "r", encoding="utf-8") as new_in:
                for line in new_in:
                    combined_out.write(line)
                    reps_out.write(line)
        os.replace(combined_tmp, combined_path)
        os.replace(reps_tmp, reps_path)
    except Exception:
        for temporary in (combined_tmp, reps_tmp):
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass
        raise

    logger.info("Original sequences retained: %s", f"{original_sequence_count:,}")
    logger.info("Original representatives retained: %s", f"{len(found_representatives):,}")
    logger.info("New genes added to incremental clustering input: %s", f"{new_count:,}")
    logger.info("New genomes detected: %s", f"{len(new_genomes):,}")
    return {
        "combined_fasta": combined_path,
        "representatives_fasta": reps_path,
        "new_genes_fasta": normalised_new_fasta,
        "representatives_dna_fasta": reps_path,
        "representatives_aa_fasta": reps_path,
        "new_genes_dna_fasta": normalised_new_fasta,
        "original_sequence_count": original_sequence_count,
        "original_cluster_count": original_cluster_count,
        "representative_count": len(found_representatives),
        "new_gene_count": new_count,
        "new_genomes": new_genomes,
        "original_is_dna": original_is_dna,
    }
