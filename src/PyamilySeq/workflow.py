"""Shared Species/Genus pangenome workflow built on the canonical family ledger."""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import subprocess
import logging
import sqlite3
import tempfile
from collections import Counter, OrderedDict

try:
    from .families import FamilyInputError, build_family_ledger, read_fasta_records
    from .utils import detect_sequence_type, translate_frame, wrap_sequence
    from .visualisations import write_ledger_visualisations

except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from families import FamilyInputError, build_family_ledger, read_fasta_records
    from utils import detect_sequence_type, translate_frame, wrap_sequence
    from visualisations import write_ledger_visualisations



GPA_METADATA_COLUMNS = [
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
]


def _normalise_species_groups(value):
    try:
        groups = [int(item.strip()) for item in str(value).split(",") if item.strip()]
    except ValueError:
        raise FamilyInputError("Species groups must be comma-separated integer percentages") from None
    if not groups:
        raise FamilyInputError("At least one species group is required")
    if any(group < 0 or group > 100 for group in groups):
        raise FamilyInputError("Species groups must be between 0 and 100")
    groups = sorted(set(groups), reverse=True)
    if 0 not in groups:
        groups.append(0)
    return groups


def _normalise_genus_groups(value):
    raw = [item.strip() for item in str(value).split(",") if item.strip() and item.strip() != ">"]
    try:
        groups = set()
        for item in raw:
            if "-" in item:
                start_text, end_text = item.split("-", 1)
                start, end = int(start_text), int(end_text)
                if start > end:
                    raise ValueError
                groups.update(range(start, end + 1))
            else:
                groups.add(int(item))
        groups = sorted(groups)
    except ValueError:
        raise FamilyInputError(
            "Genus groups must be positive integers or inclusive ranges such as 1-10"
        ) from None
    if not groups or any(group < 1 for group in groups):
        raise FamilyInputError("Genus groups must contain positive integers")
    return groups


def _prevalence_bin(count, taxa_total, mode, groups):
    if mode == "Species":
        for group in groups:
            threshold = 0 if group == 0 else math.ceil(group * taxa_total / 100.0)
            if count >= threshold:
                return str(group)
        return "0"
    for group in groups:
        if count <= group:
            return str(group)
    return ">"


def _prevalence_bin_label(group, groups, mode):
    """Return a human-readable label for one mutually exclusive prevalence bin."""
    if mode == "Species":
        value = int(group)
        positive = sorted({int(item) for item in groups if int(item) > 0}, reverse=True)
        if value == 0:
            return f"<{min(positive)}%" if positive else "0% threshold"
        higher = [item for item in positive if item > value]
        if not higher:
            return f">={value}%"
        return f"{value}-<{min(higher)}%"

    numeric = sorted({int(item) for item in groups})
    if str(group) == ">":
        return f">{numeric[-1]} genera" if numeric else "above configured genus bins"
    value = int(group)
    lower_values = [item for item in numeric if item < value]
    lower = (max(lower_values) + 1) if lower_values else 1
    return f"<={value} genera" if lower == 1 else f"{lower}-{value} genera"


def _stage_prefix(stage, mode):
    suffix = "core" if mode == "Species" else "genera"
    names = {
        "First": "First",
        "extended": "extended",
        "combined": "combined",
        "only_Second": "only_Second",
    }
    return f"{names.get(stage, stage)}_{suffix}"

# Human-readable labels for reconciliation provenance.
# These labels do not describe biological core/accessory status.
RECONCILIATION_STAGE_LABELS = {
    "First": "Original-only family",
    "extended": "Original family + supplement",
    "combined": "Merged original families",
    "only_Second": "Supplement-only family",
}


def _reconciliation_stage_label(stage):
    """Return a readable label for an internal reconciliation stage."""
    return RECONCILIATION_STAGE_LABELS.get(
        str(stage or ""),
        str(stage or "Unknown")
    )

def classify_families(ledger, mode, group_text):
    groups = (
        _normalise_species_groups(group_text)
        if mode == "Species"
        else _normalise_genus_groups(group_text)
    )
    taxa_total = len(ledger.taxa)
    categories = OrderedDict()
    stage_order = ("First", "extended", "combined", "only_Second")
    bin_order = [str(group) for group in groups]
    if mode == "Genus":
        bin_order.append(">")
    for stage in stage_order:
        for group in bin_order:
            categories[f"{_stage_prefix(stage, mode)}_{group}"] = []
    for family in ledger.families.values():
        group = _prevalence_bin(family.prevalence, taxa_total, mode, groups)
        family.category = f"{_stage_prefix(family.stage, mode)}_{group}"
        categories.setdefault(family.category, []).append(family.family_id)
    return categories, groups


class _GeneAnnotationStore:
    """Disk-backed lookup for per-gene GFF annotations during GPA writing.

    The sidecar can contain millions of rows. Keeping it in a Python dictionary
    would add substantial memory pressure to a workflow that otherwise streams
    sequence storage through SQLite.
    """

    def __init__(self, annotation_files=None, temp_dir=None):
        if annotation_files is None:
            annotation_files = []
        elif isinstance(annotation_files, (str, os.PathLike)):
            annotation_files = [os.fspath(annotation_files)]
        else:
            annotation_files = [os.fspath(item) for item in annotation_files if item]

        self.annotation_files = list(dict.fromkeys(annotation_files))
        missing = [path for path in self.annotation_files if not os.path.isfile(path)]
        if missing:
            raise FamilyInputError(
                "Gene-annotation sidecar file(s) not found: " + ", ".join(missing)
            )
        if not self.annotation_files:
            self.path = None
            self.connection = None
            self.lookup_cursor = None
            self.count = 0
            return

        handle = tempfile.NamedTemporaryFile(
            prefix="pyamilyseq_gene_annotations_", suffix=".sqlite3",
            dir=temp_dir, delete=False
        )
        self.path = handle.name
        handle.close()
        self.connection = sqlite3.connect(self.path)
        self.connection.execute(
            "CREATE TABLE annotations ("
            "sequence_id TEXT PRIMARY KEY, gene_name TEXT NOT NULL, annotation TEXT NOT NULL)"
        )
        self.count = 0
        self.lookup_cursor = None
        try:
            self._load()
            self.lookup_cursor = self.connection.cursor()
        except sqlite3.IntegrityError as error:
            self.close()
            raise FamilyInputError(
                "Duplicate sequence identifiers were found in gene-annotation sidecar files. "
                "Check for overlapping or duplicated GFF inputs."
            ) from error
        except Exception:
            self.close()
            raise

    def _load(self):
        batch = []
        cursor = self.connection.cursor()
        for path in self.annotation_files:
            opener = gzip.open if str(path).lower().endswith(".gz") else open
            with opener(path, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                required = {"sequence_id", "gene_name", "annotation"}
                if not reader.fieldnames or not required.issubset(set(reader.fieldnames)):
                    raise FamilyInputError(
                        f"Gene-annotation file {path} must contain columns: "
                        "sequence_id, gene_name, annotation."
                    )
                for line_number, row in enumerate(reader, start=2):
                    identifier = (row.get("sequence_id") or "").strip()
                    if not identifier:
                        raise FamilyInputError(
                            f"Empty sequence_id in gene-annotation file {path} at line {line_number}."
                        )
                    gene_name = (row.get("gene_name") or "").strip()
                    annotation = (row.get("annotation") or "").strip()
                    existing = cursor.execute(
                        "SELECT gene_name, annotation FROM annotations WHERE sequence_id = ?",
                        (identifier,),
                    ).fetchone()
                    if existing:
                        if existing != (gene_name, annotation):
                            raise FamilyInputError(
                                f"Conflicting GFF annotations for sequence {identifier!r} "
                                f"in {path}. Check for duplicate genome/gene identifiers."
                            )
                        continue
                    batch.append((identifier, gene_name, annotation))
                    if len(batch) >= 10000:
                        cursor.executemany(
                            "INSERT INTO annotations VALUES (?, ?, ?)", batch
                        )
                        self.count += len(batch)
                        batch.clear()
        if batch:
            cursor.executemany("INSERT INTO annotations VALUES (?, ?, ?)", batch)
            self.count += len(batch)
        try:
            self.connection.commit()
        except sqlite3.IntegrityError as error:
            raise FamilyInputError(
                "Duplicate sequence identifiers were found in gene-annotation sidecar files. "
                "Check for overlapping or duplicated GFF inputs."
            ) from error

    def get(self, sequence_id):
        if self.lookup_cursor is None:
            return "", ""
        try:
            self.lookup_cursor.execute(
                "SELECT gene_name, annotation FROM annotations WHERE sequence_id = ?",
                (sequence_id,),
            )
            return self.lookup_cursor.fetchone() or ("", "")
        except sqlite3.Error as error:
            raise FamilyInputError(
                f"Could not retrieve GFF annotation for sequence {sequence_id!r}: {error}"
            ) from error

    def close(self):
        connection = getattr(self, "connection", None)
        if connection is not None:
            connection.close()
            self.connection = None
        path = getattr(self, "path", None)
        if path:
            for suffix in ("", "-journal", "-wal", "-shm"):
                candidate = path + suffix
                try:
                    os.remove(candidate)
                except FileNotFoundError:
                    pass
                except OSError:
                    logging.getLogger("PyamilySeq.PyamilySeq").warning(
                        "Could not remove temporary gene-annotation database %s", candidate
                    )

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


def _family_metadata_consensus(values):
    """Return the most frequent non-empty value, breaking ties lexically."""
    counts = Counter(value for value in values if value)
    if not counts:
        return ""
    maximum = max(counts.values())
    return sorted(value for value, count in counts.items() if count == maximum)[0]


def _write_gpa(path, ledger, family_ids, mode, gene_annotation_files=None, temp_dir=None):
    columns = list(ledger.genomes if mode == "Species" else ledger.taxa)
    column_resolver = ledger.resolver.genome if mode == "Species" else ledger.resolver.taxon
    annotation_store = _GeneAnnotationStore(gene_annotation_files, temp_dir=temp_dir)
    try:
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(GPA_METADATA_COLUMNS + columns)
            for family in ledger.sorted_families(family_ids):
                occupied = family.prevalence
                sequence_count = len(family.members)
                cells = {column: [] for column in columns}
                lengths = []
                gene_names = []
                annotations = []
                for member in family.members:
                    cell_value = member
                    delimiter = ledger.resolver.sequence_delimiter
                    if mode == "Species" and delimiter and delimiter in member:
                        cell_value = member.split(delimiter, 1)[1]
                    cells[column_resolver(member)].append(cell_value)
                    lengths.append(ledger.sequences.sequence_length(member))
                    gene_name, annotation = annotation_store.get(member)
                    if gene_name:
                        gene_names.append(gene_name)
                    if annotation:
                        annotations.append(annotation)
                metadata = [
                    family.family_id,
                    _family_metadata_consensus(gene_names),
                    _family_metadata_consensus(annotations),
                    occupied,
                    sequence_count,
                    f"{sequence_count / occupied:.6f}" if occupied else "0.000000",
                    "",
                    "",
                    "",
                    "",
                    "",
                    min(lengths) if lengths else "",
                    max(lengths) if lengths else "",
                    f"{sum(lengths) / len(lengths):.2f}" if lengths else "",
                ]
                writer.writerow(metadata + ["  ".join(cells[column]) for column in columns])
    finally:
        annotation_store.close()


def _write_gene_group_table(path, ledger, family_ids):
    selected = set(family_ids)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for member in ledger.sequences:
            family_id = ledger.assignment[member]
            if family_id in selected:
                writer.writerow([member, family_id])


def _write_core_group_table(path, ledger, family_ids):
    selected = set(family_ids)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for member in ledger.sequences:
            family = ledger.families[ledger.assignment[member]]
            if family.family_id in selected:
                writer.writerow([member, family.category])


def _write_family_summary(path, ledger, family_ids, mode="Species", groups=None):
    """Write family-level data with stage provenance and prevalence kept separate."""
    if groups is None:
        # Backwards-compatible fallback for external callers. Normal workflow
        # passes the exact groups returned by classify_families().
        raw_groups = {
            family.category.rsplit("_", 1)[-1]
            for family in ledger.families.values()
            if getattr(family, "category", "")
        }
        if mode == "Species":
            groups = sorted((int(item) for item in raw_groups), reverse=True)
        else:
            groups = sorted((int(item) for item in raw_groups if item != ">"))
            if ">" in raw_groups:
                groups.append(">")

    prevalence_denominator = len(ledger.taxa)
    fields = [
        "family_id",
        "stage",  # Original internal identifier
        "reconciliation_stage_label",  # Human-readable provenance
        "prevalence_bin",
        "prevalence_count",
        "prevalence_denominator",
        "prevalence_percent",
        "category",  # Retain for downstream compatibility
        "representative",
        "taxa_count",
        "genome_count",
        "sequence_count",
        "mean_copies_per_occupied_genome",
        "multicopy_genomes",
        "source_family_count",
        "source_families",
        "source_clusters",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for family in ledger.sorted_families(family_ids):
            # Use the same member-by-member genome resolution as the
            # single-copy filter, including all First + Second members.
            genome_counts = ledger.family_copy_number(family.family_id)
            raw_bin = family.category.rsplit("_", 1)[-1]
            prevalence_percent = (
                100.0 * family.prevalence / prevalence_denominator
                if prevalence_denominator else 0.0
            )
            writer.writerow({
                "family_id": family.family_id,
                "stage": family.stage,
                "reconciliation_stage_label": _reconciliation_stage_label(
                    family.stage
                ),
                "prevalence_bin": _prevalence_bin_label(
                    raw_bin, groups, mode
                ),
                "prevalence_count": family.prevalence,
                "prevalence_denominator": prevalence_denominator,
                "prevalence_percent": f"{prevalence_percent:.4f}",
                # Retain the legacy compound category for downstream compatibility.
                "category": family.category,
                "representative": family.representative,
                "taxa_count": family.prevalence,
                "genome_count": len(genome_counts),
                "sequence_count": len(family.members),
                "mean_copies_per_occupied_genome": (
                    f"{len(family.members) / len(genome_counts):.6f}" if genome_counts else "0.000000"
                ),
                "multicopy_genomes": sum(count > 1 for count in genome_counts.values()),
                "source_family_count": len(family.source_families),
                "source_families": ";".join(family.source_families),
                "source_clusters": ";".join(family.source_clusters),
            })


def _write_reference(path, ledger, family_ids):
    with open(path, "w", encoding="utf-8") as handle:
        for family in ledger.sorted_families(family_ids):
            representative = family.representative
            handle.write(
                f">{family.family_id} representative={representative} category={family.category} "
                f"members={len(family.members)}\n"
            )
            handle.write(wrap_sequence(ledger.sequences[representative], 60) + "\n")


def _requested_family_ids(options, ledger, selected_ids):
    requested = getattr(options, "write_groups", None)
    if not requested:
        return list(selected_ids)
    tokens = {token.strip() for token in str(requested).split(",") if token.strip()}
    return [
        family_id for family_id in selected_ids
        if ledger.families[family_id].category.rsplit("_", 1)[-1] in tokens
        or ledger.families[family_id].category in tokens
    ]


def _write_family_fastas(options, ledger, family_ids):
    output_dir = os.path.join(options.output_dir, "Gene_Groups_Output")
    os.makedirs(output_dir, exist_ok=True)
    write_individual = bool(getattr(options, "write_individual_groups", False))
    sequence_type = getattr(options, "sequence_type", None)
    input_is_dna = detect_sequence_type(options.fasta)
    requested = getattr(options, "write_groups", None) or "all"
    safe_requested = "_".join(token.strip() for token in str(requested).split(",") if token.strip())
    alphabet = "dna" if input_is_dna else "aa"
    combined_path = os.path.join(output_dir, f"combined_group_sequences_{safe_requested}_{alphabet}.fasta")
    translated_path = combined_path.replace("_dna.fasta", "_aa.fasta")
    with open(combined_path, "w", encoding="utf-8") as combined:
        translated = (
            open(translated_path, "w", encoding="utf-8")
            if input_is_dna and sequence_type == "AA"
            else None
        )
        try:
            for family in ledger.sorted_families(family_ids):
                individual = None
                individual_translated = None
                if write_individual:
                    base = f"{family.category}_{family.family_id}"
                    individual = open(
                        os.path.join(output_dir, f"{base}_{alphabet}.fasta"), "w", encoding="utf-8"
                    )
                    if translated:
                        individual_translated = open(
                            os.path.join(output_dir, base + "_aa.fasta"), "w", encoding="utf-8"
                        )
                try:
                    for member in family.members:
                        sequence = ledger.sequences[member]
                        header = f"{family.family_id}|{member}"
                        combined.write(f">{header}\n{wrap_sequence(sequence, 60)}\n")
                        if individual:
                            individual.write(f">{member}\n{wrap_sequence(sequence, 60)}\n")
                        if translated:
                            amino_acids = translate_frame(sequence)
                            translated.write(f">{header}\n{wrap_sequence(amino_acids, 60)}\n")
                            if individual_translated:
                                individual_translated.write(
                                    f">{member}\n{wrap_sequence(amino_acids, 60)}\n"
                                )
                finally:
                    if individual:
                        individual.close()
                    if individual_translated:
                        individual_translated.close()
        finally:
            if translated:
                translated.close()
    return output_dir


def _run_mafft(options, input_path, output_path):
    command = ["mafft", "--thread", str(options.threads), "--auto", input_path]
    with open(output_path, "w", encoding="utf-8") as output:
        result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"MAFFT failed for {input_path} with exit code {result.returncode}: {result.stderr.strip()}"
        )


def _write_alignment(options, ledger, family_ids):
    """Align strictly complete, single-copy families and concatenate them."""

    output_dir = os.path.join(options.output_dir, "Gene_Groups_Output")
    os.makedirs(output_dir, exist_ok=True)
    genomes = list(ledger.genomes)
    input_is_dna = detect_sequence_type(options.fasta)
    concatenated = {genome: "" for genome in genomes}
    alignment_manifest = []
    families = ledger.sorted_families(family_ids)
    if not families:
        raise FamilyInputError("No families matched the requested alignment groups")
    invalid = []
    family_members_by_genome = {}
    for family in families:
        by_genome = OrderedDict()
        for member in family.members:
            genome = ledger.resolver.genome(member)
            by_genome.setdefault(genome, []).append(member)
        missing = sorted(set(genomes) - set(by_genome))
        multicopy = sorted(genome for genome, members in by_genome.items() if len(members) != 1)
        if missing or multicopy:
            invalid.append((family.family_id, missing, multicopy))
        else:
            family_members_by_genome[family.family_id] = {
                genome: members[0] for genome, members in by_genome.items()
            }
    if invalid:
        details = []
        for family_id, missing, multicopy in invalid[:10]:
            parts = []
            if missing:
                parts.append(f"missing={','.join(missing[:5])}")
            if multicopy:
                parts.append(f"multicopy={','.join(multicopy[:5])}")
            details.append(f"{family_id} ({'; '.join(parts)})")
        raise FamilyInputError(
            f"Core alignment requires exactly one sequence from every genome; "
            f"{len(invalid)} selected family/families failed: {'; '.join(details)}"
        )

    for family in families:
        by_genome = family_members_by_genome[family.family_id]
        input_path = os.path.join(output_dir, f".{family.family_id}.mafft_input.fasta")
        aligned_path = os.path.join(output_dir, f".{family.family_id}.mafft_aligned.fasta")
        with open(input_path, "w", encoding="utf-8") as handle:
            for genome, member in sorted(by_genome.items()):
                sequence = ledger.sequences[member]
                if getattr(options, "align_aa", False) and input_is_dna:
                    sequence = translate_frame(sequence)
                handle.write(f">{genome}\n{wrap_sequence(sequence, 60)}\n")
        try:
            if len(by_genome) > 1:
                _run_mafft(options, input_path, aligned_path)
                aligned = read_fasta_records(aligned_path)
            else:
                aligned = read_fasta_records(input_path)
            width = max(len(sequence) for sequence in aligned.values())
            for genome in genomes:
                concatenated[genome] += aligned.get(genome, "-" * width)
            alignment_manifest.append((family.family_id, width, family.category))
        finally:
            for temporary in (input_path, aligned_path):
                try:
                    os.remove(temporary)
                except FileNotFoundError:
                    pass
    output_path = os.path.join(options.output_dir, "core_gene_alignment.aln")
    with open(output_path, "w", encoding="utf-8") as handle:
        for genome, sequence in concatenated.items():
            handle.write(f">{genome}\n{wrap_sequence(sequence, 60)}\n")
    manifest_path = os.path.join(options.output_dir, "core_gene_alignment_manifest.tsv")
    with open(manifest_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["family_id", "alignment_columns", "category"])
        writer.writerows(alignment_manifest)


def _write_statistics(
    path, ledger, selected_ids, categories, mode, total_sequences,
    groups=None, single_copy_filter=None,
):
    """Write clearly separated reconciliation-stage and prevalence summaries.

    Stage describes how a final family arose during supplementation. Prevalence
    describes how many genomes/taxa contain that final family. The two are
    deliberately reported as separate dimensions.
    """
    selected = set(selected_ids)
    selected_families = [ledger.families[family_id] for family_id in selected]
    stage_order = ("First", "extended", "combined", "only_Second")
    stage_counts = Counter(family.stage for family in selected_families)
    taxa_total = len(ledger.taxa)

    if groups is None:
        raw_groups = {
            family.category.rsplit("_", 1)[-1]
            for family in ledger.families.values()
            if getattr(family, "category", "")
        }
        if mode == "Species":
            groups = sorted((int(item) for item in raw_groups), reverse=True)
        else:
            groups = sorted((int(item) for item in raw_groups if item != ">"))
            if ">" in raw_groups:
                groups.append(">")

    group_keys = [str(group) for group in groups]
    stage_group_counts = {stage: Counter() for stage in stage_order}
    for family in selected_families:
        group = str(_prevalence_bin(family.prevalence, taxa_total, mode, groups))
        stage_group_counts.setdefault(family.stage, Counter())[group] += 1

    overall_group_counts = Counter()
    for counts in stage_group_counts.values():
        overall_group_counts.update(counts)

    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f"Analysis mode: {mode}\n")
        handle.write(f"Number of genomes: {len(ledger.genomes)}\n")
        handle.write(f"Number of taxa: {taxa_total}\n")
        handle.write(f"Input sequences: {total_sequences}\n")
        handle.write(f"Sequences assigned exactly once: {len(ledger.assignment)}\n")
        handle.write(f"Final gene families: {len(ledger.families)}\n")
        handle.write(f"Families retained after filters: {len(selected)}\n")

        filter_details = single_copy_filter or {
            "enabled": False,
            "scope": "final_reconciled_families",
            "tolerance_percent": None,
            "retained_families": len(selected),
            "excluded_families": len(ledger.families) - len(selected),
        }
        handle.write("\nSINGLE-COPY FILTER\n")
        handle.write(
            "Scope: final reconciled families, counting every member sequence "
            "from original and supplemental inputs\n"
        )
        handle.write(
            f"Applied: {'yes' if filter_details.get('enabled') else 'no'}\n"
        )
        if filter_details.get("enabled"):
            handle.write(
                f"Multi-copy tolerance: {filter_details.get('tolerance_percent'):.4f}% "
                "of occupied genomes\n"
            )
            handle.write(
                "Criterion: number of genomes with more than one family member "
                "must not exceed floor(occupied genomes x tolerance / 100)\n"
            )
        handle.write(
            f"Families retained: {filter_details.get('retained_families', len(selected)):,}\n"
        )
        handle.write(
            f"Families excluded: {filter_details.get('excluded_families', len(ledger.families) - len(selected)):,}\n"
        )

        handle.write(
            "\nRECONCILIATION PROVENANCE (technical classification)\n"
        )
        handle.write(
            "The following labels describe how final families were formed "
            "during clustering and reconciliation. They do not indicate "
            "whether a family is biologically core or accessory.\n"
        )

        for stage in stage_order:
            display_label = _reconciliation_stage_label(stage)

            handle.write(
                f"{display_label}: {stage_counts.get(stage, 0)}\n"
            )

        if not any(
            stage_counts.get(stage, 0)
            for stage in ("extended", "combined", "only_Second")
        ):
            handle.write(
                "This is a standalone baseline pangenome analysis; no "
                "supplemental genomes have been reconciled into this output.\n"
            )

        handle.write(
            "\nRECONCILIATION PROVENANCE BY FINAL FAMILY PREVALENCE\n"
        )

        if mode == "Species":
            handle.write(
                f"Prevalence denominator: {taxa_total} genomes. "
                "Each family appears in exactly one mutually exclusive bin.\n"
            )
        else:
            handle.write(
                f"Prevalence denominator: {taxa_total} taxa. "
                "Each family appears in exactly one configured genus bin.\n"
            )

        labels = [
            _prevalence_bin_label(group, groups, mode)
            for group in group_keys
        ]

        handle.write(
            "Reconciliation process\t"
            + "\t".join(labels)
            + "\tTotal\n"
        )

        for stage in stage_order:
            display_label = _reconciliation_stage_label(stage)

            values = [
                stage_group_counts[stage].get(group, 0)
                for group in group_keys
            ]

            handle.write(
                display_label
                + "\t"
                + "\t".join(str(value) for value in values)
                + "\t"
                + str(sum(values))
                + "\n"
            )

        overall_values = [
            overall_group_counts.get(group, 0)
            for group in group_keys
        ]

        handle.write(
            "ALL_FAMILIES\t"
            + "\t".join(str(value) for value in overall_values)
            + "\t"
            + str(sum(overall_values))
            + "\n"
        )

        if mode == "Species":
            handle.write("\nCumulative prevalence across all family stages\n")
            handle.write(
                "These are cumulative counts, unlike the mutually exclusive bins above.\n"
            )
            for threshold in groups:
                threshold = int(threshold)
                if threshold <= 0:
                    continue
                minimum_genomes = math.ceil(threshold * taxa_total / 100.0)
                count = sum(
                    family.prevalence >= minimum_genomes
                    for family in selected_families
                )
                handle.write(
                    f"Families present in >= {threshold}% of genomes "
                    f"(at least {minimum_genomes}/{taxa_total}): {count}\n"
                )


def _read_tsv_dicts(path):
    """Read a tab-separated PyamilySeq family-summary file."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Required family summary was not found: {path}")
    with open(path, "r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _read_json_dict(path):
    if not os.path.isfile(path):
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        value = json.load(handle)
    return value if isinstance(value, dict) else {}


def _source_family_tokens(row):
    """Return stable source-family IDs, excluding supplemental singleton labels."""
    raw = row.get("source_families", "") or ""
    tokens = {item.strip() for item in raw.split(";") if item.strip().startswith("source_")}
    # Older family-summary files still carry stable group_<hash> family IDs.
    if not tokens:
        family_id = str(row.get("family_id", "") or "")
        if family_id.startswith("group_"):
            tokens.add("source_" + family_id[len("group_"):])
        elif family_id.startswith("source_"):
            tokens.add(family_id)
    return tokens


def _row_prevalence_count(row):
    for key in ("prevalence_count", "taxa_count", "genome_count"):
        value = row.get(key)
        if value not in (None, ""):
            try:
                return int(float(value))
            except (TypeError, ValueError):
                pass
    return None


def _row_prevalence_denominator(row, fallback):
    value = row.get("prevalence_denominator")
    if value not in (None, ""):
        try:
            return int(float(value))
        except (TypeError, ValueError):
            pass
    return int(fallback) if fallback else None


def _manifest_species_count(manifest):
    # In Species mode, each taxon is one genome. Prefer the explicit taxa field.
    for key in ("taxa", "genomes"):
        value = manifest.get(key)
        if value not in (None, ""):
            try:
                return int(value)
            except (TypeError, ValueError):
                pass
    return None


def _configured_core_threshold(manifest, fallback_groups=None):
    parameters = manifest.get("parameters", {}) or {}
    group_text = parameters.get("species_groups")
    if group_text:
        try:
            values = [int(item.strip()) for item in str(group_text).split(",") if item.strip()]
        except ValueError:
            values = []
    else:
        values = [int(item) for item in (fallback_groups or []) if str(item).isdigit()]
    positive = [value for value in values if value > 0]
    return max(positive) if positive else 99


def _core_threshold_met(count, denominator, threshold):
    if count is None or denominator is None or denominator <= 0:
        return False
    return count >= math.ceil(float(threshold) * denominator / 100.0)


def _prevalence_band_for_count(count, denominator, groups):
    if count is None or denominator is None or denominator <= 0:
        return "unknown"
    group_values = sorted({int(item) for item in groups if str(item).isdigit() and int(item) > 0}, reverse=True)
    if not group_values:
        group_values = [99, 95, 15]
    selected = "0"
    for group in group_values:
        if _core_threshold_met(count, denominator, group):
            selected = str(group)
            break
    return _prevalence_bin_label(selected, group_values + [0], "Species")


def _append_supplementation_impact(
    output_dir, baseline_output_dir, final_manifest, mode, groups
):
    """Append a baseline-core lineage comparison to the supplemented report.

    Original families are tracked into final families with their stable
    ``source_families`` identifiers. The report distinguishes baseline
    lineages from distinct final families because several source families can
    merge into one final family during reclustering.
    """
    if mode != "Species":
        raise ValueError("Baseline core lineage reporting is currently implemented for Species mode only")

    baseline_rows = _read_tsv_dicts(os.path.join(baseline_output_dir, "family_summary.tsv"))
    final_rows = _read_tsv_dicts(os.path.join(output_dir, "family_summary.tsv"))
    baseline_manifest = _read_json_dict(os.path.join(baseline_output_dir, "run_manifest.json"))

    baseline_genomes = _manifest_species_count(baseline_manifest)
    final_genomes = _manifest_species_count(final_manifest)
    if baseline_genomes is None:
        baseline_genomes = next((
            _row_prevalence_denominator(row, None) for row in baseline_rows
            if _row_prevalence_denominator(row, None) is not None
        ), None)
    if final_genomes is None:
        final_genomes = next((
            _row_prevalence_denominator(row, None) for row in final_rows
            if _row_prevalence_denominator(row, None) is not None
        ), None)
    if not baseline_genomes or not final_genomes:
        raise ValueError("Could not determine original and supplemented genome denominators from manifests/family summaries")

    core_threshold = _configured_core_threshold(baseline_manifest, groups)
    comparison_groups = sorted(
        {int(item) for item in groups if str(item).isdigit() and int(item) > 0} | {core_threshold},
        reverse=True,
    ) + [0]

    baseline_source_to_family = {}
    all_baseline_source_tokens = set()
    baseline_core_rows = []
    for row in baseline_rows:
        family_id = str(row.get("family_id", "") or "")
        count = _row_prevalence_count(row)
        denominator = _row_prevalence_denominator(row, baseline_genomes)
        tokens = _source_family_tokens(row)
        for token in tokens:
            all_baseline_source_tokens.add(token)
            baseline_source_to_family.setdefault(token, set()).add(family_id)
        if _core_threshold_met(count, denominator, core_threshold):
            baseline_core_rows.append({
                "family_id": family_id,
                "source_tokens": tokens,
                "count": count,
                "denominator": denominator,
            })

    baseline_core_tokens = set()
    for item in baseline_core_rows:
        baseline_core_tokens.update(item["source_tokens"])
    baseline_noncore_tokens = all_baseline_source_tokens - baseline_core_tokens

    final_by_id = {}
    final_rows_by_source = {}
    final_core_rows = []
    for row in final_rows:
        family_id = str(row.get("family_id", "") or "")
        if not family_id:
            continue
        final_by_id[family_id] = row
        tokens = _source_family_tokens(row)
        for token in tokens & all_baseline_source_tokens:
            final_rows_by_source.setdefault(token, set()).add(family_id)
        count = _row_prevalence_count(row)
        denominator = _row_prevalence_denominator(row, final_genomes)
        if _core_threshold_met(count, denominator, core_threshold):
            final_core_rows.append(row)

    added_genome_names = str(
        ((final_manifest.get("parameters", {}) or {}).get("added_genomes")) or ""
    )
    added_genomes = [item for item in added_genome_names.split(",") if item]

    # Read actual genome columns when GPA output exists. Do not infer which
    # supplemental genomes are present from aggregate prevalence arithmetic.
    gpa_presence = {}
    gpa_columns = set()
    gpa_path = os.path.join(output_dir, "gene_presence_absence.csv")
    if added_genomes and os.path.isfile(gpa_path):
        with open(gpa_path, "r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            gpa_columns = set(reader.fieldnames or [])
            if all(genome in gpa_columns for genome in added_genomes):
                for row in reader:
                    family_id = str(row.get("Gene", "") or "")
                    if family_id:
                        gpa_presence[family_id] = {
                            genome: bool(str(row.get(genome, "") or "").strip())
                            for genome in added_genomes
                        }

    fate_counts = Counter()
    lineage_outcome_counts = Counter()

    merged_lineages = 0
    mapped_one = 0
    unique_final_descendants = set()

    # These are counts of ORIGINAL CORE LINEAGES, not genome counts.
    original_core_presence_counts = Counter()

    detailed_rows = []

    for item in baseline_core_rows:
        source_tokens = item["source_tokens"]
        final_ids = set()
        for token in source_tokens:
            final_ids.update(final_rows_by_source.get(token, set()))
        final_ids = sorted(final_ids)
        unique_final_descendants.update(final_ids)

        final_stage_text = ""
        final_stage_label_text = ""
        supplementation_outcome = ""
        final_prevalence_text = ""
        final_bin_text = ""

        final_count = None
        final_denominator = None



        merge_count = 0

        added_present = ""
        added_absent = ""
        added_presence_text = ""

        if not source_tokens or not final_ids:
            fate = "unmapped"
            mapping_status = "unmapped"
            supplementation_outcome = "Unmapped"
        elif len(final_ids) > 1:
            fate = "split_across_multiple_final_families"
            mapping_status = "split"
            final_stage_text = ";".join(str(final_by_id[fid].get("stage", "")) for fid in final_ids)
            supplementation_outcome = "Split across final families"

            final_stage_label_text = ";".join(
                _reconciliation_stage_label(
                    final_by_id[fid].get("stage", "")
                )
                for fid in final_ids
            )
            final_prevalence_text = ";".join(
                f"{_row_prevalence_count(final_by_id[fid])}/{_row_prevalence_denominator(final_by_id[fid], final_genomes)}"
                for fid in final_ids
            )
            final_bin_text = ";".join(
                _prevalence_band_for_count(
                    _row_prevalence_count(final_by_id[fid]),
                    _row_prevalence_denominator(final_by_id[fid], final_genomes),
                    comparison_groups,
                ) for fid in final_ids
            )
        else:
            mapping_status = "single_final_family"
            mapped_one += 1
            final_id = final_ids[0]
            final_row = final_by_id[final_id]
            final_stage_text = str(final_row.get("stage", "") or "")
            final_stage_label_text = _reconciliation_stage_label(
                final_stage_text
            )
            final_count = _row_prevalence_count(final_row)
            final_denominator = _row_prevalence_denominator(final_row, final_genomes)
            if (
                final_count is None
                or final_denominator is None
                or final_denominator <= 0
            ):
                supplementation_outcome = "Prevalence unavailable"

            elif _core_threshold_met(
                    final_count,
                    final_denominator,
                    core_threshold,
            ):
                supplementation_outcome = "Retained core"

            else:
                supplementation_outcome = "Core downgraded"
            final_prevalence_text = (
                f"{final_count}/{final_denominator}" if final_count is not None and final_denominator else ""
            )
            final_bin_text = _prevalence_band_for_count(
                final_count, final_denominator, comparison_groups
            )
            fate = final_bin_text if final_bin_text != "unknown" else "prevalence_unavailable"

            final_tokens = _source_family_tokens(final_row) & all_baseline_source_tokens
            merge_count = len({
                original_family_id
                for token in final_tokens
                for original_family_id in baseline_source_to_family.get(token, set())
            })
            if merge_count > 1:
                merged_lineages += 1

            if added_genomes:
                statuses = gpa_presence.get(final_id)

                if statuses is None:
                    original_core_presence_counts["unavailable"] += 1

                else:
                    present = [
                        genome for genome, is_present in statuses.items()
                        if is_present
                    ]
                    absent = [
                        genome for genome, is_present in statuses.items()
                        if not is_present
                    ]

                    # These two TSV fields represent counts of added genomes
                    # within the final family corresponding to this lineage.
                    added_present = str(len(present))
                    added_absent = str(len(absent))

                    added_presence_text = ";".join(
                        f"{genome}={'present' if statuses[genome] else 'absent'}"
                        for genome in added_genomes
                    )

                    # The Counter represents original core lineages.
                    if not absent:
                        original_core_presence_counts[
                            "all_added_genomes_present"
                        ] += 1

                    elif not present:
                        original_core_presence_counts[
                            "no_added_genome_present"
                        ] += 1

                    else:
                        original_core_presence_counts[
                            "some_added_genomes_present"
                        ] += 1

        # Presence cannot currently be assigned to one final family when an
        # original lineage is unmapped or split across multiple final families.
        if added_genomes and mapping_status in {"unmapped", "split"}:
            original_core_presence_counts["not_assessed"] += 1

        fate_counts[fate] += 1
        lineage_outcome_counts[supplementation_outcome] += 1
        detailed_rows.append({
            "baseline_family_id": item["family_id"],
            "baseline_status": "Original core",
            "baseline_source_ids": ";".join(sorted(source_tokens)),
            "baseline_prevalence": f"{item['count']}/{item['denominator']}",

            "mapping_status": mapping_status,
            "final_family_ids": ";".join(final_ids),

            "final_stage": final_stage_text,
            "final_stage_label": final_stage_label_text,

            "final_prevalence": final_prevalence_text,
            "final_prevalence_bin": final_bin_text,

            "fate_of_baseline_core_family": fate,
            "supplementation_outcome": supplementation_outcome,

            "number_of_original_families_in_final_destination": merge_count,

            # Explicitly named columns: counts of added genomes present/absent
            # in this lineage's mapped final family (not lineage counts).
            "added_genome_count_present": added_present,
            "added_genome_count_absent": added_absent,
            # Retain legacy column names for compatibility with existing readers.
            "added_genomes_present": added_present,
            "added_genomes_absent": added_absent,
            "added_genome_presence_by_id": added_presence_text,
        })

    # Classify each final core family by the original lineage(s) it contains.
    final_core_origin_counts = Counter()
    final_core_with_baseline_core_ids = set()
    for row in final_core_rows:
        tokens = _source_family_tokens(row)
        has_core = bool(tokens & baseline_core_tokens)
        has_original_noncore = bool(tokens & baseline_noncore_tokens)
        if has_core:
            final_core_origin_counts["contains_original_core_lineage"] += 1
            final_core_with_baseline_core_ids.add(str(row.get("family_id", "")))
        elif has_original_noncore:
            final_core_origin_counts["contains_original_noncore_lineage_only"] += 1
        else:
            final_core_origin_counts["no_original_source_lineage_identified"] += 1

    fate_path = os.path.join(output_dir, "supplementation_core_lineage.tsv")
    fate_fields = [
        "baseline_family_id",
        "baseline_status",
        "baseline_source_ids",
        "baseline_prevalence",
        "mapping_status",
        "final_family_ids",
        "final_stage",
        "final_stage_label",
        "final_prevalence",
        "final_prevalence_bin",
        "fate_of_baseline_core_family",
        "supplementation_outcome",
        "number_of_original_families_in_final_destination",
        "added_genome_count_present",
        "added_genome_count_absent",
        # Legacy aliases retained for existing downstream readers.
        "added_genomes_present",
        "added_genomes_absent",
        "added_genome_presence_by_id",
    ]
    with open(fate_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fate_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(detailed_rows)

    baseline_core_total = len(baseline_core_rows)
    final_core_total = len(final_core_rows)
    delta = final_core_total - baseline_core_total
    delta_pct = (100.0 * delta / baseline_core_total) if baseline_core_total else 0.0
    summary_path = os.path.join(output_dir, "summary_statistics.txt")
    existing_summary = ""
    if os.path.isfile(summary_path):
        with open(summary_path, "r", encoding="utf-8") as handle:
            existing_summary = handle.read()
    marker = "SUPPLEMENTATION IMPACT: TRACKING THE ORIGINAL PANGENOME"
    marker_pos = existing_summary.find(marker)
    if marker_pos >= 0:
        section_start = existing_summary.rfind("\n\n", 0, marker_pos)
        if section_start >= 0:
            existing_summary = existing_summary[:section_start]
        else:
            existing_summary = existing_summary[:marker_pos]
    with open(summary_path, "w", encoding="utf-8") as handle:
        if existing_summary.strip():
            handle.write(existing_summary.rstrip() + "\n\n")
        handle.write("SUPPLEMENTATION IMPACT: TRACKING THE ORIGINAL PANGENOME\n")
        handle.write("======================================================\n")
        handle.write("BASELINE PANGENOME (before any supplemental genomes were added)\n")
        handle.write(f"Original genomes: {baseline_genomes:,}\n")
        handle.write(f"Original gene families: {len(baseline_rows):,}\n")
        handle.write(
            f"Original core definition: >= {core_threshold}% prevalence "
            f"(at least {math.ceil(core_threshold * baseline_genomes / 100.0)}/{baseline_genomes} genomes)\n"
        )
        handle.write(f"Original core families: {baseline_core_total:,}\n")

        handle.write("\nSUPPLEMENTED PANGENOME\n")
        handle.write(f"Total genomes after supplementation: {final_genomes:,}\n")
        handle.write(f"Additional genomes: {len(added_genomes):,}\n")
        if added_genomes:
            handle.write("Additional genome IDs: " + ", ".join(added_genomes) + "\n")
        handle.write(f"Final gene families: {len(final_rows):,}\n")
        handle.write(
            f"Final core families at the same >= {core_threshold}% threshold: {final_core_total:,}\n"
        )
        handle.write(
            f"Net change in the total number of core families: {delta:+,} "
            f"({delta_pct:+.2f}% relative to the original core count)\n"
        )

        handle.write("\nFATE OF THE ORIGINAL CORE FAMILY SET\n")
        handle.write(f"Original core family lineages tracked: {baseline_core_total:,}\n")
        handle.write(f"Mapped to exactly one final family: {mapped_one:,}\n")
        handle.write(f"Split across multiple final families: {fate_counts.get('split_across_multiple_final_families', 0):,}\n")
        handle.write(f"Not mapped to a final family: {fate_counts.get('unmapped', 0):,}\n")
        handle.write(
            f"Distinct final families containing original core lineages (all prevalence levels): "
            f"{len(unique_final_descendants):,}\n"
        )
        for label in [
            _prevalence_bin_label(str(group), comparison_groups, "Species")
            for group in comparison_groups if int(group) > 0
        ]:
            handle.write(
                f"Original core lineages whose final family is in {label}: "
                f"{fate_counts.get(label, 0):,}\n"
            )
        if merged_lineages:
            handle.write(
                f"Original core lineages whose final family also contains other original families: "
                f"{merged_lineages:,}\n"
            )
        handle.write(
            "\nBIOLOGICAL OUTCOMES OF THE ORIGINAL CORE LINEAGES\n"
        )

        for outcome in (
                "Retained core",
                "Core downgraded",
                "Split across final families",
                "Unmapped",
                "Prevalence unavailable",
        ):
            handle.write(
                f"{outcome}: {lineage_outcome_counts.get(outcome, 0):,}\n"
            )

        handle.write(
            f"Total original core lineages classified: "
            f"{sum(lineage_outcome_counts.values()):,}\n"
        )
        if added_genomes:
            handle.write(
                "\nADDED-GENOME PRESENCE IN FINAL FAMILIES "
                "DESCENDED FROM ORIGINAL CORE\n"
            )
            handle.write(
                "Counts below represent original core family lineages, not "
                "genome counts or unique final-family counts. A lineage is "
                "counted once according to the presence of the added genome(s) "
                "in its mapped final family.\n"
            )
            handle.write(
                f"Number of added genomes in this run: {len(added_genomes):,}\n"
            )

            if not all(genome in gpa_columns for genome in added_genomes):
                handle.write(
                    "WARNING: One or more added-genome columns were unavailable "
                    "in gene_presence_absence.csv; affected lineages are reported "
                    "as presence unavailable.\n"
                )

            handle.write(
                "Original core lineages whose final family contains all added genomes: "
                f"{original_core_presence_counts.get('all_added_genomes_present', 0):,}\n"
            )
            handle.write(
                "Original core lineages whose final family contains some but not all "
                "added genomes: "
                f"{original_core_presence_counts.get('some_added_genomes_present', 0):,}\n"
            )
            handle.write(
                "Original core lineages whose final family contains none of the "
                "added genomes: "
                f"{original_core_presence_counts.get('no_added_genome_present', 0):,}\n"
            )
            handle.write(
                "Presence unavailable: "
                f"{original_core_presence_counts.get('unavailable', 0):,}\n"
            )
            handle.write(
                "Presence not assessed (split or unmapped lineages): "
                f"{original_core_presence_counts.get('not_assessed', 0):,}\n"
            )

            presence_total = sum(original_core_presence_counts.values())
            handle.write(
                f"Total original core lineages accounted for: "
                f"{presence_total:,}/{baseline_core_total:,}\n"
            )

        handle.write(
            "\nORIGIN OF THE FINAL SUPPLEMENTED CORE SET\n"
        )

        handle.write(
            f"Final core families containing at least one original-core lineage: "
            f"{final_core_origin_counts.get('contains_original_core_lineage', 0):,}\n"
        )
        handle.write(
            f"Newly core families derived exclusively from original non-core lineages: "
            f"{final_core_origin_counts.get('contains_original_noncore_lineage_only', 0):,}\n"
        )
        handle.write(
            f"Final core families with unresolved original source provenance: "
            f"{final_core_origin_counts.get('no_original_source_lineage_identified', 0):,}\n"
        )
        handle.write("\nINTERPRETATION\n")
        handle.write(
            "The lineage counts refer to original core families, tracked by stable source-family IDs. "
            "Several original families can merge into one final family; those cases are counted once "
            "per original lineage in the lineage sections, while final-core origin counts represent "
            "distinct final families. A split or unmapped lineage is not assigned an inferred final "
            "prevalence outcome.\n"
        )
        handle.write(
            "Added-genome presence categories count original core lineages, not genomes or unique "
            "final families. Presence is read from gene_presence_absence.csv when the required "
            "added-genome columns are available. Reconciliation provenance labels describe how "
            "families were formed; they are not biological core/accessory classifications.\n"
        )

    return {
        "baseline_core_families": baseline_core_total,
        "final_core_families": final_core_total,
        "baseline_core_lineages_mapped_to_one_final_family": mapped_one,
        "distinct_final_descendants": len(unique_final_descendants),
        "lineage_fate_counts": dict(fate_counts),
        "lineage_outcome_counts": dict(lineage_outcome_counts),
        "added_genome_presence_lineage_counts": dict(original_core_presence_counts),
        "lineage_tsv": fate_path,
    }


def _write_run_manifest(
    path, options, ledger, selected_ids, mode, total_sequences,
    single_copy_filter=None,
):
    manifest = {
        'schema_version': 1,
        'analysis_mode': mode,
        'cluster_format': options.cluster_format,
        'cluster_file': os.path.abspath(options.clusters),
        'fasta_file': os.path.abspath(options.fasta) if getattr(options, 'fasta', None) else None,
        'sequence_db': os.path.abspath(options.sequence_db) if getattr(options, 'sequence_db', None) else None,
        'sequence_db_out': os.path.abspath(options.sequence_db_out) if getattr(options, 'sequence_db_out', None) else None,
        'sequence_db_path': os.path.abspath(ledger.sequences.path),
        'sequence_db_retained': bool(getattr(ledger.sequences, 'persistent', False)),
        'gene_annotation_files': (
            list(getattr(options, 'gene_annotation_files', []) or [])
            if not isinstance(getattr(options, 'gene_annotation_files', None), str)
            else [getattr(options, 'gene_annotation_files')]
        ),
        'reclustered_file': (
            os.path.abspath(options.reclustered) if getattr(options, 'reclustered', None) else None
        ),
        'input_sequences': total_sequences,
        'assigned_sequences': len(ledger.assignment),
        'exactly_once_assignment': len(ledger.assignment) == total_sequences,
        'final_families': len(ledger.families),
        'selected_families': len(selected_ids),
        'single_copy_filter': single_copy_filter or {
            'enabled': False,
            'scope': 'final_reconciled_families',
            'tolerance_percent': None,
            'retained_families': len(selected_ids),
            'excluded_families': len(ledger.families) - len(selected_ids),
        },
        'genomes': len(ledger.genomes),
        'taxa': len(ledger.taxa),
        'stage_counts': dict(Counter(family.stage for family in ledger.families.values())),
        'parameters': {
            key: value for key, value in vars(options).items()
            if isinstance(value, (str, int, float, bool, type(None)))
        },
    }
    with open(path, 'w', encoding='utf-8') as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)


def run_pangenome_workflow(options, mode):
    logger = logging.getLogger("PyamilySeq.PyamilySeq")
    logger.info("Starting %s workflow", mode)
    ledger = build_family_ledger(
        fasta_path=options.fasta,
        cluster_path=options.clusters,
        cluster_format=options.cluster_format,
        mode=mode,
        recluster_path=getattr(options, "reclustered", None),
        sequence_tag=getattr(options, "sequence_tag", None) or "StORF",
        taxon_map=getattr(options, "taxon_map", None),
        sequence_delimiter=getattr(options, "sequence_delimiter", "|") or "|",
        options=options,
        logger=logger,
        sequence_db_path=getattr(options, "sequence_db", None),
        supplemental_sequence_ids=getattr(options, "supplemental_sequence_ids", None),
        sequence_db_out_path=getattr(options, "sequence_db_out", None),
    )
    try:
        total_sequences = len(ledger.sequences)
        logger.info("Family ledger ready: %s sequences, %s final families", total_sequences, len(ledger.families))
        group_text = options.species_groups if mode == "Species" else options.core_groups
        categories, groups = classify_families(ledger, mode, group_text)
        logger.info("Classified families into %s categories", len(categories))

        all_family_ids = list(ledger.families)
        single_copy_enabled = bool(getattr(options, "single_copy_only", False))
        raw_tolerance = getattr(options, "single_copy_tolerance", 0.0)
        try:
            tolerance = float(raw_tolerance or 0.0)
        except (TypeError, ValueError):
            raise FamilyInputError(
                "Single-copy tolerance must be a finite number between 0 and 100"
            ) from None
        if not math.isfinite(tolerance) or tolerance < 0.0 or tolerance > 100.0:
            raise FamilyInputError(
                "Single-copy tolerance must be a finite number between 0 and 100"
            )
        single_copy_enabled = single_copy_enabled or tolerance > 0.0

        if single_copy_enabled:
            # This operates on the final ledger, after First and Second
            # families have been reconciled.  Do not move this filter into the
            # individual clustering parsers: doing so would miss cross-stage
            # and cross-source-family copy-number conflicts.
            selected_ids = ledger.selected_single_copy(tolerance)
        else:
            selected_ids = all_family_ids

        single_copy_filter = {
            "enabled": single_copy_enabled,
            "scope": "final_reconciled_families",
            "tolerance_percent": tolerance if single_copy_enabled else None,
            "counting_unit": "member_sequences_per_genome",
            "includes_original_and_supplemental_members": True,
            "retained_families": len(selected_ids),
            "excluded_families": len(all_family_ids) - len(selected_ids),
        }
        if single_copy_enabled:
            logger.info(
                "Single-copy filter evaluated %s final reconciled families at %.4f%% tolerance: "
                "%s retained, %s excluded",
                len(all_family_ids),
                tolerance,
                len(selected_ids),
                len(all_family_ids) - len(selected_ids),
            )

        output_dir = os.path.abspath(options.output_dir)
        os.makedirs(output_dir, exist_ok=True)
        logger.info("Writing assignment and summary outputs to %s", output_dir)
        ledger.write_assignments(os.path.join(output_dir, "gene_family_assignments.tsv"))
        _write_family_summary(
            os.path.join(output_dir, "family_summary.tsv"), ledger, selected_ids, mode, groups
        )
        _write_statistics(
            os.path.join(output_dir, "summary_statistics.txt"),
            ledger,
            selected_ids,
            categories,
            mode,
            total_sequences,
            groups,
            single_copy_filter,
        )
        _write_run_manifest(
            os.path.join(output_dir, 'run_manifest.json'),
            options,
            ledger,
            selected_ids,
            mode,
            total_sequences,
            single_copy_filter,
        )
        if getattr(options, "gene_presence_absence_out", True):
            logger.info("Writing gene_presence_absence.csv")
            _write_gpa(
                os.path.join(output_dir, "gene_presence_absence.csv"),
                ledger,
                selected_ids,
                mode,
                gene_annotation_files=getattr(options, "gene_annotation_files", None),
                temp_dir=getattr(options, "tmp_dir", None),
            )
        baseline_output_dir = getattr(options, "baseline_output_dir", None)
        if baseline_output_dir:
            if mode == "Species":
                logger.info("Writing original-core supplementation impact report")
                _append_supplementation_impact(
                    output_dir, baseline_output_dir,
                    _read_json_dict(os.path.join(output_dir, "run_manifest.json")),
                    mode, groups,
                )
            else:
                logger.warning(
                    "Original-core lineage report was not generated: it currently supports Species mode only."
                )
        if getattr(options, "gene_group_table", False):
            logger.info("Writing gene_group_table.tsv")
            _write_gene_group_table(os.path.join(output_dir, "gene_group_table.tsv"), ledger, selected_ids)
        if getattr(options, "core_group_table", False):
            logger.info("Writing core_group_table.tsv")
            _write_core_group_table(os.path.join(output_dir, "core_group_table.tsv"), ledger, selected_ids)
        logger.info("Writing pan_genome_reference.fa")
        _write_reference(os.path.join(output_dir, "pan_genome_reference.fa"), ledger, selected_ids)

        requested_ids = _requested_family_ids(options, ledger, selected_ids)
        if getattr(options, "write_groups", None) or getattr(options, "write_individual_groups", False):
            _write_family_fastas(options, ledger, requested_ids)
        if getattr(options, "align_core", False):
            logger.info("Writing core_gene_alignment.aln")
            _write_alignment(options, ledger, requested_ids)
        if getattr(options, "visualise", False):
            logger.info("Writing visual outputs")
            write_ledger_visualisations(options, ledger, selected_ids, categories, mode)
        if getattr(options, "added_genomes", None):
            if not getattr(options, "reclustered", None):
                raise FamilyInputError("--added-genomes requires -reclustered with the combined reclustering output")
            try:
                from .genome_reports import write_extension_report
            except ImportError:
                from genome_reports import write_extension_report
            write_extension_report(
                output_dir,
                ledger,
                options.added_genomes,
                groups=groups,
                mode=mode,
            )

        print(
            f"Built {len(ledger.families):,} final families from {total_sequences:,} sequences; "
            f"all sequences assigned exactly once."
        )
        if len(selected_ids) != len(ledger.families):
            print(f"Retained {len(selected_ids):,} families after single-copy filtering.")
        return ledger
    finally:
        try:
            ledger.sequences.close()
        except Exception:
            logger.exception("Failed to close SQLite sequence store")