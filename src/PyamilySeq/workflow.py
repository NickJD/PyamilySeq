"""Shared Species/Genus pangenome workflow built on the canonical family ledger."""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import logging
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


def _stage_prefix(stage, mode):
    suffix = "core" if mode == "Species" else "genera"
    names = {
        "First": "First",
        "extended": "extended",
        "combined": "combined",
        "only_Second": "only_Second",
    }
    return f"{names.get(stage, stage)}_{suffix}"


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


def _write_gpa(path, ledger, family_ids, mode):
    columns = list(ledger.genomes if mode == "Species" else ledger.taxa)
    column_resolver = ledger.resolver.genome if mode == "Species" else ledger.resolver.taxon
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(GPA_METADATA_COLUMNS + columns)
        for family in ledger.sorted_families(family_ids):
            occupied = family.prevalence
            sequence_count = len(family.members)
            cells = {column: [] for column in columns}
            lengths = []
            for member in family.members:
                cell_value = member
                delimiter = ledger.resolver.sequence_delimiter
                if mode == "Species" and delimiter and delimiter in member:
                    cell_value = member.split(delimiter, 1)[1]
                cells[column_resolver(member)].append(cell_value)
                lengths.append(len(ledger.sequences[member]))
            metadata = [
                family.family_id,
                "",
                "",
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


def _write_family_summary(path, ledger, family_ids):
    fields = [
        "family_id", "category", "stage", "representative", "taxa_count",
        "genome_count", "sequence_count", "mean_copies_per_occupied_genome",
        "multicopy_genomes", "source_family_count", "source_families", "source_clusters",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for family in ledger.sorted_families(family_ids):
            genome_counts = family.copy_number
            writer.writerow({
                "family_id": family.family_id,
                "category": family.category,
                "stage": family.stage,
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


def _write_statistics(path, ledger, selected_ids, categories, mode, total_sequences):
    selected = set(selected_ids)
    stage_counts = Counter(ledger.families[family_id].stage for family_id in selected)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f"Analysis mode: {mode}\n")
        handle.write(f"Number of genomes: {len(ledger.genomes)}\n")
        handle.write(f"Number of taxa: {len(ledger.taxa)}\n")
        handle.write(f"Input sequences: {total_sequences}\n")
        handle.write(f"Sequences assigned exactly once: {len(ledger.assignment)}\n")
        handle.write(f"Final gene families: {len(ledger.families)}\n")
        handle.write(f"Families retained after filters: {len(selected)}\n")
        handle.write("\nFinal family stages\n")
        for stage in ("First", "extended", "combined", "only_Second"):
            handle.write(f"{stage}: {stage_counts.get(stage, 0)}\n")
        handle.write("\nMutually exclusive final categories\n")
        for category, family_ids in categories.items():
            handle.write(f"{category}: {sum(family_id in selected for family_id in family_ids)}\n")


def _write_run_manifest(path, options, ledger, selected_ids, mode, total_sequences):
    manifest = {
        'schema_version': 1,
        'analysis_mode': mode,
        'cluster_format': options.cluster_format,
        'cluster_file': os.path.abspath(options.clusters),
        'fasta_file': os.path.abspath(options.fasta) if getattr(options, 'fasta', None) else None,
        'sequence_db': os.path.abspath(options.sequence_db) if getattr(options, 'sequence_db', None) else None,
        'reclustered_file': (
            os.path.abspath(options.reclustered) if getattr(options, 'reclustered', None) else None
        ),
        'input_sequences': total_sequences,
        'assigned_sequences': len(ledger.assignment),
        'exactly_once_assignment': len(ledger.assignment) == total_sequences,
        'final_families': len(ledger.families),
        'selected_families': len(selected_ids),
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
    )
    total_sequences = len(ledger.sequences)
    logger.info("Family ledger ready: %s sequences, %s final families", total_sequences, len(ledger.families))
    group_text = options.species_groups if mode == "Species" else options.core_groups
    categories, _ = classify_families(ledger, mode, group_text)
    logger.info("Classified families into %s categories", len(categories))

    all_family_ids = list(ledger.families)
    if getattr(options, "single_copy_only", False) or float(
        getattr(options, "single_copy_tolerance", 0.0) or 0.0
    ) > 0:
        tolerance = float(getattr(options, "single_copy_tolerance", 0.0) or 0.0)
        selected_ids = ledger.selected_single_copy(tolerance)
    else:
        selected_ids = all_family_ids

    output_dir = os.path.abspath(options.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    logger.info("Writing assignment and summary outputs to %s", output_dir)
    ledger.write_assignments(os.path.join(output_dir, "gene_family_assignments.tsv"))
    _write_family_summary(os.path.join(output_dir, "family_summary.tsv"), ledger, selected_ids)
    _write_statistics(
        os.path.join(output_dir, "summary_statistics.txt"),
        ledger,
        selected_ids,
        categories,
        mode,
        total_sequences,
    )
    _write_run_manifest(
        os.path.join(output_dir, 'run_manifest.json'), options, ledger, selected_ids, mode, total_sequences,
    )
    if getattr(options, "gene_presence_absence_out", True):
        logger.info("Writing gene_presence_absence.csv")
        _write_gpa(os.path.join(output_dir, "gene_presence_absence.csv"), ledger, selected_ids, mode)
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
        write_extension_report(output_dir, ledger, options.added_genomes)

    print(
        f"Built {len(ledger.families):,} final families from {total_sequences:,} sequences; "
        f"all sequences assigned exactly once."
    )
    if len(selected_ids) != len(ledger.families):
        print(f"Retained {len(selected_ids):,} families after single-copy filtering.")
    try:
        ledger.sequences.close()
    except Exception:
        logger.debug("Sequence store close failed", exc_info=True)
    return ledger
