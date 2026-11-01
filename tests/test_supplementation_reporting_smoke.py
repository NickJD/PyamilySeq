#!/usr/bin/env python3
"""Synthetic smoke tests for the Supplementation reporting changes.

Run from the PyamilySeq package directory, for example:
    cd src/PyamilySeq
    python test_supplementation_reporting_smoke.py
"""
from __future__ import annotations

import csv
import json
import tempfile
from collections import OrderedDict
from pathlib import Path
from types import SimpleNamespace

from workflow import _append_supplementation_impact, _write_statistics
from genome_reports import write_extension_report


def write_tsv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def test_lineage_accounting():
    with tempfile.TemporaryDirectory() as temp_dir:
        root = Path(temp_dir)
        original = root / "original"
        supplemented = root / "supplemented"
        original.mkdir()
        supplemented.mkdir()

        baseline_rows = [
            {"family_id": "B_A", "stage": "First", "prevalence_count": "9",
             "prevalence_denominator": "9", "source_families": "source_A",
             "source_clusters": "101"},
            {"family_id": "B_B", "stage": "First", "prevalence_count": "8",
             "prevalence_denominator": "9", "source_families": "source_B",
             "source_clusters": "102"},
            {"family_id": "B_C", "stage": "First", "prevalence_count": "9",
             "prevalence_denominator": "9", "source_families": "source_C",
             "source_clusters": "103"},
            {"family_id": "B_D", "stage": "First", "prevalence_count": "9",
             "prevalence_denominator": "9", "source_families": "source_D",
             "source_clusters": "104"},
        ]
        final_rows = [
            {"family_id": "F_A", "stage": "extended", "prevalence_count": "10",
             "prevalence_denominator": "10", "source_families": "source_A",
             "source_clusters": "101"},
            {"family_id": "F_C", "stage": "extended", "prevalence_count": "9",
             "prevalence_denominator": "10", "source_families": "source_C",
             "source_clusters": "103"},
            {"family_id": "F_B", "stage": "combined", "prevalence_count": "10",
             "prevalence_denominator": "10", "source_families": "source_B",
             "source_clusters": "102"},
        ]
        fields = [
            "family_id", "stage", "prevalence_count",
            "prevalence_denominator", "source_families", "source_clusters",
        ]
        write_tsv(original / "family_summary.tsv", baseline_rows, fields)
        write_tsv(supplemented / "family_summary.tsv", final_rows, fields)

        original_manifest = {
            "analysis_mode": "Species",
            "genomes": 9,
            "taxa": 9,
            "parameters": {"species_groups": "99,95,15,0"},
        }
        final_manifest = {
            "analysis_mode": "Species",
            "genomes": 10,
            "taxa": 10,
            "parameters": {
                "species_groups": "99,95,15,0",
                "added_genomes": "NewGenome",
            },
        }
        (original / "run_manifest.json").write_text(
            json.dumps(original_manifest), encoding="utf-8"
        )
        (supplemented / "run_manifest.json").write_text(
            json.dumps(final_manifest), encoding="utf-8"
        )

        with open(supplemented / "gene_presence_absence.csv", "w", newline="",
                  encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["Gene", "NewGenome"])
            writer.writerow(["F_A", "gene_A"])
            writer.writerow(["F_C", ""])
            writer.writerow(["F_B", "gene_B"])

        (supplemented / "summary_statistics.txt").write_text(
            "Analysis mode: Species\nFinal gene families: 3\n",
            encoding="utf-8",
        )

        result = _append_supplementation_impact(
            str(supplemented), str(original), final_manifest,
            "Species", [99, 95, 15, 0],
        )
        summary = (supplemented / "summary_statistics.txt").read_text(
            encoding="utf-8"
        )
        with open(supplemented / "supplementation_core_lineage.tsv",
                  newline="", encoding="utf-8") as handle:
            lineage_rows = list(csv.DictReader(handle, delimiter="\t"))

        assert "Retained core: 1" in summary
        assert "Core downgraded: 1" in summary
        assert "Unmapped: 1" in summary
        assert "Original core lineages whose final family contains all added genomes: 1" in summary
        assert "Total original core lineages accounted for: 3/3" in summary
        assert result["lineage_outcome_counts"]["Retained core"] == 1
        assert len(lineage_rows) == 3
        assert lineage_rows[0]["added_genome_count_present"] in {"0", "1"}
        assert "added_genome_count_absent" in lineage_rows[0]


def test_statistics_has_one_overall_row():
    class DummyFamily:
        def __init__(self, family_id, stage, prevalence, category):
            self.family_id = family_id
            self.stage = stage
            self.prevalence = prevalence
            self.category = category

    class DummyLedger:
        def __init__(self):
            self.families = OrderedDict([
                ("F1", DummyFamily("F1", "First", 9, "First_core_99")),
                ("F2", DummyFamily("F2", "First", 9, "First_core_99")),
                ("F3", DummyFamily("F3", "First", 2, "First_core_15")),
            ])
            self.genomes = tuple(f"G{i}" for i in range(9))
            self.taxa = self.genomes
            self.assignment = {"s1": "F1", "s2": "F2", "s3": "F3"}

    with tempfile.TemporaryDirectory() as temp_dir:
        path = Path(temp_dir) / "summary_statistics.txt"
        ledger = DummyLedger()
        _write_statistics(
            str(path), ledger, list(ledger.families), {}, "Species", 3,
            [99, 95, 15, 0],
        )
        text = path.read_text(encoding="utf-8")
        assert "ALL_FAMILIES" in text
        assert "ALL_STAGES" not in text
        assert "standalone baseline pangenome analysis" in text


def test_added_genome_report_uses_gene_units_and_display_labels():
    genome = "NewGenome"
    added_ids = [f"{genome}|gene{i}" for i in range(1, 4)]

    class FakeResolver:
        mode = "Species"

        @staticmethod
        def genome(sequence_id):
            return sequence_id.split("|", 1)[0]

    linked_core = SimpleNamespace(
        family_id="F_LINKED",
        stage="extended",
        category="extended_core_99",
        prevalence=10,
        members=(added_ids[0], "Original1|g1"),
        source_clusters=("123",),
        source_families=("source_abc",),
    )
    supplement_only = SimpleNamespace(
        family_id="F_ONLY",
        stage="only_Second",
        category="only_Second_core_0",
        prevalence=1,
        members=(added_ids[1],),
        source_clusters=("singleton:gene2",),
        source_families=(),
    )
    linked_merged = SimpleNamespace(
        family_id="F_LINKED2",
        stage="combined",
        category="combined_core_15",
        prevalence=5,
        members=(added_ids[2], "Original2|g1"),
        source_clusters=("234", "235"),
        source_families=("source_def",),
    )
    ledger = SimpleNamespace(
        genomes=(genome,),
        resolver=FakeResolver(),
        sequences=added_ids,
        assignment={
            added_ids[0]: "F_LINKED",
            added_ids[1]: "F_ONLY",
            added_ids[2]: "F_LINKED2",
        },
        families={
            "F_LINKED": linked_core,
            "F_ONLY": supplement_only,
            "F_LINKED2": linked_merged,
        },
    )

    with tempfile.TemporaryDirectory() as temp_dir:
        write_extension_report(
            temp_dir, ledger, genome, groups=[99, 95, 15, 0], mode="Species"
        )
        root = Path(temp_dir)
        report = (root / "added_genome_report.html").read_text(encoding="utf-8")
        svg = (root / "added_genome_summary.svg").read_text(encoding="utf-8")
        summary = (root / "added_genome_summary.tsv").read_text(encoding="utf-8")
        details = (root / "added_genome_family_links.tsv").read_text(encoding="utf-8")

        assert "Novel family" not in report
        assert "Supplement-only family" in report
        assert "gene counts, not family counts" in report
        assert "Original family + supplement — &gt;=99% prevalence" in report
        assert "Supplement-only family — &lt;15% prevalence" in report
        assert "Genes in supplement-only families" in svg
        assert "genes_in_supplement_only_families" in summary
        assert "added_gene_counts_by_display_family_category" in summary
        assert "family_status_code" in details
        assert "supplement_only" in details


if __name__ == "__main__":
    test_lineage_accounting()
    test_statistics_has_one_overall_row()
    test_added_genome_report_uses_gene_units_and_display_labels()
    print("PASS: all three Supplementation reporting smoke tests")