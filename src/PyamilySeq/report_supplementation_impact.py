#!/usr/bin/env python3
"""Generate/update the baseline-core impact report from an existing Supplementation run.

Example:
    python report_supplementation_impact.py \
        --original-dir Results/original \
        --supplemented-dir Results/supplemented

This is a post-processing step: it does not rerun CD-HIT or rebuild the pangenome.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

try:
    from .workflow import _append_supplementation_impact, _normalise_species_groups
except ImportError:
    from workflow import _append_supplementation_impact, _normalise_species_groups


def read_json(path: str) -> dict:
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Manifest not found: {path}")
    with open(path, "r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Manifest is not a JSON object: {path}")
    return value


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare the original core families with their descendants in a supplemented PyamilySeq pangenome."
    )
    parser.add_argument("--original-dir", required=True,
                        help="Original-stage output directory containing family_summary.tsv and run_manifest.json.")
    parser.add_argument("--supplemented-dir", required=True,
                        help="Supplemented-stage output directory containing family_summary.tsv and run_manifest.json.")
    args = parser.parse_args()

    original_dir = os.path.abspath(args.original_dir)
    supplemented_dir = os.path.abspath(args.supplemented_dir)
    final_manifest = read_json(os.path.join(supplemented_dir, "run_manifest.json"))
    mode = str(final_manifest.get("analysis_mode", "Species"))
    if mode != "Species":
        parser.error("Baseline-core lineage reporting currently supports Species mode only.")

    parameters = final_manifest.get("parameters", {}) or {}
    groups_text = parameters.get("species_groups", "99,95,15,0")
    groups = _normalise_species_groups(groups_text)
    result = _append_supplementation_impact(
        supplemented_dir, original_dir, final_manifest, mode, groups
    )
    print(f"Updated: {os.path.join(supplemented_dir, 'summary_statistics.txt')}")
    print(f"Created: {result['lineage_tsv']}")
    print(f"Original core families: {result['baseline_core_families']:,}")
    print(f"Final core families: {result['final_core_families']:,}")
    print("Original core lineage fates:")
    for label, count in sorted(result["lineage_fate_counts"].items()):
        print(f"  {label}: {count:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
