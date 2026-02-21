#!/usr/bin/env python3
"""
Produce diagnostics for occurrences where the focal family appears in neighbor positions.

Usage:
  python3 scripts/diagnose_self_neighbors.py --outdir /path/to/inspect/out --family group_3

This script looks for these files in --outdir:
 - <family>_occurrences.csv (CSV)
 - <family>_neighbors_contexts.tsv (TSV)

It writes:
 - <family>_self_diagnostics.tsv (TSV) listing occurrences where the nearest upstream neighbor (offset -1)
   maps to the same family as the focal family. Columns include raw neighbor token (from occurrences),
   mapped neighbor family (from contexts), self_positions and self_reason.

"""
import argparse
import csv
import os
from pathlib import Path


def load_occurrences(path):
    occ = {}
    with open(path, newline='') as fh:
        # occurrences.csv uses comma-delimited CSV
        r = csv.DictReader(fh)
        for row in r:
            key = row.get('gene_id')
            occ[key] = row
    return occ


def load_contexts(path):
    ctx = {}
    with open(path, newline='') as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            key = row.get('gene_id')
            ctx[key] = row
    return ctx


def main():
    parser = argparse.ArgumentParser(description='Diagnose self-neighbor occurrences for a family')
    parser.add_argument('--outdir', required=True, help='Output directory produced by inspect_family')
    parser.add_argument('--family', required=True, help='Family id, e.g. group_3')
    parser.add_argument('--write-all', action='store_true', help='Write all occurrences (mark where self==True)')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    occ_path = outdir / f"{args.family}_occurrences.csv"
    ctx_path = outdir / f"{args.family}_neighbors_contexts.tsv"

    if not occ_path.exists():
        raise SystemExit(f"Occurrences file not found: {occ_path}")
    if not ctx_path.exists():
        raise SystemExit(f"Contexts file not found: {ctx_path}")

    occ = load_occurrences(str(occ_path))
    ctx = load_contexts(str(ctx_path))

    out_path = outdir / f"{args.family}_self_diagnostics.tsv"
    with open(out_path, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['gene_id', 'contig', 'start', 'end', 'strand', 'focal_family', 'first_up_raw_token', 'first_up_mapped_family', 'self_positions', 'self_reason'])
        for gid, crow in ctx.items():
            focal_family = crow.get('family_id')
            # mapped neighbor families (first token) from contexts
            up_fams = (crow.get('neighbors_up') or '').split(';')
            first_up_mapped = up_fams[0] if up_fams and up_fams[0] else ''
            # raw tokens from occurrences
            orow = occ.get(gid) or {}
            up_raw_tokens = (orow.get('neighbors_up') or '').split(';') if orow else ['']
            first_up_raw = up_raw_tokens[0] if up_raw_tokens and up_raw_tokens[0] else ''
            self_positions = crow.get('self_positions') or ''
            self_reason = crow.get('self_reason') or ''

            is_self = (first_up_mapped == focal_family and focal_family != '')
            if is_self or args.write_all:
                writer.writerow([gid, crow.get('contig') or orow.get('contig'), crow.get('start') or orow.get('start'), crow.get('end') or orow.get('end'), crow.get('strand') or orow.get('strand'), focal_family, first_up_raw, first_up_mapped, self_positions, self_reason])

    print(f'Wrote diagnostics to: {out_path}')


if __name__ == '__main__':
    main()

