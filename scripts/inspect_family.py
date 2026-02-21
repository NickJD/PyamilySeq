#!/usr/bin/env python3
"""
Inspect a single gene family in detail using the existing location_tracking utilities.

Usage example:
  python3 scripts/inspect_family.py --family group_3 \
    --gff-dir /path/to/gffs --clusters /path/to/gene_presence_absence.csv \
    --outdir out/inspect_group_3 --k 5

This script loads `src/PyamilySeq/location_tracking.py` at runtime (runpy) and calls
its parsing/mapping functions so we don't duplicate parsing logic.

Outputs (in --outdir):
 - occurrences.csv: detailed per-occurrence rows (same as family_locations but only for this family)
 - summary.csv: single-row family summary
 - neighbors_contexts.tsv: one line per occurrence with neighbor lists and simple context

Extra feature suggestions (printed at the end):
 - neighbor family aggregation: map neighbor gene IDs to families and show family-level synteny
 - interactive HTML visualization of synteny with links to genome browser coordinates
 - sequence extraction: produce a FASTA of each gene occurrence if genome FASTA is available
 - alignment-based comparison: run pairwise nucleotide alignment of occurrences to measure divergence
 - phylogenetic context: overlay family presence on species tree
 - optional output of neighbor family IDs (requires clusters->member->family reverse map)

The script is intentionally dependency-free (only standard library) and reuses the
project's parsing code. It expects `location_tracking.py` to be located at ../src/PyamilySeq/
relative to this script (standard project layout).
"""

import argparse
import csv
import os
import runpy
import sys
from pathlib import Path
from collections import defaultdict, Counter


def write_csv_rows(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def main():
    parser = argparse.ArgumentParser(description='Inspect a single gene family in detail')
    parser.add_argument('--family', required=True, help='Family id to inspect, e.g. group_3')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--gff-dir', help='Directory containing per-genome GFF files')
    group.add_argument('--gff', help='Single combined GFF file')
    parser.add_argument('--clusters', required=True, help='Clusters file or gene_presence_absence.csv')
    parser.add_argument('--outdir', required=True, help='Output directory for inspection files')
    parser.add_argument('--k', type=int, default=3, help='Number of neighbors upstream/downstream to include')
    parser.add_argument('--rel_delta', type=float, default=0.05, help='Relative index delta for conserved position')
    parser.add_argument('--id-map', help='Optional id-map file mapping cluster token -> gff id')
    parser.add_argument('--exclude-self', action='store_true', help='Exclude the focal family from neighbor-family aggregation and heatmaps')
    parser.add_argument('--plot', action='store_true', help='Also produce PNG visualizations (requires matplotlib)')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()

    # locate the location_tracking module file relative to this script
    script_dir = Path(__file__).resolve().parent
    lt_path = script_dir.joinpath('..', 'src', 'PyamilySeq', 'location_tracking.py').resolve()
    if not lt_path.exists():
        print(f'ERROR: location_tracking.py not found at expected path: {lt_path}', file=sys.stderr)
        sys.exit(1)

    # load the module functions using runpy
    mod = runpy.run_path(str(lt_path))
    parse_gff = mod.get('parse_gff')
    parse_gff_dir = mod.get('parse_gff_dir')
    parse_clusters = mod.get('parse_clusters')
    resolve_genome_qualifiers = mod.get('resolve_genome_qualifiers')
    filter_singletons = mod.get('filter_singletons')
    map_family_members_to_features = mod.get('map_family_members_to_features')
    compute_neighbors_for_contigs = mod.get('compute_neighbors_for_contigs')
    compute_family_metrics = mod.get('compute_family_metrics')

    # sanity check functions
    for fn in ('parse_gff', 'parse_gff_dir', 'parse_clusters', 'map_family_members_to_features', 'compute_family_metrics'):
        if fn not in mod:
            print(f'ERROR: {fn} not found in location_tracking module', file=sys.stderr)
            sys.exit(1)

    # parse GFF(s)
    if args.gff_dir:
        contigs, id_to_feature, gff_genomes = parse_gff_dir(args.gff_dir)
    else:
        contigs, id_to_feature = parse_gff(args.gff)
        gff_genomes = []

    # parse clusters and resolve qualifiers if possible
    families = parse_clusters(args.clusters)
    if gff_genomes:
        families = resolve_genome_qualifiers(families, gff_genomes)

    # Build reverse mapping: canonical gene id -> set of family ids
    # Use helpers from location_tracking if available
    normalise_member_to_gid = mod.get('normalise_member_to_gid')
    sanitize_id = mod.get('sanitize_id')
    gid_to_families = defaultdict(set)
    for fam, members in families.items():
        for m in members:
            # m may be 'genome|gene' or similar
            try:
                gid = normalise_member_to_gid(m) if normalise_member_to_gid else (m.split('|')[-1] if '|' in m else m)
            except Exception:
                gid = m
            if sanitize_id:
                gid = sanitize_id(gid)
            gid_to_families[gid].add(fam)

    if args.family not in families:
        print(f'ERROR: family {args.family} not found in clusters (available families: {len(families)})', file=sys.stderr)
        # provide partial hint: maybe family is present without genome qualifier
        # attempt to find close names
        candidates = [f for f in families.keys() if args.family in f]
        if candidates:
            print('Did you mean one of:', file=sys.stderr)
            for c in candidates[:20]:
                print('  ', c, file=sys.stderr)
        sys.exit(1)

    # restrict families dict to just the requested family
    fam_subset = {args.family: families[args.family]}

    # map family members to features
    family_occ, unmapped = map_family_members_to_features(fam_subset, id_to_feature)

    # compute neighbors on contigs
    compute_neighbors_for_contigs(contigs, k=args.k)

    # compute family metrics for this family
    occ_rows, sum_rows = compute_family_metrics(family_occ, k=args.k, rel_delta=args.rel_delta)

    # prepare outputs
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    occ_out = outdir / f'{args.family}_occurrences.csv'
    sum_out = outdir / f'{args.family}_summary.csv'
    ctx_out = outdir / f'{args.family}_neighbors_contexts.tsv'

    # write occurrences and summary
    occ_fieldnames = ['family_id', 'gene_id', 'raw_id', 'genome', 'contig', 'start', 'end', 'strand', 'ordinal_index', 'relative_index', 'contig_gene_count', 'neighbors_up', 'neighbors_down']
    write_csv_rows(str(occ_out), occ_rows, occ_fieldnames)
    # summary may have different fieldnames
    if sum_rows:
        sum_fieldnames = list(sum_rows[0].keys())
        write_csv_rows(str(sum_out), sum_rows, sum_fieldnames)

    # neighbors contexts: one line per occurrence with human-friendly neighbor listing
    with open(ctx_out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['family_id', 'gene_id', 'genome', 'contig', 'start', 'end', 'strand', 'neighbors_up', 'neighbors_down', 'neighbor_count', 'self_neighbor', 'self_reason', 'self_positions'])
        # Also create neighbor-family mappings and aggregate counts
        neighbor_family_counter = Counter()
        occurrences_with_family = defaultdict(int)
        neighbor_family_by_genome = defaultdict(lambda: defaultdict(int))
        genomes_set = set()
        total_occurrences = len(occ_rows)
        total_neighbor_tokens = 0
        for r in occ_rows:
            up = r.get('neighbors_up') or ''
            down = r.get('neighbors_down') or ''
            genome = r.get('genome') or 'unknown'
            genomes_set.add(genome)
            ncount = 0
            # track which neighbor families appear at least once in this occurrence
            occ_fams = set()
            up_fams = []
            down_fams = []
            # track self-trigger positions for exact token equality and cluster-family matches
            self_positions = []
            self_reasons = set()

            if up:
                tokens = [t for t in up.split(';') if t]
                ncount += len(tokens)
                for t in tokens:
                    total_neighbor_tokens += 1
                    mapped = gid_to_families.get(t)
                    if not mapped and sanitize_id:
                        mapped = gid_to_families.get(sanitize_id(t))
                    if mapped:
                        # always record family membership for detection of self-neighbors
                        for mf in sorted(mapped):
                            occ_fams.add(mf)
                            # only count in aggregations if not excluding self or mf != focal
                            if not (args.exclude_self and mf == args.family):
                                neighbor_family_counter[mf] += 1
                                neighbor_family_by_genome[mf][genome] += 1
                            # if the mapped family equals focal family, note cluster-based self match
                            if mf == r.get('family_id'):
                                self_reasons.add('SELF_CLUSTER')

                        # check exact token equality with focal gene id
                        focal_id = r.get('gene_id')
                        if focal_id and t == focal_id:
                            self_reasons.add('SELF_EXACT')
                            self_positions.append(f'up{len(up_fams)+1}')
                            # append placeholder for this neighbor position so offsets align, but do not count the family
                            up_fams.append('NA')
                        else:
                            up_fams.append('|'.join(sorted(mapped)))
                    else:
                        # use explicit placeholder for unmapped neighbor tokens
                        up_fams.append('NA')

            if down:
                tokens = [t for t in down.split(';') if t]
                ncount += len(tokens)
                for t in tokens:
                    total_neighbor_tokens += 1
                    mapped = gid_to_families.get(t)
                    if not mapped and sanitize_id:
                        mapped = gid_to_families.get(sanitize_id(t))
                    if mapped:
                        for mf in sorted(mapped):
                            occ_fams.add(mf)
                            if not (args.exclude_self and mf == args.family):
                                neighbor_family_counter[mf] += 1
                                neighbor_family_by_genome[mf][genome] += 1
                            if mf == r.get('family_id'):
                                self_reasons.add('SELF_CLUSTER')

                        focal_id = r.get('gene_id')
                        if focal_id and t == focal_id:
                            self_reasons.add('SELF_EXACT')
                            # down positions count from 1; append placeholder to preserve offset alignment
                            self_positions.append(f'down{len(down_fams)+1}')
                            down_fams.append('NA')
                        else:
                            down_fams.append('|'.join(sorted(mapped)))
                    else:
                        # explicit placeholder for unmapped neighbor tokens
                        down_fams.append('NA')
            # determine self-neighbor flag and increment occurrences_with_family
            self_neighbor_flag = False
            filtered_occ_fams = set()
            for mf in occ_fams:
                if mf == args.family:
                    self_neighbor_flag = True
                    if args.exclude_self:
                        # do not include this family in aggregated occurrence counts
                        continue
                filtered_occ_fams.add(mf)
            for mf in filtered_occ_fams:
                occurrences_with_family[mf] += 1
            # join family lists aligned to neighbor tokens (semicolon-separated)
            up_fams_s = ';'.join(up_fams) if up_fams else ''
            down_fams_s = ';'.join(down_fams) if down_fams else ''
            # decide self_reason string
            if self_reasons:
                # prefer EXACT over CLUSTER if both present
                if 'SELF_EXACT' in self_reasons:
                    reason = 'SELF_EXACT'
                else:
                    reason = 'SELF_CLUSTER'
            else:
                reason = ''
            positions_s = ';'.join(self_positions) if self_positions else ''
            w.writerow([r.get('family_id'), r.get('gene_id'), r.get('genome'), r.get('contig'), r.get('start'), r.get('end'), r.get('strand'), up_fams_s, down_fams_s, ncount, int(self_neighbor_flag), reason, positions_s])

    # Write aggregated neighbor-family summary CSV (clearer columns)
    summary_out = outdir / f'{args.family}_neighbor_family_summary.csv'
    with open(summary_out, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['neighbor_family', 'count', 'occurrences_with_family', 'fraction_of_neighbor_tokens', 'fraction_of_occurrences'])
        for fam_id, cnt in neighbor_family_counter.most_common():
            occ_with = occurrences_with_family.get(fam_id, 0)
            frac_tokens = cnt / total_neighbor_tokens if total_neighbor_tokens else 0.0
            frac_occ = occ_with / total_occurrences if total_occurrences else 0.0
            writer.writerow([fam_id, cnt, occ_with, f'{frac_tokens:.6f}', f'{frac_occ:.6f}'])

    # Write heatmap CSV: rows = neighbor_family, cols = genome, values = counts
    genomes = sorted(genomes_set)
    heatmap_out = outdir / f'{args.family}_neighbor_family_heatmap.csv'
    with open(heatmap_out, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(['neighbor_family'] + genomes)
        for fam_id in sorted(neighbor_family_by_genome.keys()):
            row = [fam_id] + [neighbor_family_by_genome[fam_id].get(g, 0) for g in genomes]
            writer.writerow(row)

    # Write a simple self-contained HTML heatmap for quick visual inspection
    html_out = outdir / f'{args.family}_neighbor_family_heatmap.html'
    max_count = 0
    for fam_counts in neighbor_family_by_genome.values():
        for v in fam_counts.values():
            if v > max_count:
                max_count = v
    # safe guard
    if max_count == 0:
        max_count = 1
    with open(html_out, 'w') as fh:
        fh.write('<!doctype html>\n<html><head><meta charset="utf-8"><title>Neighbor family heatmap</title>\n')
        fh.write('<style>table{border-collapse:collapse;font-family:Arial,Helvetica,sans-serif}td,th{border:1px solid #ccc;padding:4px;text-align:center}\n')
        fh.write('.label{font-size:12px;padding:6px;text-align:left}</style>\n')
        fh.write('</head><body>\n')
        fh.write(f'<h2>Neighbor family heatmap for {args.family}</h2>\n')
        fh.write('<p>Cell color intensity corresponds to count of neighbor tokens for this neighbor family in the given genome.</p>\n')
        # summary table
        fh.write('<h3>Summary (top neighbor families)</h3>\n')
        fh.write('<table><tr><th>neighbor_family</th><th>count</th><th>occurrences_with_family</th><th>fraction_of_neighbor_tokens</th><th>fraction_of_occurrences</th></tr>\n')
        for fam_id, cnt in neighbor_family_counter.most_common(50):
            occ_with = occurrences_with_family.get(fam_id, 0)
            frac_tokens = cnt / total_neighbor_tokens if total_neighbor_tokens else 0.0
            frac_occ = occ_with / total_occurrences if total_occurrences else 0.0
            fh.write(f'<tr><td class="label">{fam_id}</td><td>{cnt}</td><td>{occ_with}</td><td>{frac_tokens:.6f}</td><td>{frac_occ:.6f}</td></tr>\n')
        fh.write('</table>\n')
        # heatmap
        fh.write('<h3>Heatmap (neighbor families x genomes)</h3>\n')
        fh.write('<table>\n<tr><th>neighbor_family</th>')
        for g in genomes:
            fh.write(f'<th>{g}</th>')
        fh.write('</tr>\n')
        for fam_id in sorted(neighbor_family_by_genome.keys()):
            fh.write(f'<tr><td class="label">{fam_id}</td>')
            for g in genomes:
                v = neighbor_family_by_genome[fam_id].get(g, 0)
                # color ramp from white to red based on v/max_count
                ratio = v / max_count
                red = int(255 * ratio)
                green = 255 - int(150 * ratio)
                blue = 255 - int(200 * ratio)
                color = f'#{red:02x}{green:02x}{blue:02x}'
                fh.write(f'<td style="background:{color}">{v}</td>')
            fh.write('</tr>\n')
        fh.write('</table>\n')
        fh.write('</body></html>\n')

    # Optionally produce WMG-style PNG visualizations using WMG.make_heatmap and make_barplots
    if args.plot:
        try:
            # build WMG-style records from occ_rows
            def build_wmg_records(occ_rows, k):
                # offsets: -k..0..k (0 reserved for focal gene)
                offsets = list(range(-k, k + 1))
                records = []
                for i, r in enumerate(occ_rows):
                    rec = {
                        'genome': r.get('genome') or 'unknown',
                        'gff_file': r.get('contig'),
                        'contig': r.get('contig'),
                        'focal_geneid': r.get('gene_id'),
                        'focal_family': r.get('family_id'),
                        'focal_index': i,
                        'focal_length': None,
                        'neighbors': {}
                    }
                    up_tokens = [t for t in (r.get('neighbors_up') or '').split(';') if t]
                    down_tokens = [t for t in (r.get('neighbors_down') or '').split(';') if t]

                    # ensure lengths: pad with 'NA' to k
                    while len(up_tokens) < k:
                        up_tokens.append('NA')
                    while len(down_tokens) < k:
                        down_tokens.append('NA')

                    # insert focal gene at offset 0 so visualizations can center on focal family
                    focal_token = r.get('gene_id') if r.get('gene_id') else None
                    focal_family = r.get('family_id') if r.get('family_id') else None
                    rec['neighbors'][0] = {
                        'geneid': focal_token,
                        'geneid_raw': focal_token,
                        'family': focal_family,
                        'length': None,
                        'start': r.get('start'),
                        'end': r.get('end'),
                        'strand': r.get('strand')
                    }

                    # map tokens: up nearest-first corresponds to offsets -1, -2, ...
                    for idx, tok in enumerate(up_tokens):
                        off = -(idx + 1)
                        token = tok if tok and tok != 'NA' else None
                        # If the neighbor token exactly equals the focal token, treat it as missing
                        if token and focal_token and token == focal_token:
                            fam = None
                            geneid_val = None
                        else:
                            fams = gid_to_families.get(token) if token else None
                            fam = sorted(fams)[0] if fams else None
                            geneid_val = token
                        rec['neighbors'][off] = {
                            'geneid': geneid_val,
                            'geneid_raw': token,
                            'family': fam,
                            'length': None,
                            'start': None,
                            'end': None,
                            'strand': None
                        }

                    # map downstream tokens to offsets 1..k
                    for idx, tok in enumerate(down_tokens):
                        off = idx + 1
                        token = tok if tok and tok != 'NA' else None
                        # If the neighbor token exactly equals the focal token, treat it as missing
                        if token and focal_token and token == focal_token:
                            fam = None
                            geneid_val = None
                        else:
                            fams = gid_to_families.get(token) if token else None
                            fam = sorted(fams)[0] if fams else None
                            geneid_val = token
                        rec['neighbors'][off] = {
                            'geneid': geneid_val,
                            'geneid_raw': token,
                            'family': fam,
                            'length': None,
                            'start': None,
                            'end': None,
                            'strand': None
                        }

                    records.append(rec)
                return records, offsets

            records_wmg, offsets_wmg = build_wmg_records(occ_rows, args.k)
            # try to load WMG plotting functions
            wmg_path = os.path.join(os.path.dirname(__file__), '..', '..', 'WMG', 'src', 'WMG', 'WMG.py')
            # fallback to repo path
            if not os.path.exists(wmg_path):
                wmg_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'WMG', 'src', 'WMG', 'WMG.py')
            if os.path.exists(wmg_path):
                wmg_mod = runpy.run_path(wmg_path)
                make_heatmap = wmg_mod.get('make_heatmap')
                make_barplots = wmg_mod.get('make_barplots')
                if make_heatmap:
                    try:
                        ok = make_heatmap(records_wmg, offsets_wmg, str(outdir), f"{args.family}_wmg_neighborhood")
                        if ok:
                            print('WMG-style heatmap generated.')
                    except Exception as e:
                        print('WMG heatmap generation failed:', e)
                if make_barplots:
                    try:
                        stats = wmg_mod.get('detailed_offset_stats')({off: [rec['neighbors'][off].get('family') for rec in records_wmg] for off in offsets_wmg})
                        ok2 = make_barplots(stats, str(outdir), f"{args.family}_wmg_neighborhood")
                        if ok2:
                            print('WMG-style barplots generated.')
                    except Exception as e:
                        print('WMG barplots generation failed:', e)
            else:
                print('WMG module not found; skipping WMG-style PNG visualizations')
        except Exception as e:
            print('Error while producing WMG-style visualizations:', e)

    # Print a short report
    print(f'Inspection for family: {args.family}')
    if sum_rows:
        s = sum_rows[0]
        print('\nSummary metrics:')
        for k, v in s.items():
            print(f'  {k}: {v}')
    print('\nOccurrences written to:', occ_out)
    print('Summary written to:', sum_out)
    print('Neighbor contexts written to:', ctx_out)

    if unmapped and unmapped.get(args.family):
        print('\nUnmapped cluster members for this family (sample):')
        for token in unmapped.get(args.family, [])[:50]:
            print('  ', token)

    # Suggest extra improvements/features
    print('\nSuggested extra features:')
    suggestions = [
        'Map neighbor gene IDs to family IDs and report neighbor-family conservation (requires reverse mapping family->members)',
        'If genome FASTA files available, extract sequences for each occurrence and output a multi-FASTA for the family',
        'Produce an interactive HTML synteny visualization (SVG or D3) with links to coordinates',
        'Compute pairwise nucleotide similarity between occurrences (e.g., using sowhat like BLAST or simple pairwise alignment) to detect paralogs vs orthologs',
        'Report neighbor orientation consistency and optionally plot distributions across genomes',
        'Provide an option to collapse very similar neighbor tokens (aliasing / synonyms) via an id-map file',
    ]
    for s in suggestions:
        print(' -', s)


if __name__ == '__main__':
    main()

