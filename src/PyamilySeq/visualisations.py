import html
import json
import os
import random
from collections import Counter, OrderedDict


PALETTE = {
    'First': '#2f7ed8',
    'extended': '#50b87f',
    'combined': '#f2a93b',
    'Second': '#9b6bd3',
    'only': '#d95763',
    'default': '#607d8b',
}


def _escape(value):
    return html.escape(str(value), quote=True)


def _write(path, text):
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write(text)


def _category_colour(category):
    if category.startswith('First'):
        return PALETTE['First']
    if category.startswith('extended'):
        return PALETTE['extended']
    if category.startswith('combined'):
        return PALETTE['combined']
    if category.startswith('only_Second'):
        return PALETTE['only']
    if category.startswith('Second'):
        return PALETTE['Second']
    return PALETTE['default']


def _nice_label(label):
    return str(label).replace('_', ' ')


def _compact_label(label, max_chars=18):
    value = str(label)
    if len(value) <= max_chars:
        return value
    suffix_length = max_chars - 9
    return f'{value[:6]}...{value[-suffix_length:]}'


def _svg_header(width, height):
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img">',
        '<style>'
        'text{font-family:Arial,Helvetica,sans-serif;fill:#1f2933}'
        '.title{font-size:24px;font-weight:700}'
        '.subtitle{font-size:13px;fill:#52606d}'
        '.axis{font-size:12px;fill:#52606d}'
        '.small{font-size:10px;fill:#52606d}'
        '.label{font-size:12px}'
        '</style>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
    ]


def _render_category_svg(category_counts, output_path, mode_label, taxa_label):
    width = 980
    row_h = 34
    top = 82
    left = 230
    right = 70
    bottom = 40
    rows = max(1, len(category_counts))
    height = top + rows * row_h + bottom
    max_count = max(category_counts.values()) if category_counts else 1
    chart_w = width - left - right

    svg = _svg_header(width, height)
    svg.append(f'<text class="title" x="28" y="38">{_escape(mode_label)} Pangenome Categories</text>')
    svg.append(f'<text class="subtitle" x="28" y="62">Gene family counts by {_escape(taxa_label)}-presence category</text>')

    if not category_counts:
        svg.append('<text class="label" x="28" y="110">No families available to plot.</text>')
    else:
        for idx, (category, count) in enumerate(category_counts.items()):
            y = top + idx * row_h
            bar_w = 0 if max_count == 0 else int(chart_w * (count / max_count))
            colour = _category_colour(category)
            svg.append(f'<text class="label" x="28" y="{y + 20}">{_escape(_nice_label(category))}</text>')
            svg.append(f'<rect x="{left}" y="{y + 6}" width="{bar_w}" height="20" rx="3" fill="{colour}"/>')
            svg.append(f'<text class="axis" x="{left + bar_w + 8}" y="{y + 21}">{count}</text>')

    svg.append('</svg>\n')
    _write(output_path, '\n'.join(svg))


def _histogram_bins(values, taxa_total):
    if taxa_total <= 60:
        counts = Counter(values)
        return [(str(i), counts.get(i, 0)) for i in range(1, taxa_total + 1)]

    bin_count = 20
    width = max(1, (taxa_total + bin_count - 1) // bin_count)
    bins = OrderedDict()
    start = 1
    while start <= taxa_total:
        end = min(taxa_total, start + width - 1)
        label = f'{start}-{end}' if start != end else str(start)
        bins[label] = 0
        start = end + 1
    for value in values:
        idx = min((max(1, value) - 1) // width, len(bins) - 1)
        label = list(bins.keys())[idx]
        bins[label] += 1
    return list(bins.items())


def _render_histogram_svg(cluster_taxa, output_path, mode_label, taxa_label):
    taxa_counts = [len(set(values)) for values in cluster_taxa.values()]
    taxa_total = max(taxa_counts) if taxa_counts else 1
    bins = _histogram_bins(taxa_counts, taxa_total)

    width = 980
    height = 430
    left = 58
    right = 28
    top = 84
    bottom = 78
    chart_w = width - left - right
    chart_h = height - top - bottom
    max_count = max((count for _, count in bins), default=1)
    bar_gap = 3
    bar_w = max(2, (chart_w - bar_gap * max(0, len(bins) - 1)) / max(1, len(bins)))

    svg = _svg_header(width, height)
    svg.append(f'<text class="title" x="28" y="38">{_escape(mode_label)} Family Prevalence</text>')
    svg.append(f'<text class="subtitle" x="28" y="62">Number of families observed across {_escape(taxa_label)}</text>')
    svg.append(f'<line x1="{left}" y1="{top + chart_h}" x2="{width - right}" y2="{top + chart_h}" stroke="#9aa5b1"/>')
    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + chart_h}" stroke="#9aa5b1"/>')

    for idx, (label, count) in enumerate(bins):
        x = left + idx * (bar_w + bar_gap)
        h = 0 if max_count == 0 else chart_h * (count / max_count)
        y = top + chart_h - h
        svg.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_w:.2f}" height="{h:.2f}" fill="#2f7ed8">')
        svg.append(f'<title>{_escape(label)} {_escape(taxa_label)}: {count} families</title></rect>')
        if len(bins) <= 30 or idx % max(1, len(bins) // 20) == 0:
            svg.append(f'<text class="small" x="{x + bar_w / 2:.2f}" y="{top + chart_h + 18}" text-anchor="middle" transform="rotate(45 {x + bar_w / 2:.2f},{top + chart_h + 18})">{_escape(label)}</text>')

    svg.append(f'<text class="axis" x="{left}" y="{height - 18}">{_escape(taxa_label.title())} per family</text>')
    svg.append(f'<text class="axis" x="20" y="{top - 8}">Families</text>')
    svg.append('</svg>\n')
    _write(output_path, '\n'.join(svg))


def _cluster_sort_items(cluster_taxa, cluster_sequences):
    return sorted(
        cluster_taxa.items(),
        key=lambda item: (len(set(item[1])), len(cluster_sequences.get(item[0], [])), str(item[0])),
        reverse=True,
    )


def _render_heatmap_svg(cluster_taxa, cluster_sequences, taxa, output_path, mode_label, taxa_label, limit):
    selected = _cluster_sort_items(cluster_taxa, cluster_sequences)[:limit]
    taxa = list(taxa)
    if len(taxa) <= 15:
        cell_w = 64
    elif len(taxa) <= 30:
        cell_w = 34
    elif len(taxa) <= 80:
        cell_w = 14
    else:
        cell_w = max(4, min(10, 1100 / max(1, len(taxa))))
    row_h = 14
    left = 190
    top = 190 if len(taxa) <= 80 else 128
    right = 30
    bottom = 86
    width = int(max(980, left + len(taxa) * cell_w + right))
    height = int(top + max(1, len(selected)) * row_h + bottom)

    svg = _svg_header(width, height)
    svg.append(f'<text class="title" x="28" y="38">{_escape(mode_label)} Presence/Absence Heatmap</text>')
    svg.append(f'<text class="subtitle" x="28" y="62">Top {len(selected)} families sorted by {_escape(taxa_label)} prevalence and sequence count</text>')
    if len(cluster_taxa) > limit:
        svg.append(f'<text class="subtitle" x="28" y="82">Showing {len(selected)} of {len(cluster_taxa)} total families.</text>')

    for idx, taxon in enumerate(taxa):
        x = left + idx * cell_w + cell_w / 2
        if len(taxa) <= 80 or idx % max(1, len(taxa) // 80) == 0:
            visible_taxon = _compact_label(taxon)
            svg.append(
                f'<text class="small" x="{x:.2f}" y="{top - 10}" text-anchor="start" '
                f'transform="rotate(-60 {x:.2f},{top - 10})"><title>{_escape(taxon)}</title>'
                f'{_escape(visible_taxon)}</text>'
            )

    for row_idx, (cluster, present_taxa) in enumerate(selected):
        y = top + row_idx * row_h
        present = set(present_taxa)
        seq_count = len(cluster_sequences.get(cluster, []))
        cluster_label = str(cluster) if str(cluster).startswith('group_') else f'group_{cluster}'
        visible_family = cluster_label if len(cluster_label) <= 24 else cluster_label[:23] + '.'
        svg.append(
            f'<text class="small" x="{left - 8}" y="{y + 10}" text-anchor="end">'
            f'<title>{_escape(cluster_label)}</title>{_escape(visible_family)}</text>'
        )
        for col_idx, taxon in enumerate(taxa):
            x = left + col_idx * cell_w
            fill = '#2f7ed8' if taxon in present else '#e6edf3'
            state = 'present' if taxon in present else 'absent'
            tooltip = f'{cluster_label}; {taxon}; {state}; {seq_count} sequences'
            svg.append(f'<rect x="{x:.2f}" y="{y}" width="{cell_w:.2f}" height="{row_h - 1}" fill="{fill}">')
            svg.append(f'<title>{_escape(tooltip)}</title></rect>')

    svg.append(f'<text class="axis" x="{left}" y="{height - 24}">Columns: {_escape(taxa_label)}. Filled cells indicate presence.</text>')
    svg.append('</svg>\n')
    _write(output_path, '\n'.join(svg))


def _write_data_tables(output_dir, category_counts, cluster_taxa, cluster_sequences):
    category_path = os.path.join(output_dir, 'category_counts.tsv')
    with open(category_path, 'w', encoding='utf-8') as handle:
        handle.write('category\tfamily_count\n')
        for category, count in category_counts.items():
            handle.write(f'{category}\t{count}\n')

    prevalence_path = os.path.join(output_dir, 'family_prevalence.tsv')
    with open(prevalence_path, 'w', encoding='utf-8') as handle:
        handle.write('family_id\ttaxa_count\tsequence_count\n')
        for cluster, taxa in _cluster_sort_items(cluster_taxa, cluster_sequences):
            family_id = str(cluster) if str(cluster).startswith('group_') else f'group_{cluster}'
            handle.write(f'{family_id}\t{len(set(taxa))}\t{len(cluster_sequences.get(cluster, []))}\n')


def _render_index(output_dir, mode_label, taxa_label, taxa_total, family_total, sequence_total, limit):
    path = os.path.join(output_dir, 'pangenome_visualisations.html')
    html_text = f'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>{_escape(mode_label)} Pangenome Visualisations</title>
  <style>
    body {{ font-family: Arial, Helvetica, sans-serif; margin: 28px; color: #1f2933; background: #f7f9fb; }}
    main {{ max-width: 1100px; margin: 0 auto; }}
    section {{ background: white; border: 1px solid #d9e2ec; border-radius: 6px; padding: 18px; margin: 18px 0; }}
    img {{ max-width: 100%; height: auto; border: 1px solid #e6edf3; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 12px; }}
    .metric {{ background: #eef4fb; border-radius: 6px; padding: 12px; }}
    .metric strong {{ display: block; font-size: 22px; }}
    a {{ color: #1d5f99; }}
  </style>
</head>
<body>
<main>
  <h1>{_escape(mode_label)} Pangenome Visualisations</h1>
  <section class="metrics">
    <div class="metric"><strong>{taxa_total}</strong>{_escape(taxa_label.title())}</div>
    <div class="metric"><strong>{family_total}</strong>Gene families</div>
    <div class="metric"><strong>{sequence_total}</strong>Sequences in families</div>
    <div class="metric"><strong>{limit}</strong>Heatmap row limit</div>
  </section>
  <section>
    <h2>Category Counts</h2>
    <img src="core_categories.svg" alt="Pangenome category counts">
  </section>
  <section>
    <h2>Family Prevalence</h2>
    <img src="family_prevalence_histogram.svg" alt="Family prevalence histogram">
  </section>
  <section>
    <h2>Presence/Absence Heatmap</h2>
    <img src="presence_absence_heatmap.svg" alt="Presence absence heatmap">
  </section>
  <section>
    <h2>Underlying Tables</h2>
    <p><a href="category_counts.tsv">category_counts.tsv</a> and <a href="family_prevalence.tsv">family_prevalence.tsv</a></p>
  </section>
</main>
</body>
</html>
'''
    _write(path, html_text)


def write_pangenome_visualisations(options, taxa_dict, cores, cluster_taxa_sorted,
                                   cluster_sequences_sorted, mode_label='Species',
                                   taxa_label='genomes'):
    if not getattr(options, 'visualise', False):
        return None

    output_dir = os.path.join(os.path.abspath(options.output_dir), 'Visualisations')
    os.makedirs(output_dir, exist_ok=True)

    category_counts = OrderedDict((key, len(values)) for key, values in cores.items())
    taxa = sorted(taxa_dict.keys())
    limit = int(getattr(options, 'visualise_heatmap_limit', 200) or 200)
    limit = max(1, limit)
    family_total = len(cluster_taxa_sorted)
    sequence_total = sum(len(values) for values in cluster_sequences_sorted.values())

    _render_category_svg(
        category_counts,
        os.path.join(output_dir, 'core_categories.svg'),
        mode_label,
        taxa_label,
    )
    _render_histogram_svg(
        cluster_taxa_sorted,
        os.path.join(output_dir, 'family_prevalence_histogram.svg'),
        mode_label,
        taxa_label,
    )
    _render_heatmap_svg(
        cluster_taxa_sorted,
        cluster_sequences_sorted,
        taxa,
        os.path.join(output_dir, 'presence_absence_heatmap.svg'),
        mode_label,
        taxa_label,
        limit,
    )
    _write_data_tables(output_dir, category_counts, cluster_taxa_sorted, cluster_sequences_sorted)
    _render_index(output_dir, mode_label, taxa_label, len(taxa), family_total, sequence_total, limit)

    return output_dir


def _copy_number_matrix(ledger, family_ids, taxa, mode):
    resolver = ledger.resolver.genome if mode == 'Species' else ledger.resolver.taxon
    rows = []
    for family in ledger.sorted_families(family_ids):
        counts = Counter(resolver(member) for member in family.members)
        rows.append((family.family_id, [counts.get(taxon, 0) for taxon in taxa]))
    return rows


def _render_copy_number_svg(rows, taxa, path, limit):
    selected = rows[:limit]
    if len(taxa) <= 15:
        cell_w = 64
    elif len(taxa) <= 30:
        cell_w = 34
    elif len(taxa) <= 80:
        cell_w = 14
    else:
        cell_w = max(4, min(10, 1100 / max(1, len(taxa))))
    row_h = 14
    left, top, right, bottom = 190, (190 if len(taxa) <= 80 else 128), 30, 86
    width = int(max(980, left + len(taxa) * cell_w + right))
    height = int(top + max(1, len(selected)) * row_h + bottom)
    colours = {0: '#edf1f4', 1: '#2f7ed8', 2: '#50b87f', 3: '#f2a93b'}
    svg = _svg_header(width, height)
    svg.append('<text class="title" x="28" y="38">Copy-Number Heatmap</text>')
    svg.append('<text class="subtitle" x="28" y="62">0 copies grey, 1 blue, 2 green, 3 orange, 4+ red</text>')
    if len(rows) > limit:
        svg.append(f'<text class="subtitle" x="28" y="82">Showing {len(selected)} of {len(rows)} families.</text>')
    for index, taxon in enumerate(taxa):
        x = left + index * cell_w + cell_w / 2
        if len(taxa) <= 80 or index % max(1, len(taxa) // 80) == 0:
            visible_taxon = _compact_label(taxon)
            svg.append(
                f'<text class="small" x="{x:.2f}" y="{top - 10}" text-anchor="start" '
                f'transform="rotate(-60 {x:.2f},{top - 10})"><title>{_escape(taxon)}</title>'
                f'{_escape(visible_taxon)}</text>'
            )
    for row_index, (family_id, counts) in enumerate(selected):
        y = top + row_index * row_h
        visible_family = family_id if len(family_id) <= 24 else family_id[:23] + '.'
        svg.append(
            f'<text class="small" x="{left - 8}" y="{y + 10}" text-anchor="end">'
            f'<title>{_escape(family_id)}</title>{_escape(visible_family)}</text>'
        )
        for column_index, count in enumerate(counts):
            x = left + column_index * cell_w
            fill = '#d95763' if count >= 4 else colours[count]
            svg.append(
                f'<rect x="{x:.2f}" y="{y}" width="{cell_w:.2f}" height="{row_h - 1}" fill="{fill}">'
                f'<title>{_escape(family_id)}; {_escape(taxa[column_index])}; copies={count}</title></rect>'
            )
    svg.append('</svg>\n')
    _write(path, '\n'.join(svg))


def _write_copy_number_table(path, rows, taxa):
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('family_id\t' + '\t'.join(taxa) + '\n')
        for family_id, counts in rows:
            handle.write(f'{family_id}\t' + '\t'.join(str(value) for value in counts) + '\n')


def _accumulation_rows(rows, taxa, permutations=100, seed=1):
    if not taxa:
        return []
    presence = [set(taxon for taxon, count in zip(taxa, counts) if count) for _, counts in rows]
    rng = random.Random(seed)
    pangenome = [[] for _ in taxa]
    core = [[] for _ in taxa]
    permutations = max(1, min(permutations, 1000))
    for _ in range(permutations):
        order = list(taxa)
        rng.shuffle(order)
        selected = set()
        for index, taxon in enumerate(order):
            selected.add(taxon)
            pangenome[index].append(sum(bool(family_taxa & selected) for family_taxa in presence))
            core[index].append(sum(selected <= family_taxa for family_taxa in presence))
    result = []
    for index in range(len(taxa)):
        pan_values = sorted(pangenome[index])
        core_values = sorted(core[index])
        low_index = int(0.025 * (len(pan_values) - 1))
        high_index = int(0.975 * (len(pan_values) - 1))
        result.append({
            'taxa_sampled': index + 1,
            'pangenome_mean': sum(pan_values) / len(pan_values),
            'pangenome_low': pan_values[low_index],
            'pangenome_high': pan_values[high_index],
            'core_mean': sum(core_values) / len(core_values),
            'core_low': core_values[low_index],
            'core_high': core_values[high_index],
        })
    return result


def _write_accumulation_table(path, rows):
    fields = [
        'taxa_sampled', 'pangenome_mean', 'pangenome_low', 'pangenome_high',
        'core_mean', 'core_low', 'core_high',
    ]
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('\t'.join(fields) + '\n')
        for row in rows:
            handle.write('\t'.join(
                f'{row[field]:.6f}' if isinstance(row[field], float) else str(row[field])
                for field in fields
            ) + '\n')


def _polyline(values, left, top, chart_w, chart_h, max_value):
    if not values:
        return ''
    denominator = max(1, len(values) - 1)
    return ' '.join(
        f'{left + index * chart_w / denominator:.2f},{top + chart_h - value * chart_h / max_value:.2f}'
        for index, value in enumerate(values)
    )


def _render_accumulation_svg(path, rows):
    width, height = 980, 470
    left, right, top, bottom = 72, 32, 82, 64
    chart_w, chart_h = width - left - right, height - top - bottom
    pan = [row['pangenome_mean'] for row in rows]
    core = [row['core_mean'] for row in rows]
    max_value = max(pan + core + [1])
    svg = _svg_header(width, height)
    svg.append('<text class="title" x="28" y="38">Pangenome Accumulation</text>')
    svg.append('<text class="subtitle" x="28" y="62">Deterministic permutation mean (seed 1); pangenome and strict-core trajectories</text>')
    svg.append(f'<line x1="{left}" y1="{top + chart_h}" x2="{width - right}" y2="{top + chart_h}" stroke="#9aa5b1"/>')
    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + chart_h}" stroke="#9aa5b1"/>')
    for tick in range(5):
        value = max_value * tick / 4
        y = top + chart_h - chart_h * tick / 4
        svg.append(f'<line x1="{left - 4}" y1="{y:.2f}" x2="{width - right}" y2="{y:.2f}" stroke="#e6edf3"/>')
        svg.append(f'<text class="small" x="{left - 8}" y="{y + 4:.2f}" text-anchor="end">{value:.0f}</text>')
    if rows:
        svg.append(f'<polyline points="{_polyline(pan, left, top, chart_w, chart_h, max_value)}" fill="none" stroke="#2f7ed8" stroke-width="3"/>')
        svg.append(f'<polyline points="{_polyline(core, left, top, chart_w, chart_h, max_value)}" fill="none" stroke="#d95763" stroke-width="3"/>')
    svg.append('<rect x="628" y="50" width="14" height="4" fill="#2f7ed8"/><text class="axis" x="650" y="57">Pangenome</text>')
    svg.append('<rect x="760" y="50" width="14" height="4" fill="#d95763"/><text class="axis" x="782" y="57">Strict core</text>')
    if rows:
        svg.append(f'<text class="small" x="{left}" y="{top + chart_h + 18}">1</text>')
        svg.append(f'<text class="small" x="{width - right}" y="{top + chart_h + 18}" text-anchor="end">{len(rows)}</text>')
    svg.append(f'<text class="axis" x="{left + chart_w / 2}" y="{height - 20}" text-anchor="middle">Taxa sampled</text>')
    svg.append('<text class="axis" x="18" y="74">Families</text>')
    svg.append('</svg>\n')
    _write(path, '\n'.join(svg))


def _render_family_size_svg(path, rows):
    sizes = [sum(counts) for _, counts in rows]
    counts = Counter(sizes)
    labels = sorted(counts)
    if len(labels) > 40:
        labels = labels[:39] + [max(labels)]
    width, height = 980, 430
    left, right, top, bottom = 62, 28, 82, 72
    chart_w, chart_h = width - left - right, height - top - bottom
    max_count = max([counts[label] for label in labels] + [1])
    bar_w = chart_w / max(1, len(labels))
    svg = _svg_header(width, height)
    svg.append('<text class="title" x="28" y="38">Family Size Distribution</text>')
    svg.append('<text class="subtitle" x="28" y="62">Number of sequences assigned to each final family</text>')
    for index, label in enumerate(labels):
        value = counts[label]
        x = left + index * bar_w
        bar_h = chart_h * value / max_count
        svg.append(
            f'<rect x="{x:.2f}" y="{top + chart_h - bar_h:.2f}" width="{max(1, bar_w - 2):.2f}" '
            f'height="{bar_h:.2f}" fill="#50b87f"><title>size={label}; families={value}</title></rect>'
        )
        if len(labels) <= 20 or index % max(1, len(labels) // 15) == 0:
            svg.append(f'<text class="small" x="{x + bar_w / 2:.2f}" y="{top + chart_h + 18}" text-anchor="middle">{label}</text>')
    svg.append('</svg>\n')
    _write(path, '\n'.join(svg))


def _write_interactive_explorer(path, rows, taxa, ledger):
    data = []
    for (family_id, counts) in rows:
        family = ledger.families[family_id]
        data.append({
            'id': family_id,
            'category': family.category,
            'stage': family.stage,
            'prevalence': family.prevalence,
            'sequences': len(family.members),
            'counts': counts,
        })
    payload = json.dumps({'taxa': taxa, 'families': data}, separators=(',', ':')).replace('</', '<\\/')
    content = '''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>PyamilySeq Pangenome Explorer</title>
<style>
body{font-family:Arial,Helvetica,sans-serif;margin:24px;color:#1f2933;background:#f7f9fb}main{max-width:1400px;margin:auto}
.controls{display:flex;gap:16px;align-items:end;flex-wrap:wrap;background:#fff;border:1px solid #d9e2ec;padding:14px;border-radius:6px}
label{font-size:12px;color:#52606d;display:grid;gap:4px}input,select{font:inherit;padding:6px;border:1px solid #9aa5b1;border-radius:4px}
#viewport{margin-top:16px;background:#fff;border:1px solid #d9e2ec;overflow:auto;max-height:75vh}canvas{display:block}.status{margin:10px 0;color:#52606d}
</style></head><body><main><h1>Pangenome Explorer</h1>
<div class="controls"><label>Family search<input id="search" type="search"></label><label>Stage<select id="stage"><option value="">All</option><option>First</option><option>extended</option><option>combined</option><option>only_Second</option></select></label><label>Minimum prevalence<input id="minimum" type="number" min="0" value="0"></label><label>Rows<input id="limit" type="number" min="1" max="2000" value="300"></label></div>
<p class="status" id="status"></p><div id="viewport"><canvas id="matrix"></canvas></div>
<script>const DATA=__DATA__;const canvas=document.getElementById('matrix'),ctx=canvas.getContext('2d');
const controls=['search','stage','minimum','limit'].map(id=>document.getElementById(id));controls.forEach(c=>c.addEventListener('input',draw));
function draw(){const q=controls[0].value.toLowerCase(),stage=controls[1].value,min=+controls[2].value||0,limit=Math.max(1,+controls[3].value||300);
const rows=DATA.families.filter(f=>(!q||f.id.toLowerCase().includes(q)||f.category.toLowerCase().includes(q))&&(!stage||f.stage===stage)&&f.prevalence>=min).slice(0,limit);
const cw=8,rh=12,left=170,top=120;canvas.width=Math.max(700,left+DATA.taxa.length*cw+20);canvas.height=Math.max(180,top+rows.length*rh+20);ctx.clearRect(0,0,canvas.width,canvas.height);ctx.font='10px Arial';ctx.fillStyle='#52606d';
DATA.taxa.forEach((t,i)=>{const short=t.length>18?t.slice(0,6)+'...'+t.slice(-9):t;ctx.save();ctx.translate(left+i*cw+5,top-8);ctx.rotate(-Math.PI/3);ctx.fillText(short,0,0);ctx.restore()});
const colours=['#edf1f4','#2f7ed8','#50b87f','#f2a93b','#d95763'];rows.forEach((f,r)=>{ctx.fillStyle='#52606d';ctx.textAlign='right';ctx.fillText(f.id,left-7,top+r*rh+9);ctx.textAlign='left';f.counts.forEach((v,c)=>{ctx.fillStyle=colours[Math.min(v,4)];ctx.fillRect(left+c*cw,top+r*rh,cw-1,rh-1)})});
document.getElementById('status').textContent=`Showing ${rows.length} of ${DATA.families.length} families across ${DATA.taxa.length} taxa`;}
draw();</script></main></body></html>'''.replace('__DATA__', payload)
    _write(path, content)


def write_ledger_visualisations(options, ledger, family_ids, categories, mode_label='Species'):
    """Write static summaries plus a scalable canvas-based final-family explorer."""

    output_dir = os.path.join(os.path.abspath(options.output_dir), 'Visualisations')
    os.makedirs(output_dir, exist_ok=True)
    selected = set(family_ids)
    cluster_taxa = OrderedDict(
        (family.family_id, list(family.taxa))
        for family in ledger.sorted_families(selected)
    )
    cluster_sequences = OrderedDict(
        (family.family_id, list(family.members))
        for family in ledger.sorted_families(selected)
    )
    selected_categories = OrderedDict(
        (category, [family_id for family_id in ids if family_id in selected])
        for category, ids in categories.items()
    )
    taxa = list(ledger.genomes if mode_label == 'Species' else ledger.taxa)
    taxa_label = 'genomes' if mode_label == 'Species' else 'genera'
    limit = max(1, int(getattr(options, 'visualise_heatmap_limit', 200) or 200))

    category_counts = OrderedDict((key, len(values)) for key, values in selected_categories.items())
    _render_category_svg(category_counts, os.path.join(output_dir, 'core_categories.svg'), mode_label, taxa_label)
    _render_histogram_svg(cluster_taxa, os.path.join(output_dir, 'family_prevalence_histogram.svg'), mode_label, taxa_label)
    _render_heatmap_svg(
        cluster_taxa, cluster_sequences, taxa,
        os.path.join(output_dir, 'presence_absence_heatmap.svg'), mode_label, taxa_label, limit,
    )
    _write_data_tables(output_dir, category_counts, cluster_taxa, cluster_sequences)

    rows = _copy_number_matrix(ledger, selected, taxa, mode_label)
    _write_copy_number_table(os.path.join(output_dir, 'copy_number_matrix.tsv'), rows, taxa)
    _render_copy_number_svg(rows, taxa, os.path.join(output_dir, 'copy_number_heatmap.svg'), limit)
    accumulation = _accumulation_rows(
        rows,
        taxa,
        permutations=int(getattr(options, 'visualise_permutations', 100) or 100),
        seed=int(getattr(options, 'visualise_seed', 1) or 1),
    )
    _write_accumulation_table(os.path.join(output_dir, 'pangenome_accumulation.tsv'), accumulation)
    _render_accumulation_svg(os.path.join(output_dir, 'pangenome_accumulation.svg'), accumulation)
    _render_family_size_svg(os.path.join(output_dir, 'family_size_distribution.svg'), rows)
    _write_interactive_explorer(os.path.join(output_dir, 'pangenome_explorer.html'), rows, taxa, ledger)
    _render_index(output_dir, mode_label, taxa_label, len(taxa), len(rows), sum(len(ledger.families[f].members) for f in selected), limit)

    index_path = os.path.join(output_dir, 'pangenome_visualisations.html')
    with open(index_path, encoding='utf-8') as handle:
        index = handle.read()
    additions = '''
  <section><h2>Copy Number</h2><img src="copy_number_heatmap.svg" alt="Copy-number heatmap"></section>
  <section><h2>Pangenome Accumulation</h2><img src="pangenome_accumulation.svg" alt="Pangenome accumulation curves"></section>
  <section><h2>Family Sizes</h2><img src="family_size_distribution.svg" alt="Family size distribution"></section>
  <section><h2>Interactive Explorer</h2><p><a href="pangenome_explorer.html">Open the searchable canvas matrix</a></p></section>
'''
    index = index.replace('</main>', additions + '</main>')
    _write(index_path, index)
    return output_dir
