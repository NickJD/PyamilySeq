"""Reports for genome additions to, and queries against, PyamilySeq pangenomes."""

from __future__ import annotations

import argparse
import csv
import html
import os
from collections import Counter, defaultdict


RECONCILIATION_STAGE_LABELS = {
    "First": "Original-only family",
    "extended": "Original family + supplement",
    "combined": "Merged original families",
    "only_Second": "Supplement-only family",
}


def _reconciliation_stage_label(stage):
    return RECONCILIATION_STAGE_LABELS.get(
        str(stage or ""),
        str(stage or "Unknown"),
    )


def _prevalence_display_label(group, groups=None, mode="Species"):
    """Format a final-family prevalence bin without changing its stored key."""
    if mode == "Species":
        values = sorted(
            {int(value) for value in (groups or [99, 95, 15, 0]) if int(value) > 0},
            reverse=True,
        )
        value = int(group)
        if value == 0:
            return f"<{min(values)}% prevalence" if values else "below configured thresholds"
        higher = [item for item in values if item > value]
        if not higher:
            return f">={value}% prevalence"
        return f"{value}-<{min(higher)}% prevalence"

    values = sorted({int(value) for value in (groups or []) if str(value).isdigit()})
    if str(group) == ">":
        return f">{values[-1]} genera" if values else "above configured genus bins"
    value = int(group)
    lower_values = [item for item in values if item < value]
    lower = (max(lower_values) + 1) if lower_values else 1
    return f"<={value} genera" if lower == 1 else f"{lower}-{value} genera"


def _display_family_category(category, groups=None, mode="Species"):
    """Readable display label for an internal stage/prevalence category."""
    category = str(category or "")
    suffix_name = "core" if mode == "Species" else "genera"

    for stage, stage_label in RECONCILIATION_STAGE_LABELS.items():
        prefix = f"{stage}_{suffix_name}_"
        if category.startswith(prefix):
            raw_group = category[len(prefix):]
            prevalence_label = _prevalence_display_label(
                raw_group, groups=groups, mode=mode
            )
            return f"{stage_label} — {prevalence_label}"

    return category.replace("_", " ")


def write_extension_report(output_dir, ledger, genome_ids, groups=None, mode=None):
    """Describe how added genomes map onto the reconciled original family catalogue."""
    genome_ids = str(genome_ids)
    if os.path.isfile(genome_ids):
        with open(genome_ids, encoding="utf-8-sig") as handle:
            values = [line.strip() for line in handle if line.strip() and not line.lstrip().startswith("#")]
    else:
        values = genome_ids.split(",")
    wanted = {item.strip() for item in values if item.strip()}
    known = set(ledger.genomes)
    missing = sorted(wanted - known)
    if missing:
        raise ValueError("Added genome IDs not found in sequence IDs: " + ", ".join(missing))
    if not wanted:
        raise ValueError("--added-genomes requires at least one genome ID")

    mode = mode or getattr(ledger.resolver, "mode", "Species")
    if groups is None:
        groups = [99, 95, 15, 0] if mode == "Species" else [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, ">"]

    os.makedirs(output_dir, exist_ok=True)
    summary = defaultdict(Counter)
    families = {}
    details_path = os.path.join(output_dir, "added_genome_family_links.tsv")
    fields = [
        "genome",
        "gene_id",
        "final_family_id",
        "family_status_code",
        "family_status",
        "family_category",                 # Internal category; retained for scripts.
        "family_category_label",           # Human-readable category.
        "family_genome_count",
        "family_gene_count",
        "original_cluster_count",
        "original_cluster_ids",
        "source_lineage_ids",
        "family_stage",                    # Internal stage; retained for scripts.
        "reconciliation_stage_label",      # Human-readable provenance.
    ]
    with open(details_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sequence_id in ledger.sequences:
            genome = ledger.resolver.genome(sequence_id)
            if genome not in wanted:
                continue
            family = ledger.families[ledger.assignment[sequence_id]]
            # An original source family exists when the final reconciled family
            # descends from a first-round cluster. New-only families have none.
            original_clusters = [label for label in family.source_clusters if not label.startswith("singleton:")]
            status_code = "linked" if family.stage != "only_Second" else "supplement_only"
            status_label = (
                "Linked to original pangenome"
                if status_code == "linked"
                else "Supplement-only family"
            )
            category_label = _display_family_category(
                family.category, groups=groups, mode=mode
            )
            stage_label = _reconciliation_stage_label(family.stage)

            summary[genome][status_code] += 1
            summary[genome]["genes"] += 1
            # Category counts are counts of ADDED GENES, not counts of families.
            summary[genome][family.category] += 1

            writer.writerow({
                "genome": genome,
                "gene_id": sequence_id,
                "final_family_id": family.family_id,
                "family_status_code": status_code,
                "family_status": status_label,
                "family_category": family.category,
                "family_category_label": category_label,
                "family_genome_count": family.prevalence,
                "family_gene_count": len(family.members),
                "original_cluster_count": len(original_clusters),
                "original_cluster_ids": ";".join(original_clusters),
                "source_lineage_ids": ";".join(family.source_families),
                "family_stage": family.stage,
                "reconciliation_stage_label": stage_label,
            })

            family_row = families.setdefault(family.family_id, {
                "genome": genome,
                "family_id": family.family_id,
                "status_code": status_code,
                "status": status_label,
                "category": family.category,
                "category_label": category_label,
                "stage": family.stage,
                "stage_label": stage_label,
                "family_genomes": family.prevalence,
                "family_genes": len(family.members),
                "original_clusters": original_clusters,
                "added_genes": [],
                "source_ids": family.source_families,
            })
            family_row["added_genes"].append(sequence_id)
    summary_path = os.path.join(output_dir, "added_genome_summary.tsv")
    with open(summary_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        # Keep the previous summary columns (including the legacy
        # genes_in_novel_families alias) so existing TSV consumers continue
        # to work. New code should use the explicit supplement-only column.
        writer.writerow([
            "genome",
            "added_genes",  # Legacy alias for added_gene_count.
            "genes_linked_to_original_pangenome",
            "genes_in_novel_families",  # Deprecated alias; not evolutionary novelty.
            "family_categories",  # Legacy internal-category count string.
            "added_gene_count",
            "genes_in_supplement_only_families",
            "added_gene_counts_by_internal_family_category",
            "added_gene_counts_by_display_family_category",
        ])
        for genome in sorted(wanted):
            counts = summary[genome]
            category_items = [
                (key, counts[key])
                for key in sorted(counts)
                if key not in {"genes", "linked", "supplement_only"}
            ]
            internal_categories = ";".join(
                f"{key}:{count}" for key, count in category_items
            )
            display_categories = ";".join(
                f"{_display_family_category(key, groups=groups, mode=mode)}:{count}"
                for key, count in category_items
            )
            writer.writerow([
                genome,
                counts["genes"],          # Legacy added_genes column.
                counts["linked"],
                counts["supplement_only"],  # Legacy genes_in_novel_families alias.
                internal_categories,      # Legacy family_categories column.
                counts["genes"],
                counts["supplement_only"],
                internal_categories,
                display_categories,
            ])

    family_summary_path = os.path.join(output_dir, "added_genome_family_summary.tsv")
    family_fields = [
        "genome",
        "final_family_id",
        "family_status_code",
        "family_status",
        "family_category",
        "family_category_label",
        "family_stage",
        "reconciliation_stage_label",
        "genomes_in_family",
        "genes_in_family",
        "added_gene_count",
        "added_gene_ids",
        "original_cluster_count",
        "original_cluster_ids",
        "source_lineage_ids",
    ]
    with open(family_summary_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=family_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for family_id in sorted(families):
            row = families[family_id]
            writer.writerow({
                "genome": row["genome"],
                "final_family_id": family_id,
                "family_status_code": row["status_code"],
                "family_status": row["status"],
                "family_category": row["category"],
                "family_category_label": row["category_label"],
                "family_stage": row["stage"],
                "reconciliation_stage_label": row["stage_label"],
                "genomes_in_family": row["family_genomes"],
                "genes_in_family": row["family_genes"],
                "added_gene_count": len(row["added_genes"]),
                "added_gene_ids": ";".join(row["added_genes"]),
                "original_cluster_count": len(row["original_clusters"]),
                "original_cluster_ids": ";".join(row["original_clusters"]),
                "source_lineage_ids": ";".join(row["source_ids"]),
            })
    _write_extension_svg(
        os.path.join(output_dir, "added_genome_summary.svg"),
        summary,
        sorted(wanted),
    )
    _write_extension_html(
        output_dir,
        summary,
        sorted(wanted),
        details_path,
        summary_path,
        family_summary_path,
        groups=groups,
        mode=mode,
    )
    return summary_path, details_path


def _write_extension_svg(path, summary, genomes):
    width, height = 1100, max(180, 100 + 68 * len(genomes))
    colors = [
        ("linked", "#247a62", "Genes in original-linked families"),
        ("supplement_only", "#c66a32", "Genes in supplement-only families"),
    ]
    max_total = max((summary[g]["genes"] for g in genomes), default=1) or 1
    left, bar_width = 345, 650
    chunks = [f"<svg xmlns='http://www.w3.org/2000/svg' role='img' aria-labelledby='title desc' width='100%' height='{height}' viewBox='0 0 {width} {height}'>",
              "<title id='title'>Added-genome gene-family relationships</title><desc id='desc'>Gene counts assigned to families linked to the original pangenome or to supplement-only families.</desc>",
              "<style>text{font:15px system-ui,sans-serif;fill:#182321}.small{font-size:12px;fill:#51615c}</style>",
              "<text x='24' y='30' style='font:bold 21px system-ui,sans-serif'>Where the added genes land</text>"]
    for i, (key, color, label) in enumerate(colors):
        x = 24 + i * 360
        chunks.append(f"<rect x='{x}' y='48' width='14' height='14' rx='3' fill='{color}'/><text x='{x+22}' y='61'>{html.escape(label)}</text>")
    for row, genome in enumerate(genomes):
        y = 94 + row * 68
        label = genome if len(genome) <= 42 else genome[:39] + "…"
        chunks.append(f"<text x='24' y='{y+18}'>{html.escape(label)}</text>")
        if label != genome:
            chunks.append(f"<title>{html.escape(genome)}</title>")
        x, scale = left, bar_width / max_total
        for key, color, _ in colors:
            amount = summary[genome][key]
            bar = amount * scale
            if bar:
                chunks.append(f"<rect x='{x:.1f}' y='{y}' width='{bar:.1f}' height='30' rx='3' fill='{color}'/>")
                if bar > 48:
                    chunks.append(f"<text x='{x+bar/2:.1f}' y='{y+20}' text-anchor='middle' style='font-size:12px;fill:white'>{amount:,}</text>")
                x += bar
        chunks.append(f"<text x='{left}' y='{y+49}' class='small'>Total genes: {summary[genome]['genes']:,} · linked: {summary[genome]['linked']:,} · supplement-only: {summary[genome]['supplement_only']:,}</text>")
    chunks.append("</svg>")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(chunks))


def _write_extension_html(
    output_dir,
    summary,
    genomes,
    details_path,
    summary_path,
    family_summary_path,
    groups=None,
    mode="Species",
):
    total = sum(summary[genome]["genes"] for genome in genomes)
    linked = sum(summary[genome]["linked"] for genome in genomes)
    supplement_only = sum(summary[genome]["supplement_only"] for genome in genomes)
    pct = (100 * linked / total) if total else 0
    genome_options = "".join(
        f"<option value='{html.escape(genome, quote=True)}'>{html.escape(genome)}</option>"
        for genome in genomes
    )
    family_rows = _family_html_rows(
        os.path.join(output_dir, "added_genome_family_summary.tsv")
    )

    per_genome_rows = []
    for genome in genomes:
        counts = summary[genome]
        category_items = [
            (key, counts[key])
            for key in sorted(counts)
            if key not in {"genes", "linked", "supplement_only"}
        ]
        category_text = "; ".join(
            f"{_display_family_category(key, groups=groups, mode=mode)}: {count:,}"
            for key, count in category_items
        )
        per_genome_rows.append(
            "<tr>"
            f"<td>{html.escape(genome)}</td>"
            f"<td>{counts['genes']:,}</td>"
            f"<td>{counts['linked']:,}</td>"
            f"<td>{counts['supplement_only']:,}</td>"
            f"<td>{html.escape(category_text)}</td>"
            "</tr>"
        )

    body = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Added-genome impact</title><style>
:root{{--ink:#182321;--muted:#566660;--line:#dce5e1;--paper:#f5f8f6;--accent:#247a62;--new:#c66a32}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:15px/1.55 system-ui,-apple-system,sans-serif}}
main{{max-width:1440px;margin:0 auto;padding:32px 28px 64px}}h1{{font-size:clamp(28px,4vw,42px);line-height:1.12;margin:0 0 8px}}h2{{font-size:22px;margin:34px 0 8px}}p{{color:var(--muted);margin:5px 0 18px}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:12px;margin:22px 0}}.card,section{{background:white;border:1px solid var(--line);border-radius:12px;padding:18px}}
.card b{{display:block;font-size:28px;line-height:1.3}}.card span{{color:var(--muted)}}.chart{{padding:10px 16px;overflow:hidden}}.chart img{{display:block;width:100%;height:auto}}
.toolbar{{display:flex;flex-wrap:wrap;gap:10px;margin:14px 0}}input,select{{border:1px solid #bdcbc5;border-radius:8px;padding:9px 11px;background:white;font:inherit}}input{{min-width:min(340px,100%)}}
.table-wrap{{overflow:auto;max-height:660px;border:1px solid var(--line);border-radius:8px}}table{{border-collapse:collapse;width:100%;font-size:13px;background:white}}th,td{{text-align:left;vertical-align:top;padding:9px 10px;border-bottom:1px solid #e8eeeb}}th{{position:sticky;top:0;background:#edf3f0;z-index:1;white-space:nowrap}}td{{max-width:360px;overflow-wrap:anywhere}}code{{font-size:12px}}.badge{{display:inline-block;padding:2px 8px;border-radius:999px;background:#e5f2ed;color:#165b47;white-space:nowrap}}.badge.supplement-only{{background:#fbede4;color:#88451e}}
.glossary{{display:grid;grid-template-columns:minmax(180px,260px) 1fr;gap:7px 16px;margin:12px 0}}.glossary dt{{font-weight:650}}.glossary dd{{margin:0;color:var(--muted)}}.links a{{color:#17674f;margin-right:16px}}.note{{font-size:13px}}
@media(max-width:640px){{main{{padding:22px 14px}}.glossary{{grid-template-columns:1fr;gap:2px}}.glossary dd{{margin:0 0 10px}}}}
</style></head><body><main>
<h1>Added-genome impact</h1>
<p>Gene-level summary of how genes from the added genome(s) are assigned to final families in the supplemented pangenome.</p>
<div class="cards">
  <div class="card"><b>{total:,}</b><span>added genes</span></div>
  <div class="card"><b>{linked:,}</b><span>genes assigned to families linked to the original pangenome ({pct:.1f}%)</span></div>
  <div class="card"><b>{supplement_only:,}</b><span>genes assigned to supplement-only families</span></div>
  <div class="card"><b>{len(genomes):,}</b><span>added genomes summarised</span></div>
</div>
<section class="chart"><img src="added_genome_summary.svg" alt="Gene counts assigned to families linked to the original pangenome or to supplement-only families"></section>
<h2>Per-genome summary</h2>
<p>All values in this section are gene counts, not family counts. Each added gene is counted once, according to the final family to which it was assigned.</p>
<section class="table-wrap"><table><thead><tr>
<th>Genome</th><th>Added genes</th><th>Genes in families linked to original pangenome</th>
<th>Genes in supplement-only families</th><th>Added genes by final family category (gene counts)</th>
</tr></thead><tbody>
{''.join(per_genome_rows)}
</tbody></table></section>
<h2>Family-level results</h2>
<p>One row per final family containing at least one added gene. The added gene count is the number of genes from the selected added genome(s) assigned to that family.</p>
<div class="toolbar">
<input id="search" type="search" placeholder="Search family, gene, or original cluster ID" aria-label="Search family results">
<select id="status">
<option value="">All family relationships</option>
<option value="linked">Linked to original pangenome</option>
<option value="supplement_only">Supplement-only family</option>
</select>
<select id="genome"><option value="">All genomes</option>{genome_options}</select>
</div>
<div class="table-wrap"><table id="families"><thead><tr>
<th>Family relationship</th><th>Genome</th><th>Final family ID</th><th>Reconciliation provenance and prevalence</th>
<th>Genomes in family</th><th>Genes in family</th><th>Added gene count</th><th>Added gene IDs</th><th>Original cluster IDs</th>
</tr></thead><tbody>
{family_rows}
</tbody></table></div><p id="rowCount" class="note"></p>
<h2>How to read this report</h2><dl class="glossary">
<dt>Linked to original pangenome</dt><dd>The final family contains at least one family lineage from the original clustering. The added gene is therefore assigned to a family with an identifiable original-family contribution.</dd>
<dt>Supplement-only family</dt><dd>The final family has no identifiable original-family lineage. This describes its relationship to the supplied baseline pangenome; it does not establish that the sequence is evolutionarily novel.</dd>
<dt>Reconciliation provenance and prevalence</dt><dd>The first part describes how the final family was formed during reconciliation; the second part describes its prevalence across the current genome set. These are separate properties.</dd>
<dt>Original cluster IDs</dt><dd>Labels from the original clustering input that contributed to this family. These are useful for tracing a final family back to the original clusters.</dd>
<dt>Source lineage IDs</dt><dd>Internal stable identifiers used to preserve exact provenance. They are identifiers, not biological annotations.</dd>
<dt>family_summary.tsv</dt><dd>One row per final family across the pangenome, including prevalence and copy-number statistics. Source lineage columns record provenance rather than family names.</dd>
</dl>
<p class="links">
<a href="added_genome_summary.tsv">Download per-genome summary TSV</a>
<a href="added_genome_family_summary.tsv">Download family summary TSV</a>
<a href="added_genome_family_links.tsv">Download gene-level links TSV</a>
<a href="family_summary.tsv">Download full pangenome family summary</a>
</p>
<script>
const search=document.getElementById('search'),status=document.getElementById('status'),genome=document.getElementById('genome'),
rows=[...document.querySelectorAll('#families tbody tr')],count=document.getElementById('rowCount');
function filterRows(){{
  const q=search.value.toLowerCase();let shown=0;
  for(const row of rows){{
    const ok=(!q||row.dataset.search.includes(q))
      &&(!status.value||row.dataset.status===status.value)
      &&(!genome.value||row.dataset.genome===genome.value);
    row.hidden=!ok;if(ok)shown++;
  }}
  count.textContent=`Showing ${{shown.toLocaleString()}} of ${{rows.length.toLocaleString()}} final families.`;
}}
search.addEventListener('input',filterRows);status.addEventListener('change',filterRows);genome.addEventListener('change',filterRows);filterRows();
</script></main></body></html>"""
    with open(os.path.join(output_dir, "added_genome_report.html"), "w", encoding="utf-8") as handle:
        handle.write(body)


def _family_html_rows(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            status = row["family_status"]
            status_code = row.get("family_status_code", "")
            if not status_code:
                status_code = (
                    "supplement_only"
                    if status in {"Novel family", "Supplement-only family"}
                    else "linked"
                )

            category_label = row.get("family_category_label") or row.get("family_category", "")
            provenance_label = row.get("reconciliation_stage_label", "")
            if not provenance_label:
                provenance_label = _reconciliation_stage_label(row.get("family_stage", ""))

            values = [
                status,
                row["genome"],
                row["final_family_id"],
                category_label,
                row["genomes_in_family"],
                row["genes_in_family"],
                row["added_gene_count"],
                row["added_gene_ids"],
                row["original_cluster_ids"] or "—",
            ]
            cells = []
            for index, value in enumerate(values):
                escaped = html.escape(str(value))
                if index == 0:
                    badge_class = "badge supplement-only" if status_code == "supplement_only" else "badge"
                    cells.append(f'<td><span class="{badge_class}">{escaped}</span></td>')
                elif index == 3:
                    cells.append(
                        f'<td>{escaped}<div class="note">{html.escape(provenance_label)}</div></td>'
                    )
                else:
                    cells.append(f"<td>{escaped}</td>")

            search = " ".join([
                row["final_family_id"],
                row["added_gene_ids"],
                row["original_cluster_ids"],
                row["family_category"],
                category_label,
                provenance_label,
            ]).lower()

            rows.append(
                f'<tr data-status="{html.escape(status_code, quote=True)}" '
                f'data-genome="{html.escape(row["genome"], quote=True)}" '
                f'data-search="{html.escape(search, quote=True)}">'
                + "".join(cells)
                + "</tr>"
            )
    return "\n".join(rows)


def query_pangenome(gpa_path, genome_list_path, output_dir, core_percent=99.0):
    """Summarize selected genomes from a Roary/PyamilySeq GPA CSV."""
    if not 0 <= core_percent <= 100:
        raise ValueError("core_percent must be between 0 and 100")
    with open(genome_list_path, encoding="utf-8") as handle:
        requested = [line.strip().split(",", 1)[0].strip() for line in handle if line.strip() and not line.lstrip().startswith("#")]
    requested = list(dict.fromkeys(requested))
    if not requested:
        raise ValueError("Genome list is empty")
    with open(gpa_path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "Gene" not in reader.fieldnames:
            raise ValueError("Input must be a gene_presence_absence.csv with a Gene column")
        genome_cols = [c for c in reader.fieldnames if c not in {"Gene", "Non-unique Gene name", "Annotation", "No. isolates", "No. sequences", "Avg sequences per isolate", "Genome Fragment", "Order within Fragment", "Accessory Fragment", "Accessory Order with Fragment", "QC", "Min group size nuc", "Max group size nuc", "Avg group size nuc"}]
        missing = sorted(set(requested) - set(genome_cols))
        if missing:
            raise ValueError("Requested genomes absent from GPA columns: " + ", ".join(missing))
        rows = list(reader)
    if not genome_cols:
        raise ValueError("No genome columns found in GPA file")
    os.makedirs(output_dir, exist_ok=True)
    detail_fields = ["genome", "family_id", "category", "family_prevalence_percent", "gene_copies", "genes"]
    summaries = []
    with open(os.path.join(output_dir, "genome_query_families.tsv"), "w", newline="", encoding="utf-8") as detail:
        writer = csv.DictWriter(detail, fieldnames=detail_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for genome in requested:
            category_counts = Counter(); gene_counts = Counter(); total = 0
            for row in rows:
                members = [cell.strip() for cell in (row.get(genome) or "").split() if cell.strip()]
                if not members:
                    continue
                total += len(members)
                prevalence = sum(bool((row.get(col) or "").strip()) for col in genome_cols) / len(genome_cols) * 100
                category = "core" if prevalence >= core_percent else "accessory"
                category_counts[category] += 1; gene_counts[category] += len(members)
                writer.writerow({"genome": genome, "family_id": row["Gene"], "category": category,
                                 "family_prevalence_percent": f"{prevalence:.2f}", "gene_copies": len(members), "genes": ";".join(members)})
            summaries.append([genome, total, category_counts["core"], gene_counts["core"], category_counts["accessory"], gene_counts["accessory"], len(genome_cols)])
    with open(os.path.join(output_dir, "genome_query_summary.tsv"), "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome", "genes_in_families", "core_families", "genes_in_core_families", "accessory_families", "genes_in_accessory_families", "pangenome_genomes"])
        writer.writerows(summaries)
    return summaries


def main(argv=None):
    parser = argparse.ArgumentParser(description="Query genomes against a PyamilySeq gene_presence_absence.csv")
    parser.add_argument("--gpa", required=True, help="PyamilySeq/Roary gene_presence_absence.csv")
    parser.add_argument("--genomes", required=True, help="Text file with one genome ID per line")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--core-percent", type=float, default=99.0, help="Core prevalence threshold (default: 99)")
    args = parser.parse_args(argv)
    query_pangenome(args.gpa, args.genomes, args.output_dir, args.core_percent)


if __name__ == "__main__":
    main()