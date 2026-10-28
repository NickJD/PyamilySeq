"""Reports for genome additions to, and queries against, PyamilySeq pangenomes."""

from __future__ import annotations

import argparse
import csv
import html
import os
from collections import Counter, defaultdict


def write_extension_report(output_dir, ledger, genome_ids):
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
    os.makedirs(output_dir, exist_ok=True)
    summary = defaultdict(Counter)
    families = {}
    details_path = os.path.join(output_dir, "added_genome_family_links.tsv")
    fields = ["genome", "gene_id", "final_family_id", "family_status", "family_category",
              "family_genome_count", "family_gene_count", "original_cluster_count", "original_cluster_ids",
              "source_lineage_ids", "family_stage"]
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
            status_key = "linked" if family.stage != "only_Second" else "novel"
            status_label = "Linked to original pangenome" if status_key == "linked" else "Novel family"
            summary[genome][status_key] += 1
            summary[genome]["genes"] += 1
            summary[genome][family.category] += 1
            writer.writerow({"genome": genome, "gene_id": sequence_id, "final_family_id": family.family_id,
                             "family_status": status_label, "family_category": family.category,
                             "family_genome_count": family.prevalence, "family_gene_count": len(family.members),
                             "original_cluster_count": len(original_clusters),
                             "original_cluster_ids": ";".join(original_clusters),
                             "source_lineage_ids": ";".join(family.source_families), "family_stage": family.stage})
            family_row = families.setdefault(family.family_id, {
                "genome": genome, "family_id": family.family_id, "status": status_label,
                "category": family.category, "stage": family.stage, "family_genomes": family.prevalence,
                "family_genes": len(family.members), "original_clusters": original_clusters,
                "added_genes": [], "source_ids": family.source_families,
            })
            family_row["added_genes"].append(sequence_id)
    summary_path = os.path.join(output_dir, "added_genome_summary.tsv")
    with open(summary_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["genome", "added_genes", "genes_linked_to_original_pangenome", "genes_in_novel_families", "family_categories"])
        for genome in sorted(wanted):
            counts = summary[genome]
            writer.writerow([genome, counts["genes"], counts["linked"], counts["novel"],
                             ";".join(f"{key}:{counts[key]}" for key in sorted(counts) if key not in {"genes", "linked", "novel"})])
    family_summary_path = os.path.join(output_dir, "added_genome_family_summary.tsv")
    family_fields = ["genome", "final_family_id", "family_status", "family_category", "family_stage",
                     "genomes_in_family", "genes_in_family", "added_gene_count", "added_gene_ids",
                     "original_cluster_count", "original_cluster_ids", "source_lineage_ids"]
    with open(family_summary_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=family_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for family_id in sorted(families):
            row = families[family_id]
            writer.writerow({"genome": row["genome"], "final_family_id": family_id, "family_status": row["status"],
                             "family_category": row["category"], "family_stage": row["stage"],
                             "genomes_in_family": row["family_genomes"], "genes_in_family": row["family_genes"],
                             "added_gene_count": len(row["added_genes"]), "added_gene_ids": ";".join(row["added_genes"]),
                             "original_cluster_count": len(row["original_clusters"]),
                             "original_cluster_ids": ";".join(row["original_clusters"]),
                             "source_lineage_ids": ";".join(row["source_ids"])})
    _write_extension_svg(os.path.join(output_dir, "added_genome_summary.svg"), summary, sorted(wanted))
    _write_extension_html(output_dir, summary, sorted(wanted), details_path, summary_path, family_summary_path)
    return summary_path, details_path


def _write_extension_svg(path, summary, genomes):
    width, height = 1100, max(180, 100 + 68 * len(genomes))
    colors = [("linked", "#247a62", "Genes linked to original pangenome"),
              ("novel", "#c66a32", "Genes in novel families")]
    max_total = max((summary[g]["genes"] for g in genomes), default=1) or 1
    left, bar_width = 345, 650
    chunks = [f"<svg xmlns='http://www.w3.org/2000/svg' role='img' aria-labelledby='title desc' width='100%' height='{height}' viewBox='0 0 {width} {height}'>",
              "<title id='title'>Added genome family links</title><desc id='desc'>Counts of genes linked to original families and genes in novel families.</desc>",
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
        chunks.append(f"<text x='{left}' y='{y+49}' class='small'>Total genes: {summary[genome]['genes']:,} · linked: {summary[genome]['linked']:,} · novel: {summary[genome]['novel']:,}</text>")
    chunks.append("</svg>")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(chunks))


def _write_extension_html(output_dir, summary, genomes, details_path, summary_path, family_summary_path):
    total = sum(summary[genome]["genes"] for genome in genomes)
    linked = sum(summary[genome]["linked"] for genome in genomes)
    novel = sum(summary[genome]["novel"] for genome in genomes)
    pct = (100 * linked / total) if total else 0
    genome_options = "".join(f"<option value='{html.escape(genome, quote=True)}'>{html.escape(genome)}</option>" for genome in genomes)
    # Link to PyamilySeq's standard family summary alongside the focused tables.
    body = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Added genome impact</title><style>
:root{{--ink:#182321;--muted:#566660;--line:#dce5e1;--paper:#f5f8f6;--accent:#247a62;--new:#c66a32}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:15px/1.55 system-ui,-apple-system,sans-serif}}
main{{max-width:1440px;margin:0 auto;padding:32px 28px 64px}}h1{{font-size:clamp(28px,4vw,42px);line-height:1.12;margin:0 0 8px}}h2{{font-size:22px;margin:34px 0 8px}}p{{color:var(--muted);margin:5px 0 18px}}
.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:22px 0}}.card,section{{background:white;border:1px solid var(--line);border-radius:12px;padding:18px}}
.card b{{display:block;font-size:28px;line-height:1.3}}.card span{{color:var(--muted)}}.chart{{padding:10px 16px;overflow:hidden}}.chart img{{display:block;width:100%;height:auto}}
.toolbar{{display:flex;flex-wrap:wrap;gap:10px;margin:14px 0}}input,select{{border:1px solid #bdcbc5;border-radius:8px;padding:9px 11px;background:white;font:inherit}}input{{min-width:min(340px,100%)}}
.table-wrap{{overflow:auto;max-height:660px;border:1px solid var(--line);border-radius:8px}}table{{border-collapse:collapse;width:100%;font-size:13px;background:white}}th,td{{text-align:left;vertical-align:top;padding:9px 10px;border-bottom:1px solid #e8eeeb}}th{{position:sticky;top:0;background:#edf3f0;z-index:1;white-space:nowrap}}td{{max-width:360px;overflow-wrap:anywhere}}code{{font-size:12px}}.badge{{display:inline-block;padding:2px 8px;border-radius:999px;background:#e5f2ed;color:#165b47;white-space:nowrap}}.badge.novel{{background:#fbede4;color:#88451e}}
.glossary{{display:grid;grid-template-columns:minmax(180px,260px) 1fr;gap:7px 16px;margin:12px 0}}.glossary dt{{font-weight:650}}.glossary dd{{margin:0;color:var(--muted)}}.links a{{color:#17674f;margin-right:16px}}.note{{font-size:13px}}
@media(max-width:640px){{main{{padding:22px 14px}}.glossary{{grid-template-columns:1fr;gap:2px}}.glossary dd{{margin:0 0 10px}}}}
</style></head><body><main>
<h1>Added genome impact</h1><p>Summary of how genes from the selected genome(s) map to families in the reclustered pangenome.</p>
<div class="cards"><div class="card"><b>{total:,}</b><span>genes in added genome(s)</span></div><div class="card"><b>{linked:,}</b><span>genes linked to an original family ({pct:.1f}%)</span></div><div class="card"><b>{novel:,}</b><span>genes in novel families</span></div><div class="card"><b>{sum(summary[g]['linked'] + summary[g]['novel'] > 0 for g in genomes):,}</b><span>genomes summarized</span></div></div>
<section class="chart"><img src="added_genome_summary.svg" alt="Gene counts linked to existing families or in novel families"></section>
<h2>Per-genome summary</h2><p>Counts are genes, not families. A gene is linked when its final family contains at least one first-round family.</p>
<section class="table-wrap"><table><thead><tr><th>Genome</th><th>Added genes</th><th>Linked to original pangenome</th><th>In novel families</th><th>Family category counts</th></tr></thead><tbody>
{''.join(f'<tr><td>{html.escape(g)}</td><td>{summary[g]["genes"]:,}</td><td>{summary[g]["linked"]:,}</td><td>{summary[g]["novel"]:,}</td><td>{html.escape("; ".join(f"{k}: {summary[g][k]}" for k in sorted(summary[g]) if k not in {"genes", "linked", "novel"}))}</td></tr>' for g in genomes)}
</tbody></table></section>
<h2>Family-level results</h2><p>One row per final family that contains at least one added gene. Use the search box to find a gene or family ID.</p>
<div class="toolbar"><input id="search" type="search" placeholder="Search family, gene, or cluster ID" aria-label="Search family results"><select id="status"><option value="">All family statuses</option><option>Linked to original pangenome</option><option>Novel family</option></select><select id="genome"><option value="">All genomes</option>{genome_options}</select></div>
<div class="table-wrap"><table id="families"><thead><tr><th>Family status</th><th>Genome</th><th>Final family ID</th><th>Category</th><th>Genomes in family</th><th>Genes in family</th><th>Added gene count</th><th>Added gene IDs</th><th>Original cluster IDs</th></tr></thead><tbody>
{_family_html_rows(os.path.join(output_dir, "added_genome_family_summary.tsv"))}
</tbody></table></div><p id="rowCount" class="note"></p>
<h2>How to read the files</h2><dl class="glossary">
<dt>Final family ID</dt><dd>Stable PyamilySeq identifier for the final reconciled family.</dd>
<dt>Family status</dt><dd>“Linked to original pangenome” means the final family includes a first-round family. “Novel family” means it has no first-round family member.</dd>
<dt>Original cluster IDs</dt><dd>Labels from the original clustering input that contributed to this family. These are the most useful field for tracing a family back to the original clusters.</dd>
<dt>Family category</dt><dd>Prevalence bin and reconciliation stage, such as <code>extended_core_99</code>. The stage is First, extended, combined, or only_Second; the suffix is its prevalence group.</dd>
<dt>Source lineage IDs</dt><dd>Internal stable IDs used to preserve exact provenance. They are identifiers, not biological annotations; use original cluster IDs for interpretation.</dd>
<dt>family_summary.tsv</dt><dd>One row per retained final family across the whole pangenome. It includes copy number and prevalence statistics. The lineage columns are provenance fields, not family names.</dd>
</dl>
<p class="links"><a href="added_genome_summary.tsv">Download per-genome summary TSV</a><a href="added_genome_family_summary.tsv">Download family summary TSV</a><a href="added_genome_family_links.tsv">Download gene-level links TSV</a><a href="family_summary.tsv">Download full pangenome family summary</a></p>
<script>
const search=document.getElementById('search'),status=document.getElementById('status'),genome=document.getElementById('genome'),rows=[...document.querySelectorAll('#families tbody tr')],count=document.getElementById('rowCount');
function filterRows(){{const q=search.value.toLowerCase();let shown=0;for(const row of rows){{const ok=(!q||row.dataset.search.includes(q))&&(!status.value||row.dataset.status===status.value)&&(!genome.value||row.dataset.genome===genome.value);row.hidden=!ok;if(ok)shown++}}count.textContent=`Showing ${{shown.toLocaleString()}} of ${{rows.length.toLocaleString()}} families.`}}
search.addEventListener('input',filterRows);status.addEventListener('change',filterRows);genome.addEventListener('change',filterRows);filterRows();
</script></main></body></html>'''
    with open(os.path.join(output_dir, "added_genome_report.html"), "w", encoding="utf-8") as handle:
        handle.write(body)


def _family_html_rows(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            status = row["family_status"]
            values = [status, row["genome"], row["final_family_id"], row["family_category"],
                      row["genomes_in_family"], row["genes_in_family"], row["added_gene_count"],
                      row["added_gene_ids"], row["original_cluster_ids"] or "—"]
            cells = "".join(f"<td>{'<span class=\"badge' + (' novel' if status == 'Novel family' else '') + '\">' + html.escape(value) + '</span>' if i == 0 else html.escape(value)}</td>" for i, value in enumerate(values))
            search = " ".join([row["final_family_id"], row["added_gene_ids"], row["original_cluster_ids"], row["family_category"]]).lower()
            rows.append(f'<tr data-status="{html.escape(status, quote=True)}" data-genome="{html.escape(row["genome"], quote=True)}" data-search="{html.escape(search, quote=True)}">{cells}</tr>')
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
