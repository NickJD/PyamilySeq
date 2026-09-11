#!/usr/bin/env python3
"""Write coordinate-ordered family paths for each genome.

The gene-family table is authoritative for identifiers and the GFF is
authoritative for contig/coordinate order.  FASTA order is never used as a
surrogate for genomic order.
"""

import argparse
import csv
import gzip
import os
from collections import defaultdict


RANK_INDEX = {
    "species": -1,
    "genus": -2,
    "family": -3,
    "order": -4,
    "class": -5,
    "phylum": -6,
    "domain": 0,
}


def open_text(path):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    return opener(path, 'rt', encoding='utf-8')


def load_lineages(path):
    result = {}
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            result[fields[0]] = [token for token in fields[-1].split(';') if token]
    return result


def load_genome_groups(path, lineages, rank):
    result = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'Sample', 'Species'}
        if not required <= set(reader.fieldnames or ()):
            raise ValueError(f'{path} must contain Sample and Species columns')
        for row in reader:
            lineage = lineages.get(row['Species'])
            if not lineage:
                continue
            index = RANK_INDEX[rank]
            result[row['Sample']] = lineage[index] if abs(index) <= len(lineage) else lineage[-1]
    return result

def read_gene_families(path, delimiter='|'):
    gene_to_family = defaultdict(dict)
    gene_counts = defaultdict(int)

    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter='\t')

        for line_number, row in enumerate(reader, start=1):
            if not row or row[0].startswith('#'):
                continue
            if len(row) < 2:
                raise ValueError(f"Gene map line {line_number} has fewer than two columns")
            if line_number == 1 and row[0].lower() in {'sequence_id', 'gene_id'}:
                continue

            gene, family = row[:2]
            genome = gene.split(delimiter, 1)[0]
            genome_map = gene_to_family[genome]

            previous = genome_map.get(gene)
            if previous is not None and previous != family:
                raise ValueError(f"Gene {gene!r} is assigned to multiple families")

            genome_map[gene] = family
            gene_counts[genome] += 1
            # if gene in gene_to_family and gene_to_family[gene] != family:
            #     raise ValueError(
            #         f'Gene {gene!r} is assigned to multiple families'
            #     )
            # 
            # gene_to_family[gene] = family
            # 
            # genome = gene.split(delimiter, 1)[0]
            # gene_counts[genome] += 1

    return gene_to_family, gene_counts
# def read_gene_families(path, delimiter='|'):
#     gene_to_family = {}
#     with open_text(path) as handle:
#         reader = csv.reader(handle, delimiter='\t')
#         for line_number, row in enumerate(reader, start=1):
#             if not row or row[0].startswith('#'):
#                 continue
#             if line_number == 1 and row[0].lower() in {'sequence_id', 'gene_id'}:
#                 continue
#             if len(row) < 2:
#                 raise ValueError(f'Gene map line {line_number} has fewer than two columns')
#             gene, family = row[:2]
#             if gene in gene_to_family and gene_to_family[gene] != family:
#                 raise ValueError(f'Gene {gene!r} is assigned to multiple families')
#             gene_to_family[gene] = family
#     return gene_to_family


# def read_fasta_ids(path, delimiter='|'):
#     by_genome = defaultdict(set)
#     with open_text(path) as handle:
#         for line_number, line in enumerate(handle, start=1):
#             if not line.startswith('>'):
#                 continue
#             identifier = line[1:].strip().split(maxsplit=1)[0]
#             if delimiter not in identifier:
#                 raise ValueError(
#                     f'FASTA identifier {identifier!r} lacks genome delimiter {delimiter!r}; '
#                     'use --sequence-delimiter or provide genome-qualified IDs'
#                 )
#             genome_ids = by_genome[identifier.split(delimiter, 1)[0]]
#             if identifier in genome_ids:
#                 raise ValueError(f'Duplicate FASTA identifier {identifier!r} at line {line_number}')
#             genome_ids.add(identifier)
#     return by_genome


def gff_path_for_genome(gff_dir, genome):
    candidates = [
        os.path.join(gff_dir, genome + suffix)
        for suffix in ('.gff', '.gff3', '.gff.gz', '.gff3.gz')
    ]
    existing = [path for path in candidates if os.path.exists(path)]
    if len(existing) > 1:
        raise ValueError(f'Multiple GFF files found for {genome}: {existing}')
    return existing[0] if existing else None


def parse_attributes(text):
    attributes = {}
    for token in (item.strip() for item in text.split(';') if item.strip()):
        if '=' in token:
            key, value = token.split('=', 1)
        elif ' ' in token:
            key, value = token.split(None, 1)
            value = value.strip('"')
        else:
            continue
        attributes[key] = value
    return attributes


def load_ordered_loci(gff_path, genome, gene_to_family, delimiter='|'):
    contigs = defaultdict(list)
    with open_text(gff_path) as handle:
        for line in handle:
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'CDS':
                continue
            attributes = parse_attributes(fields[8])
            local_ids = [
                attributes.get(key) for key in ('ID', 'locus_tag', 'protein_id', 'gene', 'Name', 'Parent')
                if attributes.get(key)
            ]
            qualified = set()
            prefix = genome + delimiter
            for local_id in local_ids:
                for token in local_id.split(","):
                    candidates = [token]

                    #if not token.startswith(genome + delimiter):
                    #    candidates.append(f"{genome}{delimiter}{token}")
                    
                    if not token.startswith(prefix):
                        candidates.append(prefix + token)

                    for candidate in candidates:
                        family = gene_to_family.get(candidate)
                        if family is not None:
                            qualified.add(candidate)
                    #for candidate in candidates:
                    #    if candidate in gene_to_family: # and candidate in fasta_ids:
                    #        qualified.add(candidate)
                    #        qualified.add(candidate)
            #qualified = list(dict.fromkeys(qualified))

            if len(qualified) > 1:
                raise ValueError(
                    f'GFF feature {fields[0]}:{fields[3]}-{fields[4]} maps to multiple genes: {qualified}'
                )
            if not qualified:
                continue
            #contigs[fields[0]].append((int(fields[3]), int(fields[4]), list(qualified)[0]))
            contigs[fields[0]].append((int(fields[3]), int(fields[4]), next(iter(qualified))))
    #for loci in contigs.values():
    #    loci.sort(key=lambda item: (item[0], item[1], item[2]))
    return contigs


def main():
    parser = argparse.ArgumentParser(description='Build GFF-coordinate-ordered PyamilySeq family paths')
    #parser.add_argument('--fasta', required=True)
    parser.add_argument('--gene-map', required=True, help='Two-column sequence_id-to-family table')
    parser.add_argument('--gff-dir', required=True)
    parser.add_argument('--species-calls', required=True)
    parser.add_argument('--lineages', required=True)
    parser.add_argument('--rank', default='genus', choices=RANK_INDEX)
    parser.add_argument('--sequence-delimiter', default='|')
    parser.add_argument('--output', required=True)
    parser.add_argument('--mapping-report', help='Optional per-genome mapping QC TSV')
    args = parser.parse_args()

    lineages = load_lineages(args.lineages)
    genome_groups = load_genome_groups(args.species_calls, lineages, args.rank)
    #gene_to_family = read_gene_families(args.gene_map, args.sequence_delimiter)
    gene_to_family, gene_counts = read_gene_families(
        args.gene_map,
        args.sequence_delimiter,
    )
    #fasta_by_genome = read_fasta_ids(args.fasta, args.sequence_delimiter)
    report_path = args.mapping_report or args.output + '.mapping.tsv'

    with open(args.output, 'w', newline='', encoding='utf-8') as output, open(
        report_path, 'w', newline='', encoding='utf-8'
    ) as report:
        writer = csv.writer(output)
        report_writer = csv.writer(report, delimiter='\t', lineterminator='\n')
        writer.writerow(['Genome', 'Group', 'Ordered_Families'])
        report_writer.writerow([
            'genome', 'fasta_genes', 'genes_in_family_map', 'genes_with_gff_coordinates',
            'unmapped_to_gff', 'gff_path',
        ])
        for genome in sorted(gene_counts):
            gff_path = gff_path_for_genome(args.gff_dir, genome)
            #mapped_genes = {gene for gene in fasta_by_genome[genome] if gene in gene_to_family}
            # mapped_count = 0
            # for gene in gene_counts[genome]:
            #     if gene in gene_to_family:
            #         mapped_count += 1
            mapped_count = gene_counts[genome]
            if not gff_path:
                report_writer.writerow([
                    genome, gene_counts[genome], mapped_count,
                    0, mapped_count, '',
                ])
                continue
            genome_gene_map = gene_to_family[genome]
            contigs = load_ordered_loci(
                gff_path,
                genome,
                genome_gene_map,
                #fasta_by_genome[genome],
                args.sequence_delimiter,
            )
            tokens = []
            coordinate_genes = 0
            for contig_index, contig in enumerate(sorted(contigs)):
                if contig_index and tokens:
                    tokens.append('|')
                for _, _, gene in contigs[contig]:
                    #tokens.append(gene_to_family[gene])
                    tokens.append(genome_gene_map[gene])
                    coordinate_genes += 1
            writer.writerow([genome, genome_groups.get(genome, 'NA'), ' '.join(tokens)])
            report_writer.writerow([
                genome,
                gene_counts[genome],
                mapped_count,
                coordinate_genes,
                mapped_count - coordinate_genes,
                gff_path,
            ])


if __name__ == '__main__':
    main()