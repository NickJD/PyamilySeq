"""Input validation helpers for PyamilySeq workflows."""

import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict

try:
    from .location_tracking import (
        as_feature_list,
        map_family_members_to_features,
        normalise_member_to_gid,
        parse_clusters,
        parse_gene_types,
        parse_gff,
        parse_gff_dir,
    )
    from .utils import read_fasta
except (ModuleNotFoundError, ImportError):
    from location_tracking import (
        as_feature_list,
        map_family_members_to_features,
        normalise_member_to_gid,
        parse_clusters,
        parse_gene_types,
        parse_gff,
        parse_gff_dir,
    )
    from utils import read_fasta


ISSUE_FIELDS = ['severity', 'code', 'subject', 'message', 'detail']


def add_issue(issues, severity, code, message, subject='', detail=''):
    issues.append({
        'severity': severity,
        'code': code,
        'subject': subject,
        'message': message,
        'detail': detail,
    })


def fasta_headers(path):
    headers = []
    opener = open
    if str(path).lower().endswith('.gz'):
        import gzip
        opener = gzip.open
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('>'):
                headers.append(line[1:].strip().split()[0])
    return headers


def member_tokens(families):
    return [member for members in families.values() for member in members]


def validate_clusters(path, issues, cluster_format='auto'):
    try:
        families = parse_clusters(path, cluster_format=cluster_format)
    except (ValueError, OSError) as error:
        add_issue(issues, 'ERROR', 'CLUSTER_PARSE_ERROR', str(error), path)
        return {}
    if not families:
        add_issue(issues, 'ERROR', 'CLUSTERS_EMPTY', 'No families could be parsed from the cluster file.', path)
        return families

    empty = [family for family, members in families.items() if not members]
    if empty:
        add_issue(
            issues,
            'WARNING',
            'EMPTY_FAMILIES',
            f'{len(empty)} families have no parsed members.',
            path,
            ','.join(empty[:20]),
        )

    member_to_families = defaultdict(set)
    for family, members in families.items():
        for member in members:
            member_to_families[member].add(family)
    reused = {member: fams for member, fams in member_to_families.items() if len(fams) > 1}
    if reused:
        sample = ';'.join(
            f'{member}:{",".join(sorted(fams)[:5])}'
            for member, fams in list(reused.items())[:20]
        )
        add_issue(
            issues,
            'ERROR',
            'MEMBER_IN_MULTIPLE_FAMILIES',
            f'{len(reused)} member IDs appear in more than one family.',
            path,
            sample,
        )
    return families


def validate_fasta(path, families, issues):
    headers = fasta_headers(path)
    header_counts = Counter(headers)
    duplicate_headers = [header for header, count in header_counts.items() if count > 1]
    if duplicate_headers:
        add_issue(
            issues,
            'ERROR',
            'FASTA_DUPLICATE_HEADERS',
            f'{len(duplicate_headers)} FASTA headers are duplicated.',
            path,
            ','.join(duplicate_headers[:20]),
        )

    sequences = read_fasta(path)
    cluster_members = set(member_tokens(families))
    fasta_ids = set(sequences)
    fasta_by_gene = defaultdict(set)
    for fasta_id in fasta_ids:
        fasta_by_gene[normalise_member_to_gid(fasta_id)].add(fasta_id)
    missing = [
        member for member in sorted(cluster_members)
        if member not in fasta_ids
        and len(fasta_by_gene.get(normalise_member_to_gid(member), ())) != 1
    ]
    matched_fasta = set()
    for member in cluster_members:
        if member in fasta_ids:
            matched_fasta.add(member)
        else:
            aliases = fasta_by_gene.get(normalise_member_to_gid(member), set())
            if len(aliases) == 1:
                matched_fasta.update(aliases)
    extra = sorted(fasta_ids - matched_fasta)
    if missing:
        add_issue(
            issues,
            'ERROR',
            'CLUSTER_MEMBERS_MISSING_FROM_FASTA',
            f'{len(missing)} cluster members were not found in the FASTA headers.',
            path,
            ','.join(missing[:30]),
        )
    if extra:
        add_issue(
            issues,
            'WARNING',
            'FASTA_SEQUENCES_NOT_IN_CLUSTERS',
            f'{len(extra)} FASTA sequences were not referenced and will require singleton backfill.',
            path,
            ','.join(extra[:30]),
        )
    return {
        'fasta_headers': len(headers),
        'fasta_unique_headers': len(fasta_ids),
        'fasta_duplicate_headers': len(duplicate_headers),
        'cluster_members_missing_from_fasta': len(missing),
        'fasta_sequences_not_in_clusters': len(extra),
        'predicted_sequences_assigned_after_backfill': len(fasta_ids),
    }


def validate_gff(path=None, gff_dir=None, gene_types=('CDS', 'gene'), families=None, issues=None,
                 feature_layer='auto', min_mapping_confidence=0.0, minimum_mapping_rate=0.0):
    issues = issues if issues is not None else []
    if gff_dir:
        contigs, id_to_feature, gff_genomes = parse_gff_dir(
            gff_dir, gene_types=gene_types, feature_layer=feature_layer,
        )
        subject = gff_dir
    else:
        contigs, id_to_feature = parse_gff(path, gene_types=gene_types, feature_layer=feature_layer)
        gff_genomes = []
        subject = path

    feature_count = sum(len(features) for features in contigs.values())
    ambiguous_aliases = [
        key for key, value in id_to_feature.items()
        if len(as_feature_list(value)) > 1
    ]
    if not feature_count:
        add_issue(issues, 'ERROR', 'GFF_NO_FEATURES', 'No matching GFF features were parsed.', subject)
    if ambiguous_aliases:
        add_issue(
            issues,
            'WARNING',
            'GFF_AMBIGUOUS_ALIASES',
            f'{len(ambiguous_aliases)} GFF ID aliases map to multiple features.',
            subject,
            ','.join(sorted(ambiguous_aliases)[:30]),
        )

    mapping_summary = {}
    if families:
        _, _, diagnostics = map_family_members_to_features(
            families,
            id_to_feature,
            return_diagnostics=True,
            min_confidence=min_mapping_confidence,
        )
        status_counts = Counter(row['status'] for row in diagnostics)
        total = len(diagnostics)
        mapped = status_counts.get('mapped', 0)
        mapping_rate = mapped / total if total else 0.0
        mapping_summary = {
            'member_mapping_total': total,
            'member_mapping_mapped': mapped,
            'member_mapping_ambiguous': status_counts.get('ambiguous', 0),
            'member_mapping_unmapped': status_counts.get('unmapped', 0),
            'member_mapping_low_confidence': status_counts.get('low_confidence', 0),
            'member_mapping_rate': mapping_rate,
            'mean_mapping_confidence': (
                sum(float(row.get('confidence') or 0) for row in diagnostics if row['status'] == 'mapped') / mapped
                if mapped else 0.0
            ),
        }
        if status_counts.get('ambiguous', 0):
            sample = [row['member_id'] for row in diagnostics if row['status'] == 'ambiguous'][:30]
            add_issue(
                issues,
                'WARNING',
                'GFF_MAPPING_AMBIGUOUS',
                f"{status_counts['ambiguous']} cluster members mapped ambiguously to GFF features.",
                subject,
                ','.join(sample),
            )
        if status_counts.get('unmapped', 0):
            sample = [row['member_id'] for row in diagnostics if row['status'] == 'unmapped'][:30]
            add_issue(
                issues,
                'WARNING',
                'GFF_MAPPING_UNMAPPED',
                f"{status_counts['unmapped']} cluster members could not be mapped to GFF features.",
                subject,
                ','.join(sample),
            )
        if status_counts.get('low_confidence', 0):
            sample = [row['member_id'] for row in diagnostics if row['status'] == 'low_confidence'][:30]
            add_issue(
                issues,
                'WARNING',
                'GFF_MAPPING_LOW_CONFIDENCE',
                f"{status_counts['low_confidence']} cluster members were below the mapping-confidence threshold.",
                subject,
                ','.join(sample),
            )
        if mapping_rate < minimum_mapping_rate:
            add_issue(
                issues,
                'ERROR',
                'GFF_MAPPING_RATE_BELOW_MINIMUM',
                f'Mapping rate {mapping_rate:.6f} is below required minimum {minimum_mapping_rate:.6f}.',
                subject,
            )

    summary = {
        'gff_contigs': len(contigs),
        'gff_features': feature_count,
        'gff_genomes': len(gff_genomes),
        'gff_ambiguous_aliases': len(ambiguous_aliases),
    }
    summary.update(mapping_summary)
    return summary


def severity_counts(issues):
    counts = Counter(issue['severity'] for issue in issues)
    return {
        'ERROR': counts.get('ERROR', 0),
        'WARNING': counts.get('WARNING', 0),
        'INFO': counts.get('INFO', 0),
    }


def write_validation_outputs(outdir, summary, issues):
    os.makedirs(outdir, exist_ok=True)
    json_path = os.path.join(outdir, 'validation_report.json')
    txt_path = os.path.join(outdir, 'validation_report.txt')
    issues_path = os.path.join(outdir, 'validation_issues.tsv')

    with open(json_path, 'w', encoding='utf-8') as handle:
        json.dump({'summary': summary, 'issues': issues}, handle, indent=2, sort_keys=True)

    with open(txt_path, 'w', encoding='utf-8') as handle:
        handle.write('PyamilySeq validation report\n')
        for key, value in summary.items():
            handle.write(f'{key}\t{value}\n')
        handle.write('\nIssues\n')
        if not issues:
            handle.write('No issues detected.\n')
        for issue in issues:
            handle.write(
                f"{issue['severity']}\t{issue['code']}\t{issue['subject']}\t"
                f"{issue['message']}\t{issue['detail']}\n"
            )

    with open(issues_path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=ISSUE_FIELDS, delimiter='\t')
        writer.writeheader()
        for issue in issues:
            writer.writerow(issue)
    return json_path, txt_path, issues_path


def validate_inputs(clusters, outdir, fasta=None, gff=None, gff_dir=None, gene_types=('CDS', 'gene'),
                    cluster_format='auto', feature_layer='auto', min_mapping_confidence=0.0,
                    minimum_mapping_rate=0.0):
    issues = []
    families = validate_clusters(clusters, issues, cluster_format=cluster_format)
    summary = {
        'clusters': os.path.abspath(clusters),
        'families': len(families),
        'members': len(member_tokens(families)),
    }

    if families:
        sizes = [len(members) for members in families.values()]
        summary.update({
            'singleton_families': sum(1 for size in sizes if size == 1),
            'largest_family_size': max(sizes) if sizes else 0,
        })

    if fasta:
        summary.update(validate_fasta(fasta, families, issues))

    if gff or gff_dir:
        summary.update(validate_gff(
            gff,
            gff_dir,
            gene_types,
            families,
            issues,
            feature_layer=feature_layer,
            min_mapping_confidence=min_mapping_confidence,
            minimum_mapping_rate=minimum_mapping_rate,
        ))

    summary.update({f'issues_{key.lower()}': value for key, value in severity_counts(issues).items()})
    paths = write_validation_outputs(outdir, summary, issues)
    return summary, issues, paths


def run_from_cli():
    parser = argparse.ArgumentParser(description='Validate PyamilySeq input files before or after a run')
    parser.add_argument('--clusters', required=True, help='Cluster file or gene_presence_absence.csv')
    parser.add_argument(
        '--cluster-format', default='auto',
        choices=['auto', 'gpa', 'family-table', 'gene-table', 'csv', 'gene-table-csv', 'CD-HIT', 'MMseqs'],
        help='Explicit cluster input adapter (default: auto).',
    )
    parser.add_argument('--fasta', help='Optional FASTA file used for clustering')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--gff', help='Optional single GFF/GFF3 file')
    group.add_argument('--gff-dir', help='Optional directory of per-genome GFF/GFF3 files')
    parser.add_argument('--gene-types', default='CDS,gene',
                        help="Comma-separated GFF feature types to use, or 'all' (default: CDS,gene)")
    parser.add_argument('--feature-layer', default='auto', choices=['auto', 'CDS', 'gene', 'all'],
                        help='Locus layer to validate after parent/child collapse (default: auto).')
    parser.add_argument('--min-mapping-confidence', type=float, default=0.0)
    parser.add_argument('--minimum-mapping-rate', type=float, default=0.0)
    parser.add_argument('--outdir', required=True, help='Directory for validation outputs')
    parser.add_argument('--fail-on', choices=['error', 'warning', 'never'], default='error',
                        help='Exit non-zero on errors, warnings, or never (default: error)')
    args = parser.parse_args()
    if not 0 <= args.min_mapping_confidence <= 1:
        parser.error('--min-mapping-confidence must be between 0 and 1')
    if not 0 <= args.minimum_mapping_rate <= 1:
        parser.error('--minimum-mapping-rate must be between 0 and 1')

    summary, issues, paths = validate_inputs(
        args.clusters,
        args.outdir,
        fasta=args.fasta,
        gff=args.gff,
        gff_dir=args.gff_dir,
        gene_types=parse_gene_types(args.gene_types),
        cluster_format=args.cluster_format,
        feature_layer=args.feature_layer,
        min_mapping_confidence=args.min_mapping_confidence,
        minimum_mapping_rate=args.minimum_mapping_rate,
    )
    sys.stderr.write(f'Wrote validation outputs: {", ".join(paths)}\n')
    counts = severity_counts(issues)
    if args.fail_on == 'error' and counts['ERROR']:
        sys.exit(1)
    if args.fail_on == 'warning' and (counts['ERROR'] or counts['WARNING']):
        sys.exit(1)
    return summary, issues


if __name__ == '__main__':
    run_from_cli()
