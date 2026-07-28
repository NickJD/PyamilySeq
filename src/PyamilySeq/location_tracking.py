"""
Location tracking for gene families.

This module maps clustered family members to GFF gene coordinates and computes
per-family location and synteny metrics (ordinal position, relative position,
neighbourhood conservation, orientation consistency, etc.).

Design assumptions (reasonable defaults):
- Cluster file: tab-delimited file with two columns: family_id <TAB> member_list
  where member_list is comma/semicolon/whitespace separated. Member identifiers
  may be of the form "genome|geneid", "genome:geneid", or just "geneid".
- GFF parsing: extract gene features of types 'CDS' or 'gene' by default and
  prefer the attribute key 'ID' or 'Name' as the gene identifier.

This module is intentionally self-contained and aims to reuse the project's
GFF parsing conventions where possible.
"""

from collections import defaultdict, Counter
import csv
import gzip
import hashlib
import html
import os
import statistics
import argparse
import re
from urllib.parse import unquote


PRIMARY_ID_KEYS = (
    'ID', 'Name', 'locus_tag', 'protein_id', 'gene', 'gene_id',
    'transcript_id', 'Parent', 'Derives_from',
)
UNKNOWN_NEIGHBOR_FAMILY = 'UNASSIGNED'


def sanitize_id(s):
    """Return a cleaned identifier token from a raw attribute string.

    - Split on common separators (;, comma, whitespace) and pick the first meaningful token.
    - Remove common prefixes like 'gene:', 'transcript:', 'protein:', 'ID='.
    - Strip quotes, whitespace and semicolons.
    - Return empty string for falsy input.
    """
    if not s:
        return ''
    s0 = unquote(str(s))
    # split on common separators and pick candidate tokens
    parts = re.split(r'[,;\s]+', s0)
    # prefer a token without '=' and that contains alphanumerics
    for p in parts:
        if not p:
            continue
        if '=' in p:
            continue
        p = p.strip()
        # remove common prefixes (case-insensitive)
        pl = p.lower()
        for pref in ('gene:', 'transcript:', 'protein:', 'id=', 'ensb:', 'ensb_'):
            if pl.startswith(pref):
                # strip the prefix length from original p to preserve case in rest
                p = p[len(pref):]
                pl = p.lower()
        p = p.strip(" '\";")
        if p:
            return p
    # fallback to first non-empty trimmed part
    for p in parts:
        if p and p.strip():
            return p.strip(" '\";")
    return ''


def parse_gff_attributes(attributes):
    """Parse a GFF/GTF attributes column into a dict of decoded values."""
    attr_dict = {}
    for token in [t.strip() for t in attributes.split(';') if t.strip()]:
        if '=' in token:
            key, value = token.split('=', 1)
        elif ' ' in token:
            key, value = token.split(None, 1)
            value = value.strip('"')
        else:
            key, value = token, ''
        key = unquote(key.strip())
        value = unquote(value.strip())
        if key in attr_dict and value:
            attr_dict[key] = attr_dict[key] + ',' + value
        else:
            attr_dict[key] = value
    return attr_dict


def id_variants(identifier):
    """Return conservative ID aliases used by common GFF and clustering tools."""
    variants = []

    def add(value):
        if value is None:
            return
        value = unquote(str(value)).strip().strip(" '\";")
        if not value:
            return
        if value not in variants:
            variants.append(value)
        cleaned = sanitize_id(value)
        if cleaned and cleaned not in variants:
            variants.append(cleaned)

    add(identifier)
    for base in list(variants):
        low = base.lower()
        for prefix in ('gene:', 'transcript:', 'protein:', 'id=', 'gene-', 'cds-', 'rna-', 'transcript-', 'protein-'):
            if low.startswith(prefix):
                add(base[len(prefix):])
        if 'ENSB_' in base:
            add(base.replace('ENSB_', 'ENSB:'))
        if 'ENSB:' in base:
            add(base.replace('ENSB:', 'ENSB_'))
        if '|' in base:
            add(base.split('|')[-1])
        version_match = re.match(r'^(.+)\.\d+$', base)
        if version_match:
            add(version_match.group(1))
    return variants


def canonical_id(identifier):
    """Canonical form for last-resort unique matching."""
    if not identifier:
        return ''
    s0 = sanitize_id(identifier).lower()
    for prefix in ('gene:', 'transcript:', 'protein:', 'ensb:', 'ensb_', 'gene-', 'cds-', 'rna-', 'transcript-', 'protein-'):
        if s0.startswith(prefix):
            s0 = s0[len(prefix):]
            break
    return re.sub('[^0-9a-z]+', '', s0)


def feature_identity(feature):
    return (
        feature.get('contig'),
        feature.get('start'),
        feature.get('end'),
        feature.get('strand'),
        feature.get('id'),
    )


def as_feature_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def add_feature_alias(index, key, feature):
    """Add an ID alias while preserving duplicate/ambiguous mappings."""
    if not key:
        return
    for alias in id_variants(key):
        existing = index.get(alias)
        if existing is None:
            index[alias] = feature
            continue
        existing_features = as_feature_list(existing)
        if feature_identity(feature) in {feature_identity(item) for item in existing_features}:
            continue
        index[alias] = existing_features + [feature]


def split_family_label(label):
    if not label or label == UNKNOWN_NEIGHBOR_FAMILY:
        return []
    return [part for part in str(label).split('|') if part and part != UNKNOWN_NEIGHBOR_FAMILY]


def join_tokens(tokens):
    return ';'.join([str(x).strip() for x in tokens if str(x).strip()])


def _gff_basename(path):
    name = os.path.basename(str(path))
    for suffix in ('.gff3.gz', '.gff.gz', '.gff3', '.gff'):
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return os.path.splitext(name)[0]


def _attribute_values(attr_dict, keys=PRIMARY_ID_KEYS):
    values = []
    for key in keys:
        value = attr_dict.get(key)
        if not value:
            continue
        for token in re_split_members(str(value)):
            if token and token not in values:
                values.append(token)
    return values


def _collapse_locus_features(features, feature_layer='auto', parent_by_id=None):
    """Collapse parent/child GFF records and multi-part CDS records into loci."""

    if feature_layer == 'all':
        return list(features)
    parent_by_id = parent_by_id or {}

    def locus_root(feature, index):
        attrs = feature.get('attr_dict', {})
        parent = sanitize_id(
            attrs.get('Parent') or attrs.get('Derives_from') or attrs.get('transcript_id') or ''
        )
        own_id = sanitize_id(
            attrs.get('ID') or attrs.get('gene_id') or attrs.get('locus_tag')
            or feature.get('id') or ''
        )
        key = parent or own_id
        seen = set()
        while key and key not in seen and parent_by_id.get(key):
            seen.add(key)
            key = parent_by_id[key]
        return key or f"coordinate:{feature['start']}:{feature['end']}:{index}"

    groups = defaultdict(list)
    order = []
    for index, feature in enumerate(features):
        key = locus_root(feature, index)
        if key not in groups:
            order.append(key)
        groups[key].append(feature)

    collapsed = []
    preference = {'CDS': 0, 'gene': 1, 'mRNA': 2, 'transcript': 2}
    for key in order:
        candidates = groups[key]
        if feature_layer not in ('auto', None):
            selected = [feature for feature in candidates if feature.get('type') == feature_layer]
            if not selected:
                continue
        else:
            best_type = min(
                (feature.get('type') for feature in candidates),
                key=lambda value: (preference.get(value, 50), value),
            )
            selected = [feature for feature in candidates if feature.get('type') == best_type]
        chosen = dict(selected[0])
        chosen['start'] = min(feature['start'] for feature in selected)
        chosen['end'] = max(feature['end'] for feature in selected)
        chosen['length'] = chosen['end'] - chosen['start'] + 1
        chosen['feature_parts'] = len(selected)
        chosen['collapsed_feature_types'] = sorted({feature.get('type') for feature in candidates})
        aliases = []
        for feature in candidates:
            for alias in feature.get('possible_ids', []):
                if alias and alias not in aliases:
                    aliases.append(alias)
        chosen['possible_ids'] = aliases
        collapsed.append(chosen)
    return collapsed


def _annotate_contig_features(contigs):
    for contig, feats in contigs.items():
        feats.sort(key=lambda feature: (feature['start'], feature['end'], feature.get('id', '')))
        count = len(feats)
        for index, feature in enumerate(feats):
            feature['ordinal'] = index
            feature['contig_count'] = count
            ordinal_relative = index / (count - 1) if count > 1 else 0.5
            contig_length = feature.get('contig_length')
            if contig_length:
                midpoint = (feature['start'] + feature['end']) / 2.0
                feature['relative_index'] = midpoint / contig_length
                feature['position_basis'] = 'bp'
            else:
                feature['relative_index'] = ordinal_relative
                feature['position_basis'] = 'ordinal'
            feature['relative_ordinal_index'] = ordinal_relative


def parse_gff(path, gene_types=('CDS', 'gene'), feature_layer='auto', genome=None,
              contig_genome_map=None, circular_contigs=None):
    """Parse a GFF file and return (contigs, id_to_feature) where:
    - contigs: dict contig -> ordered list of feature dicts
    - id_to_feature: dict gene_id -> feature dict

    Feature dict keys: id, contig, source, type, start, end, strand, phase, attributes
    """
    raw_contigs = defaultdict(list)
    contigs = defaultdict(list)
    id_to_feature = {}
    sequence_regions = {}
    parent_by_id = {}
    user_circular_contigs = set(circular_contigs or ())
    topology_by_contig = {
        contig: ('circular', 'user') for contig in user_circular_contigs
    }
    contig_genome_map = contig_genome_map or {}
    inferred_genome = genome or _gff_basename(path)

    gene_types = None if gene_types in (None, 'all') else set(gene_types)

    opener = gzip.open if str(path).lower().endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8') as fh:
        for line_number, line in enumerate(fh, start=1):
            if line.startswith('##sequence-region'):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        sequence_regions[parts[1]] = int(parts[3])
                    except ValueError:
                        pass
                continue
            if line.startswith('##FASTA'):
                break
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attributes = parts[:9]
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            attr_dict = parse_gff_attributes(attributes)
            own_hierarchy_id = sanitize_id(
                attr_dict.get('ID') or attr_dict.get('transcript_id') or ''
            )
            parent_hierarchy_id = sanitize_id(
                attr_dict.get('Parent') or attr_dict.get('Derives_from')
                or attr_dict.get('gene_id') or ''
            )
            if own_hierarchy_id and parent_hierarchy_id and own_hierarchy_id != parent_hierarchy_id:
                parent_by_id[own_hierarchy_id] = parent_hierarchy_id
            transcript_id = sanitize_id(attr_dict.get('transcript_id') or '')
            gene_id = sanitize_id(attr_dict.get('gene_id') or '')
            if transcript_id and gene_id and transcript_id != gene_id:
                parent_by_id[transcript_id] = gene_id
            if ftype == 'region':
                sequence_regions.setdefault(seqid, end_i)
                raw_circular = attr_dict.get('Is_circular')
                if raw_circular is None:
                    raw_circular = attr_dict.get('is_circular')
                circular_value = str(raw_circular or '').lower()
                if circular_value in {'true', '1', 'yes'}:
                    topology_by_contig.setdefault(seqid, ('circular', 'gff'))
                elif circular_value in {'false', '0', 'no'}:
                    topology_by_contig.setdefault(seqid, ('linear', 'gff'))
            if gene_types is not None and ftype not in gene_types:
                continue
            gid = None
            for prefer in PRIMARY_ID_KEYS:
                if prefer in attr_dict and attr_dict[prefer]:
                    gid = attr_dict[prefer]
                    break
            if not gid:
                gid = f"{seqid}:{start_i}-{end_i}:{strand}"
            orig_gid = gid
            gid = sanitize_id(gid)
            feature = {
                'id': gid,
                'orig_id': orig_gid,
                'contig': seqid,
                'source': source,
                'type': ftype,
                'start': start_i,
                'end': end_i,
                'length': abs(end_i - start_i) + 1,
                'strand': strand,
                'phase': phase,
                'attributes': attributes,
                'attr_dict': attr_dict,
                'genome': contig_genome_map.get(seqid, inferred_genome),
                'contig_length': sequence_regions.get(seqid),
                'circular': topology_by_contig.get(seqid, ('unspecified', 'none'))[0] == 'circular',
                'contig_topology': topology_by_contig.get(seqid, ('unspecified', 'none'))[0],
                'topology_source': topology_by_contig.get(seqid, ('unspecified', 'none'))[1],
            }
            possible_ids = []
            for variant in id_variants(gid):
                if variant and variant not in possible_ids:
                    possible_ids.append(variant)
            # Only identifier-bearing attributes become aliases. Product and note
            # text are intentionally excluded to prevent false positive mappings.
            for token in _attribute_values(attr_dict):
                for alias in id_variants(token):
                    if alias and alias not in possible_ids:
                        possible_ids.append(alias)
            coord_id = f"{seqid}:{start_i}-{end_i}:{strand}"
            if coord_id not in possible_ids:
                possible_ids.append(coord_id)
            feature['possible_ids'] = possible_ids
            raw_contigs[seqid].append(feature)

    for contig, features in raw_contigs.items():
        topology, topology_source = topology_by_contig.get(contig, ('unspecified', 'none'))
        for feature in features:
            feature['contig_length'] = sequence_regions.get(contig, feature.get('contig_length'))
            feature['circular'] = topology == 'circular'
            feature['contig_topology'] = topology
            feature['topology_source'] = topology_source
        contigs[contig].extend(_collapse_locus_features(
            features, feature_layer=feature_layer, parent_by_id=parent_by_id,
        ))
    _annotate_contig_features(contigs)
    for features in contigs.values():
        for feature in features:
            for alias in feature.get('possible_ids', []):
                add_feature_alias(id_to_feature, alias, feature)
                add_feature_alias(id_to_feature, f"{feature.get('genome')}|{alias}", feature)
    return contigs, id_to_feature


def parse_gff_dir(gff_dir, gene_types=('CDS', 'gene'), feature_layer='auto', circular_contigs=None):
    """Parse all GFF files in a directory and return combined contigs and id map.

    Each contig key will be prefixed with genome name (basename without extension)
    as 'genome|contig' to avoid collisions. Feature IDs are stored in id_to_feature
    as 'genome|geneid' for exact mapping; a secondary plain geneid->feature is also
    included (last-seen wins) to allow looser matching.

    Returns (combined_contigs, combined_id_map, genome_basenames)
    """
    combined_contigs = defaultdict(list)
    combined_id_map = {}
    import glob
    gff_paths = []
    if os.path.isdir(gff_dir):
        # find .gff, .gff3 files
        for ext in ('*.gff', '*.gff3', '*.gff.gz', '*.gff3.gz'):
            gff_paths.extend(glob.glob(os.path.join(gff_dir, ext)))
    else:
        raise ValueError(f"gff_dir {gff_dir} is not a directory")
    if not gff_paths:
        raise ValueError(f"No GFF/GFF3 files found in {gff_dir}")

    genome_basenames = []
    for gff_path in sorted(gff_paths):
        genome = _gff_basename(gff_path)
        genome_basenames.append(genome)
        per_file_circular = set()
        for contig in circular_contigs or ():
            if '|' not in contig:
                per_file_circular.add(contig)
            elif contig.startswith(genome + '|'):
                per_file_circular.add(contig.split('|', 1)[1])
        contigs, id_map = parse_gff(
            gff_path,
            gene_types=gene_types,
            feature_layer=feature_layer,
            genome=genome,
            circular_contigs=per_file_circular,
        )
        # prefix contig keys and id keys
        for contig, feats in contigs.items():
            new_contig = genome + '|' + contig
            # copy and adjust features
            new_feats = []
            for f in feats:
                nf = dict(f)
                nf['contig'] = new_contig
                nf['genome'] = genome
                # create genome-qualified id key and normalize the feature id
                gid = sanitize_id(nf.get('id') or '')
                if not gid:
                    # fallback: use coordinate id
                    gid = f"{genome}|{nf.get('contig')}:{nf.get('start')}-{nf.get('end')}:{nf.get('strand')}"
                # update the feature dict id to the sanitized gid for consistency
                nf['id'] = gid
                # normalize possible_ids for the feature as sanitized tokens
                raw_poss = nf.get('possible_ids', set())
                norm_poss = []
                for pid in raw_poss:
                    if not pid:
                        continue
                    sp = sanitize_id(pid)
                    token = sp or pid
                    if token not in norm_poss:
                        norm_poss.append(token)
                # always include coordinate id
                coord = f"{nf.get('contig')}:{nf.get('start')}-{nf.get('end')}:{nf.get('strand')}"
                if coord not in norm_poss:
                    norm_poss.append(coord)
                nf['possible_ids'] = norm_poss
                qid = genome + '|' + gid
                add_feature_alias(combined_id_map, qid, nf)
                # also store plain sanitized id to allow fallback; duplicates remain ambiguous
                add_feature_alias(combined_id_map, gid, nf)
                # also add feature possible_ids (many variants) into combined map for robust matching
                for pid in nf.get('possible_ids', set()):
                    if not pid:
                        continue
                    for var in id_variants(pid):
                        add_feature_alias(combined_id_map, var, nf)
                        add_feature_alias(combined_id_map, genome + '|' + var, nf)
                new_feats.append(nf)
            combined_contigs[new_contig].extend(new_feats)

    _annotate_contig_features(combined_contigs)

    return combined_contigs, combined_id_map, genome_basenames


def resolve_genome_qualifiers(families, gff_genomes):
    """For families parsed from a gene_presence_absence.csv, which qualify members as
    'header_name|geneid', map header_name to the closest matching genome basename
    from the GFF directory (gff_genomes list). Returns a new families dict with
    qualified names replaced where resolvable.
    """
    resolved = {}

    def normalize_name(s):
        s0 = s.lower()
        # remove common file extensions
        for ext in ('.gff', '.gff3', '.fa', '.fna', '.fasta'):
            if s0.endswith(ext):
                s0 = s0[:-len(ext)]
        # remove known suffix tokens
        for token in ('_storf-reporter_extended', '_storf-reporter', '_storf', '_extended', '_combined'):
            if token in s0:
                s0 = s0.replace(token, '')
        # collapse punctuation
        s0 = s0.replace('-', '_')
        s0 = s0.replace('..', '.')
        s0 = s0.strip('_ .')
        return s0

    # Only exact and uniquely-normalized names are accepted. Substring matching
    # can silently map G1 to G10 and is therefore unsafe for genomic provenance.
    norm_map = {}
    for g in gff_genomes:
        norm = normalize_name(g)
        norm_map.setdefault(norm, []).append(g)

    def find_best(gen_header):
        # try exact
        if gen_header in gff_genomes:
            return gen_header
        nh = normalize_name(gen_header)
        candidates = norm_map.get(nh, [])
        return candidates[0] if len(candidates) == 1 else None

    for fam, members in families.items():
        new_members = []
        for m in members:
            if '|' in m:
                hdr, gid = m.split('|', 1)
                best = find_best(hdr)
                if best:
                    new_members.append(best + '|' + gid)
                else:
                    # keep original; resolution failed
                    new_members.append(m)
            else:
                new_members.append(m)
        resolved[fam] = new_members
    return resolved


def parse_clusters(path, cluster_format='auto'):
    """Parse a simple cluster file or a Roary/gene_presence_absence.csv file.

    Supported formats:
    - Two-column tab-delimited: family_id <TAB> members_string
    - Roary/PyamilySeq 'gene_presence_absence.csv' where the first column is
      the group name and subsequent genome columns list gene IDs present in
      each genome (cells may contain multiple IDs separated by whitespace/comma/semicolon).

    Returns dict family_id -> list of member strings
    """
    fam = defaultdict(list)
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    first = None
    with opener(path, 'rt', newline='', encoding='utf-8') as fh:
        for line in fh:
            if line.strip() and not line.startswith('#'):
                first = line
                break
    if first is None:
        return {}
    first_cells = next(csv.reader([first]))
    requested = (cluster_format or 'auto').lower().replace('_', '-')
    is_gpa = requested in {'gpa', 'gene-presence-absence', 'roary'} or (
        requested == 'auto'
        and '\t' not in first
        and len(first_cells) > 1
        and first_cells[0].strip().strip('"').lower() == 'gene'
    )
    if requested in {'cd-hit', 'cdhit'} or (requested == 'auto' and first.startswith('>Cluster')):
        try:
            from .families import iter_cdhit_clusters
        except ImportError:
            from families import iter_cdhit_clusters
        return {f"Cluster_{cluster.label}": list(cluster.members) for cluster in iter_cdhit_clusters(path)}

    if is_gpa:
        # Parse as Roary / gene_presence_absence.csv
        with opener(path, 'rt', newline='', encoding='utf-8') as fh:
            reader = csv.reader(fh)
            header = next(reader)
            # Find the index where genome columns start. Commonly Roary has 14 metadata columns
            # and genome columns start at index 14. Try to locate a known metadata column name
            genome_start_idx = None
            # attempt to automatically locate the first genome column by looking for
            # headers that resemble genome names (contains 'gca_', '.ASM', '_combined' or many underscores)
            for i in range(1, len(header)):
                h = header[i].strip()
                if not h:
                    continue
                # heuristics for genome-like header
                if 'gca_' in h.lower() or '.asm' in h.lower() or '_combined' in h.lower():
                    genome_start_idx = i
                    break
                # many underscores and digits
                if re.search(r'[A-Za-z]+_[A-Za-z0-9_]+_gca_', h):
                    genome_start_idx = i
                    break
                if h.count('_') >= 2 and any(c.isdigit() for c in h):
                    genome_start_idx = i
                    break
            # fallback heuristics
            if genome_start_idx is None:
                if len(header) > 14:
                    genome_start_idx = 14
                else:
                    genome_start_idx = 1

            genome_names = [h.strip() for h in header[genome_start_idx:]]
            # Normalize genome names to use as keys (avoid empty names)
            for row in reader:
                if not row:
                    continue
                family = row[0].strip()
                members = []
                # iterate over genome columns with their genome name
                for i, cell in enumerate(row[genome_start_idx:]):
                    genome = genome_names[i] if i < len(genome_names) else f'genome{i}'
                    if not cell:
                        continue
                    parts = re_split_members(cell)
                    parts = [p for p in parts if p]
                    # Qualify each member with its genome so mapping is precise
                    qualified = [genome + '|' + p for p in parts]
                    members.extend(qualified)
                fam[family].extend(member for member in members if member not in fam[family])
        return dict(fam)

    if requested in {'mmseqs', 'mmseqs2'}:
        try:
            from .families import iter_mmseqs_clusters
        except ImportError:
            from families import iter_mmseqs_clusters
        return {cluster.label: list(cluster.members) for cluster in iter_mmseqs_clusters(path)}

    delimiter = ',' if requested in {'csv', 'family-csv', 'gene-table-csv'} else '\t'
    if requested == 'auto':
        delimiter = '\t' if '\t' in first else ','
    with opener(path, 'rt', newline='', encoding='utf-8') as fh:
        rows = csv.reader(fh, delimiter=delimiter)
        header_row = None
        for row in rows:
            if row and not row[0].startswith('#'):
                header_row = row
                break
        if header_row is None:
            return {}
        header = [cell.strip().lower() for cell in header_row]
        reverse = requested in {'gene-table', 'gene-table-csv'}
        skip_header = False
        if len(header) >= 2:
            if header[0] in {'sequence_id', 'sequence', 'gene_id', 'member_id'} and header[1] in {
                'family_id', 'family', 'group', 'group_id'
            }:
                reverse = True
                skip_header = True
            elif header[0] in {'family_id', 'family', 'group', 'group_id'} and header[1] in {
                'sequence_id', 'sequence', 'gene_id', 'member_id', 'members'
            }:
                skip_header = True
        for row_index, row in enumerate(rows):
            if row_index == 0 and skip_header:
                continue
            if len(row) < 2:
                # Whitespace fallback is only safe when CSV/TSV parsing yielded one field.
                row = row[0].split(None, 1) if row else []
            if len(row) < 2:
                continue
            left, right = row[0].strip(), row[1].strip()
            if reverse:
                family, members = right, [left]
            else:
                family, members = left, [member for member in re_split_members(right) if member]
            for member in members:
                if member not in fam[family]:
                    fam[family].append(member)
        return dict(fam)


def re_split_members(s):
    # split by comma/semicolon or whitespace
    import re
    return re.split(r'[,;\s]+', s)


def normalise_member_to_gid(member):
    """Try to extract a plausible gene id from a cluster member id.

    Supports forms like 'genome|gene', 'genome:gene', or raw 'gene'.
    Returns the gene portion (after last pipe or colon) by default.
    """
    if '|' in member:
        return member.split('|')[-1]
    if re.match(r'^[^:]+:\d+-\d+:[+\-.?]$', member):
        return member
    if ':' in member and member.count(':') == 1:
        return member.split(':')[-1]
    if '/' in member:
        return member.split('/')[-1]
    return member


def filter_singletons(families):
    """Remove families that have only a single member.

    Input: families dict mapping family_id -> list of members
    Returns: new dict containing only families with more than one member.

    This helper prevents a NameError when the CLI option to skip singletons
    is enabled and keeps the behaviour explicit and testable.
    """
    if not families:
        return {}
    filtered = {fam: members for fam, members in families.items() if isinstance(members, (list, tuple)) and len(members) > 1}
    return filtered


def _member_genome(member):
    if '|' in member:
        genome, _ = member.split('|', 1)
        return genome
    return None


def _feature_genome(feature):
    if feature.get('genome'):
        return feature.get('genome')
    contig = feature.get('contig', '')
    if '|' in contig:
        return contig.split('|', 1)[0]
    return None


def _build_canonical_feature_index(id_to_feature):
    canonical_index = defaultdict(list)
    seen = defaultdict(set)
    for key, value in id_to_feature.items():
        c = canonical_id(key)
        if not c:
            continue
        for feature in as_feature_list(value):
            identity = feature_identity(feature)
            if identity in seen[c]:
                continue
            canonical_index[c].append((key, feature))
            seen[c].add(identity)
    return canonical_index


def _select_feature(candidates, member, matched_key):
    features = as_feature_list(candidates)
    if not features:
        return None, 'unmapped', 'not_found'
    if len(features) == 1:
        return features[0], 'mapped', ''

    genome = _member_genome(member)
    if genome:
        genome_matches = [feature for feature in features if _feature_genome(feature) == genome]
        if len(genome_matches) == 1:
            return genome_matches[0], 'mapped', ''

    note = f"ambiguous alias {matched_key} matched {len(features)} features"
    return None, 'ambiguous', note


def _diagnostic_row(family, member, status, match_type, matched_key='', feature=None, note='', confidence=0.0):
    row = {
        'family_id': family,
        'member_id': member,
        'status': status,
        'match_type': match_type,
        'matched_key': matched_key,
        'gene_id': '',
        'genome': '',
        'contig': '',
        'start': '',
        'end': '',
        'strand': '',
        'confidence': confidence,
        'note': note,
    }
    if feature:
        row.update({
            'gene_id': feature.get('id', ''),
            'genome': _feature_genome(feature) or '',
            'contig': feature.get('contig', ''),
            'start': feature.get('start', ''),
            'end': feature.get('end', ''),
            'strand': feature.get('strand', ''),
        })
    return row


def _lookup_member_feature(member, id_to_feature, canonical_index, allow_canonical=True):
    gid = normalise_member_to_gid(member)
    keys_to_try = []
    for key in [member, sanitize_id(member), gid, sanitize_id(gid)]:
        if key and key not in keys_to_try:
            keys_to_try.append(key)
    genome = _member_genome(member)
    for variant in id_variants(gid):
        if variant and variant not in keys_to_try:
            keys_to_try.append(variant)
        if genome:
            qualified = genome + '|' + variant
            if qualified not in keys_to_try:
                keys_to_try.insert(0, qualified)

    for key in keys_to_try:
        if key not in id_to_feature:
            continue
        feature, status, note = _select_feature(id_to_feature[key], member, key)
        if status == 'mapped':
            match_type = 'exact' if key == member else 'alias'
            if key == member:
                confidence = 1.0
            elif genome and key.startswith(genome + '|'):
                confidence = 0.9
            else:
                confidence = 0.8
            return feature, status, match_type, key, note, confidence
        return None, status, 'ambiguous', key, note, 0.0

    canonical = canonical_id(gid)
    if allow_canonical and canonical and canonical in canonical_index:
        canonical_matches = canonical_index[canonical]
        features = []
        seen = set()
        for key, feature in canonical_matches:
            identity = feature_identity(feature)
            if identity not in seen:
                features.append(feature)
                seen.add(identity)
        feature, status, note = _select_feature(features, member, canonical)
        if status == 'mapped':
            return feature, status, 'canonical_unique', canonical, note, 0.6
        return None, status, 'canonical_ambiguous', canonical, note, 0.0

    return None, 'unmapped', 'not_found', '', 'no matching GFF feature', 0.0


def map_family_members_to_features(families, id_to_feature, return_diagnostics=False,
                                   min_confidence=0.0, allow_canonical=True):
    """Map family members to GFF features using id matching heuristics.

    Returns: family_occurrences: dict family -> list of feature dicts (references into parsed features)
    and unmapped: dict family->list of member ids not mapped.
    """
    family_occ = {}
    unmapped = defaultdict(list)
    diagnostics = []
    canonical_index = _build_canonical_feature_index(id_to_feature)

    # Clear previous family annotations in case callers reuse parsed features.
    seen_features = set()
    for value in id_to_feature.values():
        for feature in as_feature_list(value):
            identity = feature_identity(feature)
            if identity in seen_features:
                continue
            feature['family_ids'] = []
            feature['family_mapping_confidence'] = {}
            seen_features.add(identity)

    for fam, members in families.items():
        occ = []
        seen_occurrences = set()
        for m in members:
            feature, status, match_type, matched_key, note, confidence = _lookup_member_feature(
                m,
                id_to_feature,
                canonical_index,
                allow_canonical=allow_canonical,
            )
            if feature is not None and confidence < min_confidence:
                note = (
                    f"mapping confidence {confidence:.2f} below required minimum "
                    f"{min_confidence:.2f}"
                )
                feature = None
                status = 'low_confidence'
            if feature:
                identity = feature_identity(feature)
                if identity not in seen_occurrences:
                    occ.append(feature)
                    seen_occurrences.add(identity)
                if fam not in feature.setdefault('family_ids', []):
                    feature['family_ids'].append(fam)
                feature.setdefault('family_mapping_confidence', {})[fam] = max(
                    confidence,
                    feature.get('family_mapping_confidence', {}).get(fam, 0.0),
                )
                diagnostics.append(_diagnostic_row(
                    fam, m, status, match_type, matched_key, feature, note, confidence,
                ))
            else:
                unmapped[fam].append(m)
                diagnostics.append(_diagnostic_row(
                    fam, m, status, match_type, matched_key, None, note, confidence,
                ))
        family_occ[fam] = occ
    if return_diagnostics:
        return family_occ, unmapped, diagnostics
    return family_occ, unmapped


def compute_neighbors_for_contigs(contigs, k=3, strand_aware=True):
    """Precompute neighbor family placeholders for each feature in contigs.

    This helper doesn't know families yet; instead it records neighbor indices for
    each feature: upstream_ids (list of ordinals) and downstream_ids.
    """
    if k < 0:
        raise ValueError('k must be non-negative')
    for contig, feats in contigs.items():
        n = len(feats)
        for i, f in enumerate(feats):
            circular = bool(f.get('circular'))

            def walk(step):
                neighbors = []
                cursor = i + step
                while len(neighbors) < min(k, max(0, n - 1)):
                    if circular and n:
                        cursor %= n
                    elif cursor < 0 or cursor >= n:
                        break
                    if cursor == i:
                        break
                    neighbors.append(feats[cursor])
                    cursor += step
                return neighbors

            left_features = walk(-1)
            right_features = walk(1)
            if strand_aware and f.get('strand') == '-':
                up_features, down_features = right_features, left_features
            else:
                up_features, down_features = left_features, right_features

            def feature_ids(features):
                out = []
                for feature in features:
                    sid = sanitize_id(feature.get('id'))
                    if sid:
                        out.append(sid)
                return out

            def feature_family_labels(features):
                out = []
                for feature in features:
                    families = sorted(feature.get('family_ids', []))
                    out.append('|'.join(families) if families else UNKNOWN_NEIGHBOR_FAMILY)
                return out

            def directional_gap(neighbor, direction):
                contig_length = f.get('contig_length')
                if direction == 'right':
                    if neighbor['start'] > f['end']:
                        return neighbor['start'] - f['end'] - 1
                    if circular and contig_length:
                        return contig_length - f['end'] + neighbor['start'] - 1
                else:
                    if neighbor['end'] < f['start']:
                        return f['start'] - neighbor['end'] - 1
                    if circular and contig_length:
                        return f['start'] + contig_length - neighbor['end'] - 1
                overlap = min(f['end'], neighbor['end']) - max(f['start'], neighbor['start']) + 1
                return -max(0, overlap)

            def relative_strands(features):
                focal = f.get('strand')
                return [
                    'same' if focal in ('+', '-') and feature.get('strand') == focal else
                    'opposite' if focal in ('+', '-') and feature.get('strand') in ('+', '-') else
                    'unknown'
                    for feature in features
                ]

            up_families = feature_family_labels(up_features)
            down_families = feature_family_labels(down_features)
            f['neighbors_up'] = feature_ids(up_features)
            f['neighbors_down'] = feature_ids(down_features)
            f['neighbor_families_up'] = up_families
            f['neighbor_families_down'] = down_families
            f['neighbor_relative_strands_up'] = relative_strands(up_features)
            f['neighbor_relative_strands_down'] = relative_strands(down_features)
            up_direction = 'right' if strand_aware and f.get('strand') == '-' else 'left'
            down_direction = 'left' if strand_aware and f.get('strand') == '-' else 'right'
            f['intergenic_bp_up'] = [directional_gap(feature, up_direction) for feature in up_features]
            f['intergenic_bp_down'] = [directional_gap(feature, down_direction) for feature in down_features]
            f['signed_neighbor_family_context'] = [
                f'-{index}:{label}' for index, label in enumerate(up_families, start=1)
            ] + [
                f'+{index}:{label}' for index, label in enumerate(down_families, start=1)
            ]
            f['boundary_censored_up'] = not circular and len(up_features) < k
            f['boundary_censored_down'] = not circular and len(down_features) < k


def compute_family_metrics(family_occ, k=3, rel_delta=0.05):
    """Compute per-family metrics and return two tables: occurrences_rows, summary_rows

    occurrences_rows: list of dicts for each member occurrence
    summary_rows: list of dicts for each family summary
    """
    occurrences_rows = []
    summary_rows = []

    def consensus_metrics(contexts):
        if not contexts:
            return (), 0.0, 0.0
        counter = Counter(contexts)
        consensus, consensus_count = counter.most_common(1)[0]
        conserved_fraction = consensus_count / len(contexts)
        consensus_set = set(consensus)
        jaccards = []
        for context in contexts:
            current_set = set(context)
            if not current_set and not consensus_set:
                jaccards.append(1.0)
            elif not current_set or not consensus_set:
                jaccards.append(0.0)
            else:
                jaccards.append(len(current_set & consensus_set) / len(current_set | consensus_set))
        return consensus, conserved_fraction, statistics.mean(jaccards) if jaccards else 0.0

    def clean_neighbor_families(tokens):
        families = []
        for token in tokens:
            families.extend(split_family_label(token))
        return families

    def ordered_similarity(context, consensus):
        width = max(len(context), len(consensus))
        if width == 0:
            return 1.0
        return sum(
            (context[index] if index < len(context) else None)
            == (consensus[index] if index < len(consensus) else None)
            for index in range(width)
        ) / width

    for fam, occ in family_occ.items():
        if not occ:
            summary_rows.append({
                'family_id': fam,
                'occurrences': 0,
                'genomes_with_family': 0,
                'contigs_with_family': 0,
                'mean_copies_per_genome': 0,
                'median_relative_position': '',
                'majority_strand': '',
                'orientation_consistency': '',
                'conserved_position_fraction': '',
                'neighbor_conserved_fraction': '',
                'neighbor_jaccard_mean': '',
                'neighbor_family_conserved_fraction': '',
                'neighbor_family_jaccard_mean': '',
                'neighbor_family_assignment_fraction': '',
                'self_neighbor_fraction': '',
                'dominant_neighbor_family_context': '',
                'dominant_signed_family_architecture': '',
                'signed_architecture_conserved_fraction': '',
                'signed_architecture_similarity_mean': '',
                'relative_orientation_conserved_fraction': '',
                'boundary_censored_fraction': '',
                'circular_occurrence_fraction': '',
                'mean_mapping_confidence': '',
                'minimum_mapping_confidence': '',
                'mean_nearest_upstream_gap_bp': '',
                'mean_nearest_downstream_gap_bp': '',
            })
            continue
        occ_rel_idxs = [o.get('relative_index', 0.5) for o in occ]
        median_rel = statistics.median(occ_rel_idxs)
        # strand majority
        strands = [o.get('strand', '.') for o in occ]
        strand_counts = Counter(strands)
        majority_strand, majority_count = strand_counts.most_common(1)[0]
        orientation_consistency = majority_count / len(occ)

        # conserved position fraction
        within = [1 for r in occ_rel_idxs if abs(r - median_rel) <= rel_delta]
        conserved_pos_frac = sum(within) / len(occ)

        genomes = [(_feature_genome(o) or o.get('contig')) for o in occ]

        # neighbor tuple consensus (up+down tuple) on raw gene IDs
        neighbor_tuples = []
        neighbor_family_contexts = []
        signed_family_contexts = []
        relative_orientation_contexts = []
        assigned_neighbor_family_tokens = 0
        total_neighbor_tokens = 0
        self_neighbor_occurrences = 0
        for o in occ:
            up = tuple(o.get('neighbors_up', []))
            down = tuple(o.get('neighbors_down', []))
            neighbor_tuples.append(tuple(list(up) + list(down)))

            fam_up = tuple(o.get('neighbor_families_up', []))
            fam_down = tuple(o.get('neighbor_families_down', []))
            neighbor_family_tokens = clean_neighbor_families(list(fam_up) + list(fam_down))
            neighbor_family_contexts.append(tuple(neighbor_family_tokens))
            signed_family_contexts.append(tuple(o.get('signed_neighbor_family_context', [])))
            relative_orientation_contexts.append(tuple(
                list(o.get('neighbor_relative_strands_up', []))
                + list(o.get('neighbor_relative_strands_down', []))
            ))
            total_neighbor_tokens += len(list(fam_up) + list(fam_down))
            assigned_neighbor_family_tokens += len(neighbor_family_tokens)
            if fam in set(neighbor_family_tokens):
                self_neighbor_occurrences += 1

        _, neighbor_conserved_fraction, neighbor_jaccard_mean = consensus_metrics(neighbor_tuples)
        dominant_family_context, neighbor_family_conserved_fraction, neighbor_family_jaccard_mean = consensus_metrics(neighbor_family_contexts)
        dominant_signed_context, signed_conserved_fraction, _ = consensus_metrics(signed_family_contexts)
        signed_similarity_mean = statistics.mean(
            ordered_similarity(context, dominant_signed_context) for context in signed_family_contexts
        ) if signed_family_contexts else 0.0
        _, relative_orientation_conserved_fraction, _ = consensus_metrics(relative_orientation_contexts)
        neighbor_family_assignment_fraction = (
            assigned_neighbor_family_tokens / total_neighbor_tokens if total_neighbor_tokens else 0.0
        )
        self_neighbor_fraction = self_neighbor_occurrences / len(occ)

        # produce occurrence rows
        for o in occ:
            row = {
                'family_id': fam,
                'gene_id': o.get('id'),
                'raw_id': o.get('orig_id') or o.get('attributes'),
                'genome': _feature_genome(o) or '',
                'contig': o.get('contig'),
                'start': o.get('start'),
                'end': o.get('end'),
                'length': o.get('length'),
                'feature_type': o.get('type'),
                'feature_parts': o.get('feature_parts', 1),
                'strand': o.get('strand'),
                'ordinal_index': o.get('ordinal'),
                'relative_index': o.get('relative_index'),
                'relative_ordinal_index': o.get('relative_ordinal_index'),
                'position_basis': o.get('position_basis'),
                'contig_gene_count': o.get('contig_count'),
                'contig_length_bp': o.get('contig_length') or '',
                'contig_topology': o.get('contig_topology', 'unspecified'),
                'topology_source': o.get('topology_source', 'none'),
                'circular_contig': bool(o.get('circular')),
                'boundary_censored_up': bool(o.get('boundary_censored_up')),
                'boundary_censored_down': bool(o.get('boundary_censored_down')),
                'mapping_confidence': o.get('family_mapping_confidence', {}).get(fam, ''),
                'neighbors_up': join_tokens(o.get('neighbors_up', [])),
                'neighbors_down': join_tokens(o.get('neighbors_down', [])),
                'neighbor_families_up': join_tokens(o.get('neighbor_families_up', [])),
                'neighbor_families_down': join_tokens(o.get('neighbor_families_down', [])),
                'neighbor_relative_strands_up': join_tokens(o.get('neighbor_relative_strands_up', [])),
                'neighbor_relative_strands_down': join_tokens(o.get('neighbor_relative_strands_down', [])),
                'intergenic_bp_up': join_tokens(o.get('intergenic_bp_up', [])),
                'intergenic_bp_down': join_tokens(o.get('intergenic_bp_down', [])),
                'signed_neighbor_family_context': join_tokens(o.get('signed_neighbor_family_context', [])),
            }
            occurrences_rows.append(row)

        mapping_confidences = [
            o.get('family_mapping_confidence', {}).get(fam)
            for o in occ
            if o.get('family_mapping_confidence', {}).get(fam) is not None
        ]
        nearest_up_gaps = [o['intergenic_bp_up'][0] for o in occ if o.get('intergenic_bp_up')]
        nearest_down_gaps = [o['intergenic_bp_down'][0] for o in occ if o.get('intergenic_bp_down')]
        summary = {
            'family_id': fam,
            'occurrences': len(occ),
            'genomes_with_family': len(set(genomes)),
            'contigs_with_family': len(set(o.get('contig') for o in occ)),
            'mean_copies_per_genome': len(occ) / len(set(genomes)) if genomes else 0,
            'median_relative_position': median_rel,
            'majority_strand': majority_strand,
            'orientation_consistency': orientation_consistency,
            'conserved_position_fraction': conserved_pos_frac,
            'neighbor_conserved_fraction': neighbor_conserved_fraction,
            'neighbor_jaccard_mean': neighbor_jaccard_mean,
            'neighbor_family_conserved_fraction': neighbor_family_conserved_fraction,
            'neighbor_family_jaccard_mean': neighbor_family_jaccard_mean,
            'neighbor_family_assignment_fraction': neighbor_family_assignment_fraction,
            'self_neighbor_fraction': self_neighbor_fraction,
            'dominant_neighbor_family_context': join_tokens(dominant_family_context),
            'dominant_signed_family_architecture': join_tokens(dominant_signed_context),
            'signed_architecture_conserved_fraction': signed_conserved_fraction,
            'signed_architecture_similarity_mean': signed_similarity_mean,
            'relative_orientation_conserved_fraction': relative_orientation_conserved_fraction,
            'boundary_censored_fraction': sum(
                bool(o.get('boundary_censored_up')) or bool(o.get('boundary_censored_down')) for o in occ
            ) / len(occ),
            'circular_occurrence_fraction': sum(bool(o.get('circular')) for o in occ) / len(occ),
            'mean_mapping_confidence': statistics.mean(mapping_confidences) if mapping_confidences else '',
            'minimum_mapping_confidence': min(mapping_confidences) if mapping_confidences else '',
            'mean_nearest_upstream_gap_bp': statistics.mean(nearest_up_gaps) if nearest_up_gaps else '',
            'mean_nearest_downstream_gap_bp': statistics.mean(nearest_down_gaps) if nearest_down_gaps else '',
        }
        summary_rows.append(summary)

    return occurrences_rows, summary_rows


def write_csv_rows(path, rows, fieldnames, delimiter=','):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


OCCURRENCE_FIELDS = [
    'family_id', 'gene_id', 'raw_id', 'genome', 'contig', 'start', 'end', 'length',
    'feature_type', 'feature_parts', 'strand', 'ordinal_index', 'relative_index',
    'relative_ordinal_index', 'position_basis', 'contig_gene_count', 'contig_length_bp',
    'contig_topology', 'topology_source', 'circular_contig',
    'boundary_censored_up', 'boundary_censored_down', 'mapping_confidence',
    'neighbors_up', 'neighbors_down', 'neighbor_families_up', 'neighbor_families_down',
    'neighbor_relative_strands_up', 'neighbor_relative_strands_down',
    'intergenic_bp_up', 'intergenic_bp_down', 'signed_neighbor_family_context',
]

SUMMARY_FIELDS = [
    'family_id', 'occurrences', 'genomes_with_family', 'contigs_with_family',
    'mean_copies_per_genome', 'median_relative_position', 'majority_strand',
    'orientation_consistency', 'conserved_position_fraction',
    'neighbor_conserved_fraction', 'neighbor_jaccard_mean',
    'neighbor_family_conserved_fraction', 'neighbor_family_jaccard_mean',
    'neighbor_family_assignment_fraction', 'self_neighbor_fraction',
    'dominant_neighbor_family_context',
    'dominant_signed_family_architecture', 'signed_architecture_conserved_fraction',
    'signed_architecture_similarity_mean', 'relative_orientation_conserved_fraction',
    'boundary_censored_fraction', 'circular_occurrence_fraction',
    'mean_mapping_confidence', 'minimum_mapping_confidence',
    'mean_nearest_upstream_gap_bp', 'mean_nearest_downstream_gap_bp',
]

MAPPING_DIAGNOSTIC_FIELDS = [
    'family_id', 'member_id', 'status', 'match_type', 'matched_key', 'gene_id',
    'genome', 'contig', 'start', 'end', 'strand', 'note',
]

MAPPING_DETAIL_FIELDS = MAPPING_DIAGNOSTIC_FIELDS[:-1] + ['confidence', 'note']


def write_mapping_diagnostics(outdir, diagnostics):
    mapping_path = os.path.join(outdir, 'member_mapping.tsv')
    problem_path = os.path.join(outdir, 'unmapped_members.tsv')
    write_csv_rows(mapping_path, diagnostics, MAPPING_DETAIL_FIELDS, delimiter='\t')
    problem_rows = [row for row in diagnostics if row.get('status') != 'mapped']
    legacy_problem_rows = [
        {field: row.get(field, '') for field in MAPPING_DIAGNOSTIC_FIELDS}
        for row in problem_rows
    ]
    write_csv_rows(problem_path, legacy_problem_rows, MAPPING_DIAGNOSTIC_FIELDS, delimiter='\t')
    return mapping_path, problem_path


def write_run_report(outdir, families, family_occ, diagnostics, contigs):
    report_path = os.path.join(outdir, 'location_tracking_report.txt')
    total_members = sum(len(members) for members in families.values())
    mapped_members = sum(1 for row in diagnostics if row.get('status') == 'mapped')
    ambiguous_members = sum(1 for row in diagnostics if row.get('status') == 'ambiguous')
    unmapped_members = sum(1 for row in diagnostics if row.get('status') == 'unmapped')
    low_confidence_members = sum(1 for row in diagnostics if row.get('status') == 'low_confidence')
    mapped_families = sum(1 for occ in family_occ.values() if occ)
    feature_count = sum(len(features) for features in contigs.values())
    topology_counts = Counter(
        features[0].get('contig_topology', 'unspecified')
        for features in contigs.values() if features
    )
    with open(report_path, 'w', encoding='utf-8') as handle:
        handle.write('PyamilySeq location tracking report\n')
        handle.write(f'families_input\t{len(families)}\n')
        handle.write(f'families_with_mapped_occurrences\t{mapped_families}\n')
        handle.write(f'members_input\t{total_members}\n')
        handle.write(f'members_mapped\t{mapped_members}\n')
        handle.write(f'members_ambiguous\t{ambiguous_members}\n')
        handle.write(f'members_unmapped\t{unmapped_members}\n')
        handle.write(f'members_low_confidence\t{low_confidence_members}\n')
        handle.write(f'gff_contigs\t{len(contigs)}\n')
        handle.write(f'gff_features\t{feature_count}\n')
        handle.write(f'circular_contigs\t{sum(bool(features and features[0].get("circular")) for features in contigs.values())}\n')
        handle.write(f'linear_contigs\t{topology_counts.get("linear", 0)}\n')
        handle.write(f'unspecified_topology_contigs\t{topology_counts.get("unspecified", 0)}\n')
        handle.write(f'contigs_with_bp_lengths\t{sum(bool(features and features[0].get("contig_length")) for features in contigs.values())}\n')
        if total_members:
            handle.write(f'member_mapping_rate\t{mapped_members / total_members:.6f}\n')
            confidences = [float(row.get('confidence') or 0) for row in diagnostics if row.get('status') == 'mapped']
            if confidences:
                handle.write(f'mean_mapping_confidence\t{statistics.mean(confidences):.6f}\n')
    return report_path


def _escape(value):
    return html.escape(str(value), quote=True)


def _family_colour(family_id):
    digest = hashlib.md5(str(family_id).encode('utf-8')).hexdigest()
    hue = int(digest[:6], 16) % 360
    return f'hsl({hue}, 62%, 48%)'


def _split_semicolon(value):
    if value is None:
        return []
    return [token for token in str(value).split(';') if token]


def _short_label(value, limit=18):
    value = str(value)
    if len(value) <= limit:
        return value
    return value[:limit - 1] + '.'


def _write_location_neighbor_summary(outdir, occurrence_rows):
    path = os.path.join(outdir, 'neighbor_family_summary.tsv')
    counts = Counter()
    occurrence_counts = Counter()
    for row in occurrence_rows:
        focal = row.get('family_id', '')
        seen_in_occurrence = set()
        for side in ('neighbor_families_up', 'neighbor_families_down'):
            for token in _split_semicolon(row.get(side)):
                if token == UNKNOWN_NEIGHBOR_FAMILY:
                    continue
                for family in split_family_label(token):
                    counts[(focal, family)] += 1
                    seen_in_occurrence.add(family)
        for family in seen_in_occurrence:
            occurrence_counts[(focal, family)] += 1
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('family_id\tneighbor_family\tneighbor_tokens\toccurrences_with_neighbor\n')
        for (family_id, neighbor_family), count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
            handle.write(
                f'{family_id}\t{neighbor_family}\t{count}\t'
                f'{occurrence_counts[(family_id, neighbor_family)]}\n'
            )
    return path


def _write_architecture_summary(outdir, occurrence_rows):
    path = os.path.join(outdir, 'neighborhood_architectures.tsv')
    counts = Counter(
        (row.get('family_id', ''), row.get('signed_neighbor_family_context', ''))
        for row in occurrence_rows
    )
    totals = Counter(row.get('family_id', '') for row in occurrence_rows)
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('family_id\tsigned_architecture\toccurrences\tfraction\n')
        for (family_id, architecture), count in sorted(
            counts.items(), key=lambda item: (item[0][0], -item[1], item[0][1])
        ):
            denominator = totals[family_id]
            handle.write(
                f'{family_id}\t{architecture}\t{count}\t'
                f'{count / denominator if denominator else 0:.6f}\n'
            )
    return path


def _write_neighborhood_graphml(outdir, occurrence_rows):
    path = os.path.join(outdir, 'family_neighborhood.graphml')
    nodes = set()
    edge_counts = Counter()
    for row in occurrence_rows:
        focal = row.get('family_id', '')
        if not focal:
            continue
        nodes.add(focal)
        for direction, field in (('upstream', 'neighbor_families_up'), ('downstream', 'neighbor_families_down')):
            for offset, token in enumerate(_split_semicolon(row.get(field)), start=1):
                for neighbor in split_family_label(token):
                    nodes.add(neighbor)
                    edge_counts[(focal, neighbor, direction, offset)] += 1
    node_ids = {node: f'n{index}' for index, node in enumerate(sorted(nodes))}
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        handle.write('<graphml xmlns="http://graphml.graphdrawing.org/xmlns">\n')
        handle.write('<key id="label" for="node" attr.name="label" attr.type="string"/>\n')
        handle.write('<key id="direction" for="edge" attr.name="direction" attr.type="string"/>\n')
        handle.write('<key id="offset" for="edge" attr.name="offset" attr.type="int"/>\n')
        handle.write('<key id="count" for="edge" attr.name="count" attr.type="int"/>\n')
        handle.write('<graph id="family_neighborhood" edgedefault="directed">\n')
        for node, node_id in node_ids.items():
            handle.write(f'<node id="{node_id}"><data key="label">{_escape(node)}</data></node>\n')
        for index, ((focal, neighbor, direction, offset), count) in enumerate(sorted(edge_counts.items())):
            handle.write(
                f'<edge id="e{index}" source="{node_ids[focal]}" target="{node_ids[neighbor]}">'
                f'<data key="direction">{direction}</data><data key="offset">{offset}</data>'
                f'<data key="count">{count}</data></edge>\n'
            )
        handle.write('</graph></graphml>\n')
    return path


def _render_location_svg(path, occurrence_rows, summary_rows, family_limit=40, rows_per_family=20):
    summaries = {
        row.get('family_id'): row
        for row in summary_rows
        if row.get('family_id') and int(row.get('occurrences') or 0) > 0
    }
    ranked_families = sorted(
        summaries,
        key=lambda family: (
            float(summaries[family].get('neighbor_family_conserved_fraction') or 0),
            int(summaries[family].get('occurrences') or 0),
            str(family),
        ),
        reverse=True,
    )[:max(1, family_limit)]
    rows_by_family = defaultdict(list)
    for row in occurrence_rows:
        if row.get('family_id') in ranked_families and len(rows_by_family[row.get('family_id')]) < rows_per_family:
            rows_by_family[row.get('family_id')].append(row)

    width = 1180
    row_h = 30
    section_h = 34
    top = 74
    left_label = 28
    center_x = 560
    box_w = 90
    box_h = 18
    gap = 8
    total_rows = sum(len(rows_by_family[family]) for family in ranked_families)
    height = top + len(ranked_families) * section_h + max(1, total_rows) * row_h + 44

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">',
        '<style>'
        'text{font-family:Arial,Helvetica,sans-serif;fill:#1f2933}'
        '.title{font-size:24px;font-weight:700}'
        '.subtitle{font-size:12px;fill:#52606d}'
        '.family{font-size:13px;font-weight:700}'
        '.small{font-size:10px;fill:#52606d}'
        '.gene{font-size:10px;fill:#ffffff}'
        '</style>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text class="title" x="28" y="36">Gene Family Neighborhoods</text>',
        '<text class="subtitle" x="28" y="58">Upstream families are shown to the left, the focal family in the center, downstream families to the right.</text>',
    ]
    y = top
    if not ranked_families:
        svg.append('<text class="family" x="28" y="96">No mapped family locations available.</text>')
    for family in ranked_families:
        summary = summaries[family]
        svg.append(f'<text class="family" x="{left_label}" y="{y + 18}">{_escape(family)}</text>')
        svg.append(
            f'<text class="subtitle" x="{left_label + 180}" y="{y + 18}">'
            f'occurrences={_escape(summary.get("occurrences"))}; '
            f'neighbor-family conservation={_escape(summary.get("neighbor_family_conserved_fraction"))}'
            '</text>'
        )
        y += section_h
        for row in rows_by_family[family]:
            genome = row.get('genome') or ''
            contig = row.get('contig') or ''
            coords = f"{row.get('start')}-{row.get('end')} {row.get('strand')}"
            svg.append(f'<text class="small" x="{left_label}" y="{y + 14}">{_escape(genome)} {_escape(coords)}</text>')
            upstream = list(reversed(_split_semicolon(row.get('neighbor_families_up'))))
            downstream = _split_semicolon(row.get('neighbor_families_down'))
            for idx, neighbor in enumerate(upstream[-3:]):
                x = center_x - (len(upstream[-3:]) - idx) * (box_w + gap)
                colour = '#9aa5b1' if neighbor == UNKNOWN_NEIGHBOR_FAMILY else _family_colour(neighbor)
                svg.append(f'<rect x="{x}" y="{y}" width="{box_w}" height="{box_h}" rx="3" fill="{colour}"><title>{_escape(neighbor)}</title></rect>')
                svg.append(f'<text class="gene" x="{x + box_w / 2}" y="{y + 13}" text-anchor="middle">{_escape(_short_label(neighbor))}</text>')
            colour = _family_colour(family)
            svg.append(f'<rect x="{center_x}" y="{y - 2}" width="{box_w}" height="{box_h + 4}" rx="3" fill="{colour}"><title>{_escape(family)}; {html.escape(contig)}; {html.escape(coords)}</title></rect>')
            svg.append(f'<text class="gene" x="{center_x + box_w / 2}" y="{y + 13}" text-anchor="middle">{_escape(_short_label(family))}</text>')
            for idx, neighbor in enumerate(downstream[:3]):
                x = center_x + box_w + gap + idx * (box_w + gap)
                colour = '#9aa5b1' if neighbor == UNKNOWN_NEIGHBOR_FAMILY else _family_colour(neighbor)
                svg.append(f'<rect x="{x}" y="{y}" width="{box_w}" height="{box_h}" rx="3" fill="{colour}"><title>{_escape(neighbor)}</title></rect>')
                svg.append(f'<text class="gene" x="{x + box_w / 2}" y="{y + 13}" text-anchor="middle">{_escape(_short_label(neighbor))}</text>')
            y += row_h
    svg.append('</svg>\n')
    with open(path, 'w', encoding='utf-8') as handle:
        handle.write('\n'.join(svg))


def write_location_visualisations(outdir, occurrence_rows, summary_rows, family_limit=40, rows_per_family=20):
    os.makedirs(outdir, exist_ok=True)
    svg_path = os.path.join(outdir, 'family_neighborhoods.svg')
    html_path = os.path.join(outdir, 'location_tracking.html')
    neighbor_summary_path = _write_location_neighbor_summary(outdir, occurrence_rows)
    _write_architecture_summary(outdir, occurrence_rows)
    _write_neighborhood_graphml(outdir, occurrence_rows)
    _render_location_svg(svg_path, occurrence_rows, summary_rows, family_limit, rows_per_family)

    mapped_families = sum(1 for row in summary_rows if int(row.get('occurrences') or 0) > 0)
    html_text = f'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>PyamilySeq Location Tracking</title>
  <style>
    body {{ font-family: Arial, Helvetica, sans-serif; margin: 28px; color: #1f2933; background: #f7f9fb; }}
    main {{ max-width: 1200px; margin: 0 auto; }}
    section {{ background: white; border: 1px solid #d9e2ec; border-radius: 6px; padding: 18px; margin: 18px 0; }}
    img {{ max-width: 100%; height: auto; border: 1px solid #e6edf3; }}
    .metrics {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 12px; }}
    .metric {{ background: #eef4fb; border-radius: 6px; padding: 12px; }}
    .metric strong {{ display: block; font-size: 22px; }}
    a {{ color: #1d5f99; }}
  </style>
</head>
<body>
<main>
  <h1>PyamilySeq Location Tracking</h1>
  <section class="metrics">
    <div class="metric"><strong>{mapped_families}</strong>Mapped families</div>
    <div class="metric"><strong>{len(occurrence_rows)}</strong>Mapped occurrences</div>
    <div class="metric"><strong>{family_limit}</strong>Visual family limit</div>
  </section>
  <section>
    <h2>Family Neighborhoods</h2>
    <img src="family_neighborhoods.svg" alt="Gene family neighborhood visualization">
  </section>
  <section>
    <h2>Tables</h2>
    <p><a href="family_locations.csv">family_locations.csv</a>, <a href="family_summary.csv">family_summary.csv</a>, <a href="neighbor_family_summary.tsv">neighbor_family_summary.tsv</a>, <a href="neighborhood_architectures.tsv">neighborhood_architectures.tsv</a>, <a href="family_neighborhood.graphml">family_neighborhood.graphml</a>, <a href="member_mapping.tsv">member_mapping.tsv</a>, <a href="unmapped_members.tsv">unmapped_members.tsv</a></p>
  </section>
</main>
</body>
</html>
'''
    with open(html_path, 'w', encoding='utf-8') as handle:
        handle.write(html_text)
    return html_path, svg_path, neighbor_summary_path


def load_id_map(path):
    """Load a two-column mapping file: cluster_token <TAB or ,> gff_id

    Returns dict mapping cluster_token -> gff_id
    """
    mp = {}
    if not path:
        return mp
    with open(path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if ',' in line and '\t' not in line:
                a, b = [p.strip() for p in line.split(',', 1)]
            else:
                parts = line.split() if '\t' not in line else line.split('\t')
                if len(parts) >= 2:
                    a, b = parts[0].strip(), parts[1].strip()
                else:
                    continue
            mp[a] = b
    return mp


def parse_gene_types(value):
    if not value:
        return ('CDS', 'gene')
    value = value.strip()
    if value.lower() == 'all':
        return None
    return tuple(part.strip() for part in value.split(',') if part.strip())


def parse_name_set(value):
    if not value:
        return set()
    if os.path.isfile(value):
        with open(value, encoding='utf-8') as handle:
            return {
                line.strip().split()[0]
                for line in handle
                if line.strip() and not line.startswith('#')
            }
    return {item.strip() for item in value.split(',') if item.strip()}


def run_from_cli():
    parser = argparse.ArgumentParser(description='Track family locations from clusters and GFFs')
    parser.add_argument('--clusters', required=True, help='Cluster file (family_id<TAB>members) or gene_presence_absence.csv')
    parser.add_argument(
        '--cluster-format', default='auto',
        choices=['auto', 'gpa', 'family-table', 'gene-table', 'csv', 'gene-table-csv', 'CD-HIT', 'MMseqs'],
        help='Explicit cluster input adapter (default: auto).',
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--gff', help='Single combined GFF file to parse (gene coordinates)')
    group.add_argument('--gff-dir', help='Directory containing per-genome GFF files used to build the pangenome')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--id-map', dest='id_map_file', help='Optional two-column file mapping cluster member token -> GFF id (cluster_id,gff_id)')
    parser.add_argument('--k', type=int, default=3, help='Number of neighbors upstream/downstream to consider')
    parser.add_argument('--rel_delta', type=float, default=0.05, help='Relative index delta for conserved position')
    parser.add_argument('--gene-types', default='CDS,gene',
                        help="Comma-separated GFF feature types to use, or 'all' (default: CDS,gene)")
    parser.add_argument('--feature-layer', default='auto', choices=['auto', 'CDS', 'gene', 'all'],
                        help='Locus layer to retain; auto prefers CDS and collapses parent/child records.')
    parser.add_argument('--contig-genome-map',
                        help='Optional two-column contig-to-genome map for a combined GFF.')
    parser.add_argument('--circular-contigs',
                        help='Comma-separated circular contigs or a one-contig-per-line file.')
    parser.add_argument('--genomic-neighbors', action='store_true',
                        help='Report genomic left/right as upstream/downstream instead of strand-relative context.')
    parser.add_argument('--min-mapping-confidence', type=float, default=0.0,
                        help='Minimum accepted mapping confidence from 0 to 1 (default: 0).')
    parser.add_argument('--strict-mapping', action='store_true',
                        help='Disable last-resort canonical-ID mapping.')
    parser.add_argument('--minimum-mapping-rate', type=float, default=0.0,
                        help='Exit non-zero when the mapped member fraction is below this 0-1 threshold.')
    parser.add_argument('--unmapped-family-report-limit', type=int, default=20,
                        help='Maximum families sampled in the stderr unmapped summary (default: 20).')
    parser.add_argument('--unmapped-token-report-limit', type=int, default=20,
                        help='Maximum unmapped member tokens shown on stderr (default: 20).')
    parser.add_argument('--visual-limit', type=int, default=40,
                        help='Maximum number of families to draw in the HTML/SVG neighborhood overview (default: 40)')
    parser.add_argument('--visual-rows-per-family', type=int, default=20,
                        help='Maximum occurrences to draw per family in the neighborhood overview (default: 20)')
    parser.add_argument('--no-visuals', action='store_true',
                        help='Skip HTML/SVG location-tracking visual outputs')
    parser.add_argument('--skip-singletons', action='store_true', default=False,
                        help='Skip families with one member (default: include them).')
    parser.add_argument('--no-skip-singletons', dest='skip_singletons', action='store_false',
                        help=argparse.SUPPRESS)
    args = parser.parse_args()

    if not 0 <= args.min_mapping_confidence <= 1:
        parser.error('--min-mapping-confidence must be between 0 and 1')
    if not 0 <= args.minimum_mapping_rate <= 1:
        parser.error('--minimum-mapping-rate must be between 0 and 1')
    if args.k < 0:
        parser.error('--k must be non-negative')
    if args.unmapped_family_report_limit < 0:
        parser.error('--unmapped-family-report-limit must be zero or greater')
    if args.unmapped_token_report_limit < 0:
        parser.error('--unmapped-token-report-limit must be zero or greater')

    # parse GFF(s)
    gene_types = parse_gene_types(args.gene_types)
    circular_contigs = parse_name_set(args.circular_contigs)
    if args.gff_dir:
        contigs, id_to_feature, gff_genomes = parse_gff_dir(
            args.gff_dir,
            gene_types=gene_types,
            feature_layer=args.feature_layer,
            circular_contigs=circular_contigs,
        )
    else:
        contig_genome_map = load_id_map(args.contig_genome_map) if args.contig_genome_map else None
        contigs, id_to_feature = parse_gff(
            args.gff,
            gene_types=gene_types,
            feature_layer=args.feature_layer,
            contig_genome_map=contig_genome_map,
            circular_contigs=circular_contigs,
        )
        gff_genomes = []
    import sys
    families = parse_clusters(args.clusters, cluster_format=args.cluster_format)
    # If families came from gene_presence_absence.csv and we have gff_genomes, try to resolve headers
    if gff_genomes:
        families = resolve_genome_qualifiers(families, gff_genomes)
    # Optionally drop singleton families to save compute and IO
    if getattr(args, 'skip_singletons', True):
        before = len(families)
        families = filter_singletons(families)
        after = len(families)
        sys.stderr.write(f'Skipped singleton families: {before - after} (remaining {after})\n')
    # If an id-map is provided, apply explicit remapping of cluster tokens to GFF ids
    if getattr(args, 'id_map_file', None):
        idmap = load_id_map(args.id_map_file)
        remap_count = 0
        for fam, members in list(families.items()):
            new_members = []
            for m in members:
                if m in idmap:
                    new_members.append(idmap[m])
                    remap_count += 1
                else:
                    new_members.append(m)
            families[fam] = new_members
        sys.stderr.write(f'Applied id-map replacements: {remap_count}\n')
    family_occ, unmapped, diagnostics = map_family_members_to_features(
        families,
        id_to_feature,
        return_diagnostics=True,
        min_confidence=args.min_mapping_confidence,
        allow_canonical=not args.strict_mapping,
    )
    compute_neighbors_for_contigs(contigs, k=args.k, strand_aware=not args.genomic_neighbors)
    occ_rows, sum_rows = compute_family_metrics(family_occ, k=args.k, rel_delta=args.rel_delta)

    out_occ = os.path.join(args.outdir, 'family_locations.csv')
    out_sum = os.path.join(args.outdir, 'family_summary.csv')
    write_csv_rows(out_occ, occ_rows, OCCURRENCE_FIELDS)
    write_csv_rows(out_sum, sum_rows, SUMMARY_FIELDS)
    mapping_path, problem_path = write_mapping_diagnostics(args.outdir, diagnostics)
    report_path = write_run_report(args.outdir, families, family_occ, diagnostics, contigs)
    sys.stderr.write(f'Wrote mapping diagnostics: {mapping_path}\n')
    sys.stderr.write(f'Wrote unmapped/ambiguous members: {problem_path}\n')
    sys.stderr.write(f'Wrote run report: {report_path}\n')
    if not args.no_visuals:
        visual_paths = write_location_visualisations(
            args.outdir,
            occ_rows,
            sum_rows,
            family_limit=max(1, args.visual_limit),
            rows_per_family=max(1, args.visual_rows_per_family),
        )
        sys.stderr.write(f'Wrote location visual outputs: {", ".join(visual_paths)}\n')

    # print unmapped summary
    if unmapped:
        sys.stderr.write('Unmapped members summary (per-family sample)\n')
        family_items = sorted(unmapped.items())
        family_limit = args.unmapped_family_report_limit
        for f, ms in family_items[:family_limit]:
            sys.stderr.write(f + '\t' + ','.join(ms[:10]) + ('...' if len(ms)>10 else '') + '\n')
        omitted = len(family_items) - min(len(family_items), family_limit)
        if omitted:
            sys.stderr.write(
                f'... {omitted} additional families omitted; see unmapped_members.tsv for all records.\n'
            )
        # aggregate unmapped tokens and print top occurrences
        from collections import Counter
        all_un = []
        for ms in unmapped.values():
            all_un.extend(ms)
        c = Counter(all_un)
        sys.stderr.write(f'\nTop {args.unmapped_token_report_limit} unmapped tokens overall:\n')
        for token, cnt in c.most_common(args.unmapped_token_report_limit):
            sys.stderr.write(f'{token}\t{cnt}\n')

    mapping_rate = (
        sum(row.get('status') == 'mapped' for row in diagnostics) / len(diagnostics)
        if diagnostics else 0.0
    )
    if mapping_rate < args.minimum_mapping_rate:
        sys.stderr.write(
            f'Mapping rate {mapping_rate:.6f} is below required minimum '
            f'{args.minimum_mapping_rate:.6f}\n'
        )
        raise SystemExit(2)


if __name__ == '__main__':
    run_from_cli()
