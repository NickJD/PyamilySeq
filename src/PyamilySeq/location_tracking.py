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
import os
import statistics
import argparse
import re


def sanitize_id(s):
    """Return a cleaned identifier token from a raw attribute string.

    - Split on common separators (;, comma, whitespace) and pick the first meaningful token.
    - Remove common prefixes like 'gene:', 'transcript:', 'protein:', 'ID='.
    - Strip quotes, whitespace and semicolons.
    - Return empty string for falsy input.
    """
    if not s:
        return ''
    s0 = str(s)
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


def parse_gff(path, gene_types=('CDS', 'gene')):
    """Parse a GFF file and return (contigs, id_to_feature) where:
    - contigs: dict contig -> ordered list of feature dicts
    - id_to_feature: dict gene_id -> feature dict

    Feature dict keys: id, contig, source, type, start, end, strand, phase, attributes
    """
    contigs = defaultdict(list)
    id_to_feature = {}

    def extract_id(attr_text):
        # look for common attribute keys
        for key in ('ID', 'Name', 'gene', 'locus_tag', 'protein_id'):
            token = key + '='
            if token in attr_text:
                try:
                    return [p for p in attr_text.split(';') if p.startswith(token)][0].split('=')[1]
                except Exception:
                    continue
        # fallback: return the last semicolon-separated token or the full string
        parts = [p for p in [s.strip() for s in attr_text.split(';')] if p]
        if parts:
            last = parts[-1]
            if '=' in last:
                return last.split('=')[1]
            return last
        return attr_text.strip()

    with open(path, 'r') as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attributes = parts[:9]
            if ftype not in gene_types:
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            # parse attributes into dict and collect candidate ids/values
            attr_dict = {}
            for token in [t for t in attributes.split(';') if t.strip()]:
                if '=' in token:
                    k, v = token.split('=', 1)
                    attr_dict[k] = v
                else:
                    # some GFFs may have flags without '='
                    attr_dict[token] = ''
            # prefer explicit keys from parsed attributes when available, then sanitize
            gid = None
            for prefer in ('ID', 'Name', 'locus_tag', 'protein_id', 'gene'):
                if prefer in attr_dict and attr_dict[prefer]:
                    gid = attr_dict[prefer]
                    break
            if not gid:
                gid = extract_id(attributes)
            # keep original raw id for debugging, then sanitize for canonical use
            orig_gid = gid
            gid = sanitize_id(gid)
            # build feature record and include parsed attributes and candidate ids
            feature = {
                'id': gid,
                'orig_id': orig_gid,
                'contig': seqid,
                'source': source,
                'type': ftype,
                'start': start_i,
                'end': end_i,
                'strand': strand,
                'phase': phase,
                'attributes': attributes,
                'attr_dict': attr_dict,
            }
            # collect possible id tokens (ID, Name, locus_tag, protein_id, and all attr values)
            possible_ids = []
            if gid:
                possible_ids.append(gid)
            for key in ('ID', 'Name', 'locus_tag', 'gene', 'protein_id'):
                if key in attr_dict and attr_dict[key]:
                    # split multi-valued attribute entries and add cleaned tokens
                    for token in re_split_members(str(attr_dict[key])):
                        t = sanitize_id(token)
                        if t and t not in possible_ids:
                            possible_ids.append(t)
            for v in attr_dict.values():
                if v:
                    # sometimes values contain multiple ids separated by commas/spaces
                    for token in re_split_members(str(v)):
                        t = sanitize_id(token)
                        if t and t not in possible_ids:
                            possible_ids.append(t)
            # coordinate-based synthetic id
            coord_id = f"{seqid}:{start_i}-{end_i}:{strand}"
            if coord_id not in possible_ids:
                possible_ids.append(coord_id)
            feature['possible_ids'] = possible_ids
            contigs[seqid].append(feature)
            # allow multiple features with same id; store last (we'll disambiguate by coords later if needed)
            # Store feature under its sanitized gid and also under sanitized possible ids
            if gid:
                id_to_feature[gid] = feature
            for pid in possible_ids:
                if pid:
                    sid = sanitize_id(pid)
                    if sid and sid not in id_to_feature:
                        id_to_feature[sid] = feature
                    # also keep the raw pid for looser matching
                    if pid not in id_to_feature:
                        id_to_feature[pid] = feature

    # sort each contig by start coordinate
    for contig, feats in contigs.items():
        feats.sort(key=lambda f: (f['start'], f['end']))
        # annotate ordinal and contig counts
        count = len(feats)
        for idx, f in enumerate(feats):
            f['ordinal'] = idx
            f['contig_count'] = count
            if count > 1:
                f['relative_index'] = idx / (count - 1)
            else:
                f['relative_index'] = 0.5
            # ensure id is cleaned and normalized for downstream neighbor lists
            if 'id' in f and isinstance(f['id'], str):
                f['id'] = f['id'].strip().strip(" '\";")
            else:
                # fallback to a coordinate-based id if missing
                f['id'] = f.get('id') or f"{f.get('contig')}:{f.get('start')}-{f.get('end')}:{f.get('strand')}"
    return contigs, id_to_feature


def parse_gff_dir(gff_dir, gene_types=('CDS', 'gene')):
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
        for ext in ('*.gff', '*.gff3'):
            gff_paths.extend(glob.glob(os.path.join(gff_dir, ext)))
    else:
        raise ValueError(f"gff_dir {gff_dir} is not a directory")

    genome_basenames = []
    def generate_id_variants(pid):
        vars = set()
        vars.add(pid)
        # common replacements
        if ':' in pid:
            vars.add(pid.replace(':', '_'))
        if '_' in pid:
            vars.add(pid.replace('_', ':'))
        # strip common prefixes
        for p in ('gene:', 'transcript:', 'protein:'):
            if pid.startswith(p):
                core = pid[len(p):]
                vars.add(core)
                vars.add(core.replace(':','_'))
                vars.add(core.replace('_',':'))
        # ENSB variants: ENSB_xxx <-> ENSB:xxx
        if 'ENSB_' in pid:
            vars.add(pid.replace('ENSB_', 'ENSB:'))
        if 'ENSB:' in pid:
            vars.add(pid.replace('ENSB:', 'ENSB_'))
        # also add lower/upper forms
        vars.update({v for v in list(vars)})
        return vars

    for gff_path in sorted(gff_paths):
        genome = os.path.splitext(os.path.basename(gff_path))[0]
        genome_basenames.append(genome)
        contigs, id_map = parse_gff(gff_path, gene_types=gene_types)
        # prefix contig keys and id keys
        for contig, feats in contigs.items():
            new_contig = genome + '|' + contig
            # copy and adjust features
            new_feats = []
            for f in feats:
                nf = dict(f)
                nf['contig'] = new_contig
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
                combined_id_map[qid] = nf
                # also store plain sanitized id to allow fallback
                combined_id_map[gid] = nf
                # also add feature possible_ids (many variants) into combined map for robust matching
                for pid in nf.get('possible_ids', set()):
                    if not pid:
                        continue
                    sp = sanitize_id(pid)
                    for var in generate_id_variants(sp or pid):
                        if var and var not in combined_id_map:
                            combined_id_map[var] = nf
                        gp = genome + '|' + var
                        if gp not in combined_id_map:
                            combined_id_map[gp] = nf
                new_feats.append(nf)
            combined_contigs[new_contig].extend(new_feats)

    # sort and annotate ordinals for each combined contig
    for contig, feats in combined_contigs.items():
        feats.sort(key=lambda f: (f['start'], f['end']))
        count = len(feats)
        for idx, f in enumerate(feats):
            f['ordinal'] = idx
            f['contig_count'] = count
            f['relative_index'] = idx / (count - 1) if count > 1 else 0.5

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

    # precompute normalized forms of gff_genomes
    norm_map = {}
    for g in gff_genomes:
        norm = normalize_name(g)
        norm_map.setdefault(norm, []).append(g)
        # also store shorter tokens split by '.' and '_'
        for part in [p for p in re_split_members(g) if p]:
            p_norm = normalize_name(part)
            norm_map.setdefault(p_norm, []).append(g)

    def find_best(gen_header):
        # try exact
        if gen_header in gff_genomes:
            return gen_header
        nh = normalize_name(gen_header)
        # direct normalized match
        if nh in norm_map:
            # prefer exact g that contains the header as substring if available
            candidates = norm_map[nh]
            if len(candidates) == 1:
                return candidates[0]
            # choose candidate with longest common substring
            best = sorted(candidates, key=lambda x: -len(x))[0]
            return best
        # try substring matching against normalized genome names
        for norm_g, originals in norm_map.items():
            if nh in norm_g or norm_g in nh:
                return originals[0]
        # try fuzzy: look for any gff genome whose basename contains the header token
        for g in gff_genomes:
            if gen_header in g or nh in g.lower():
                return g
        # no match
        return None

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


def parse_clusters(path):
    """Parse a simple cluster file or a Roary/gene_presence_absence.csv file.

    Supported formats:
    - Two-column tab-delimited: family_id <TAB> members_string
    - Roary/PyamilySeq 'gene_presence_absence.csv' where the first column is
      the group name and subsequent genome columns list gene IDs present in
      each genome (cells may contain multiple IDs separated by whitespace/comma/semicolon).

    Returns dict family_id -> list of member strings
    """
    fam = {}
    # Try to detect CSV (Roary) format by inspecting the header
    try:
        with open(path, 'r', newline='') as fh:
            # peek first line
            first = fh.readline()
            if not first:
                return fam
            # Heuristic: if the first line contains commas and contains 'Gene' header,
            # treat as CSV Roary-style
            is_csv = (',' in first) and ('Gene' in first or 'gene' in first)
    except Exception:
        is_csv = False

    if is_csv or path.lower().endswith('.csv'):
        # Parse as Roary / gene_presence_absence.csv
        import csv
        with open(path, 'r', newline='') as fh:
            reader = csv.reader(fh)
            header = next(reader)
            # Find the index where genome columns start. Commonly Roary has 14 metadata columns
            # and genome columns start at index 14. Try to locate a known metadata column name
            genome_start_idx = None
            # attempt to automatically locate the first genome column by looking for
            # headers that resemble genome names (contains 'gca_', '.ASM', '_combined' or many underscores)
            import re
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
                fam[family] = members
        return fam

    # Fallback: simple two-column parser (tab-separated or whitespace)
    with open(path, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) == 1:
                # try space-separated family and members
                parts = line.split(None, 1)
            if len(parts) < 2:
                continue
            family = parts[0]
            members_raw = parts[1]
            # split on commas, semicolons, or whitespace
            raw = [m.strip() for m in re_split_members(members_raw)]
            members = [m for m in raw if m]
            fam[family] = members
    return fam


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
    if ':' in member:
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


def map_family_members_to_features(families, id_to_feature):
    """Map family members to GFF features using id matching heuristics.

    Returns: family_occurrences: dict family -> list of feature dicts (references into parsed features)
    and unmapped: dict family->list of member ids not mapped.
    """
    family_occ = {}
    unmapped = defaultdict(list)
    # Build a canonical index for id_to_feature keys to speed up fuzzy matching.
    # canonical form: lower-case, remove non-alphanumeric characters, strip common prefixes like 'ensb' and 'gene' tokens
    def canonical(s):
        if not s:
            return ''
        s0 = s.lower()
        # remove common prefixes
        for p in ('gene:', 'transcript:', 'protein:', 'ensb:', 'ensb_'):
            if s0.startswith(p):
                s0 = s0[len(p):]
                break
        # remove non-alphanumeric
        import re
        s0 = re.sub('[^0-9a-z]+', '', s0)
        return s0

    canonical_index = defaultdict(list)
    for key in list(id_to_feature.keys()):
        c = canonical(key)
        if c:
            canonical_index[c].append(key)

    for fam, members in families.items():
        occ = []
        for m in members:
            # Try exact match first (supports 'genome|gid' keys when parse_gff_dir was used)
            if m in id_to_feature:
                occ.append(id_to_feature[m])
                continue
            # Normalize to plain gid and try
            gid = normalise_member_to_gid(m)
            if gid in id_to_feature:
                occ.append(id_to_feature[gid])
                continue
            # Try common ENSB underscore -> colon variant (CSV uses ENSB_xxx, GFF often has ENSB:xxx)
            if gid.startswith('ENSB_') or gid.startswith('ENSB-'):
                v = 'ENSB:' + gid.split('_',1)[1] if '_' in gid else gid.replace('ENSB-','ENSB:')
                if v in id_to_feature:
                    occ.append(id_to_feature[v])
                    continue
                # also try without ENSB prefix
                rest = gid.split('_',1)[1] if '_' in gid else gid
                if rest in id_to_feature:
                    occ.append(id_to_feature[rest])
                    continue
            # Try canonicalized matching using the prebuilt index
            cg = canonical(gid)
            if cg and cg in canonical_index:
                # pick first candidate mapped
                candidate_key = canonical_index[cg][0]
                occ.append(id_to_feature[candidate_key])
                continue
            # try fuzzy matching: exact suffix/prefix match of gid against id keys
            matched = None
            for k in id_to_feature.keys():
                if k.endswith(gid) or k.startswith(gid):
                    matched = id_to_feature[k]
                    break
            if matched:
                occ.append(matched)
            else:
                unmapped[fam].append(m)
        family_occ[fam] = occ
    return family_occ, unmapped


def compute_neighbors_for_contigs(contigs, k=3):
    """Precompute neighbor family placeholders for each feature in contigs.

    This helper doesn't know families yet; instead it records neighbor indices for
    each feature: upstream_ids (list of ordinals) and downstream_ids.
    """
    for contig, feats in contigs.items():
        n = len(feats)
        for i, f in enumerate(feats):
            up = []
            down = []
            for j in range(max(0, i - k), i):
                vid = feats[j].get('id')
                sid = sanitize_id(vid)
                if sid:
                    up.append(sid)
            for j in range(i + 1, min(n, i + 1 + k)):
                vid = feats[j].get('id')
                sid = sanitize_id(vid)
                if sid:
                    down.append(sid)
            # deduplicate preserving order
            def dedup_keep_order(seq):
                seen = set()
                out = []
                for x in seq:
                    if x not in seen:
                        seen.add(x)
                        out.append(x)
                return out
            # Remove occurrences of the feature's own id from neighbor lists (handles GFFs with duplicate IDs)
            this_id = sanitize_id(f.get('id'))
            up_clean = [x for x in dedup_keep_order(up) if x and x != this_id]
            down_clean = [x for x in dedup_keep_order(down) if x and x != this_id]
            f['neighbors_up'] = list(reversed(up_clean))  # nearest first
            f['neighbors_down'] = down_clean


def compute_family_metrics(family_occ, k=3, rel_delta=0.05):
    """Compute per-family metrics and return two tables: occurrences_rows, summary_rows

    occurrences_rows: list of dicts for each member occurrence
    summary_rows: list of dicts for each family summary
    """
    occurrences_rows = []
    summary_rows = []

    for fam, occ in family_occ.items():
        if not occ:
            summary_rows.append({
                'family_id': fam,
                'occurrences': 0,
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

        # neighbor tuple consensus (up+down tuple)
        neighbor_tuples = []
        neighbor_sets = []
        for o in occ:
            up = tuple(o.get('neighbors_up', []))
            down = tuple(o.get('neighbors_down', []))
            neighbor_tuples.append((up, down))
            neighbor_sets.append(set(list(up) + list(down)))
        nt_counter = Counter(neighbor_tuples)
        consensus_tuple, consensus_count = nt_counter.most_common(1)[0]
        neighbor_conserved_fraction = consensus_count / len(occ)

        # jaccard mean
        consensus_set = set(list(consensus_tuple[0]) + list(consensus_tuple[1]))
        jaccards = []
        for s in neighbor_sets:
            if not s and not consensus_set:
                j = 1.0
            elif not s or not consensus_set:
                j = 0.0
            else:
                j = len(s & consensus_set) / len(s | consensus_set)
            jaccards.append(j)
        neighbor_jaccard_mean = statistics.mean(jaccards) if jaccards else 0.0

        # produce occurrence rows
        for o in occ:
            row = {
                'family_id': fam,
                'gene_id': o.get('id'),
                'raw_id': o.get('orig_id') or o.get('attributes'),
                # if contig is qualified 'genome|contig', extract genome for easier grouping
                'genome': o.get('contig').split('|')[0] if '|' in o.get('contig','') else None,
                'contig': o.get('contig'),
                'start': o.get('start'),
                'end': o.get('end'),
                'strand': o.get('strand'),
                'ordinal_index': o.get('ordinal'),
                'relative_index': o.get('relative_index'),
                'contig_gene_count': o.get('contig_count'),
                'neighbors_up': ';'.join([str(x).strip() for x in o.get('neighbors_up', [])]),
                'neighbors_down': ';'.join([str(x).strip() for x in o.get('neighbors_down', [])]),
            }
            occurrences_rows.append(row)

        summary = {
            'family_id': fam,
            'occurrences': len(occ),
            'genomes_with_family': len(set([o.get('contig').split('|')[0] if '|' in o.get('contig','') else o.get('contig') for o in occ])),
            'median_relative_position': median_rel,
            'orientation_consistency': orientation_consistency,
            'conserved_position_fraction': conserved_pos_frac,
            'neighbor_conserved_fraction': neighbor_conserved_fraction,
            'neighbor_jaccard_mean': neighbor_jaccard_mean,
        }
        summary_rows.append(summary)

    return occurrences_rows, summary_rows


def write_csv_rows(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


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


def run_from_cli():
    parser = argparse.ArgumentParser(description='Track family locations from clusters and GFFs')
    parser.add_argument('--clusters', required=True, help='Cluster file (family_id<TAB>members) or gene_presence_absence.csv')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--gff', help='Single combined GFF file to parse (gene coordinates)')
    group.add_argument('--gff-dir', help='Directory containing per-genome GFF files used to build the pangenome')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--id-map', dest='id_map_file', help='Optional two-column file mapping cluster member token -> GFF id (cluster_id,gff_id)')
    parser.add_argument('--k', type=int, default=3, help='Number of neighbors upstream/downstream to consider')
    parser.add_argument('--rel_delta', type=float, default=0.05, help='Relative index delta for conserved position')
    parser.add_argument('--no-skip-singletons', dest='skip_singletons', action='store_false', default=True,
                        help='Default: skip families with only a single member for efficiency; set this flag to include singletons')
    args = parser.parse_args()

    # parse GFF(s)
    if args.gff_dir:
        contigs, id_to_feature, gff_genomes = parse_gff_dir(args.gff_dir)
    else:
        contigs, id_to_feature = parse_gff(args.gff)
        gff_genomes = []
    import sys
    families = parse_clusters(args.clusters)
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
    family_occ, unmapped = map_family_members_to_features(families, id_to_feature)
    compute_neighbors_for_contigs(contigs, k=args.k)
    occ_rows, sum_rows = compute_family_metrics(family_occ, k=args.k, rel_delta=args.rel_delta)

    out_occ = os.path.join(args.outdir, 'family_locations.csv')
    out_sum = os.path.join(args.outdir, 'family_summary.csv')
    write_csv_rows(out_occ, occ_rows, ['family_id','gene_id','raw_id','genome','contig','start','end','strand','ordinal_index','relative_index','contig_gene_count','neighbors_up','neighbors_down'])
    write_csv_rows(out_sum, sum_rows, ['family_id','occurrences','genomes_with_family','median_relative_position','orientation_consistency','conserved_position_fraction','neighbor_conserved_fraction','neighbor_jaccard_mean'])

    # print unmapped summary
    if unmapped:
        sys.stderr.write('Unmapped members summary (per-family sample)\n')
        for f, ms in unmapped.items():
            sys.stderr.write(f + '\t' + ','.join(ms[:10]) + ('...' if len(ms)>10 else '') + '\n')
        # aggregate unmapped tokens and print top occurrences
        from collections import Counter
        all_un = []
        for ms in unmapped.values():
            all_un.extend(ms)
        c = Counter(all_un)
        sys.stderr.write('\nTop 50 unmapped tokens overall:\n')
        for token, cnt in c.most_common(50):
            sys.stderr.write(f'{token}\t{cnt}\n')


if __name__ == '__main__':
    run_from_cli()

