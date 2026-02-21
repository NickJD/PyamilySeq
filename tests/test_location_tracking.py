from PyamilySeq.location_tracking import parse_gff, parse_clusters, map_family_members_to_features, compute_neighbors_for_contigs, compute_family_metrics


def test_location_tracking_small(tmp_path):
    # create a tiny GFF
    gff_text = """##gff-version 3
contig1	Ref	CDS	100	200	.	+	0	ID=geneA
contig1	Ref	CDS	300	400	.	+	0	ID=geneB
contig1	Ref	CDS	500	600	.	+	0	ID=geneC
contig2	Ref	CDS	50	150	.	-	0	ID=geneX
contig2	Ref	CDS	200	300	.	-	0	ID=geneY
"""
    gff_file = tmp_path / "test.gff"
    gff_file.write_text(gff_text)

    # simple cluster file: family -> members (members use gene IDs)
    cluster_text = """fam1\tgeneA,geneX\nfam2\tgeneB,geneY\n"""
    cluster_file = tmp_path / "clusters.txt"
    cluster_file.write_text(cluster_text)

    contigs, id_map = parse_gff(str(gff_file))
    families = parse_clusters(str(cluster_file))
    fam_occ, unmapped = map_family_members_to_features(families, id_map)
    compute_neighbors_for_contigs(contigs, k=1)
    occ_rows, sum_rows = compute_family_metrics(fam_occ, k=1, rel_delta=0.5)

    # expectations: fam1 has 2 occurrences, fam2 has 2
    fams = {r['family_id']: r for r in sum_rows}
    assert 'fam1' in fams and fams['fam1']['occurrences'] == 2
    assert 'fam2' in fams and fams['fam2']['occurrences'] == 2

    # occurrences rows should include geneA and geneX
    found_ids = set([r['gene_id'] for r in occ_rows])
    assert 'geneA' in found_ids and 'geneX' in found_ids
