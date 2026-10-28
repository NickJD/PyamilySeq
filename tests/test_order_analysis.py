from pathlib import Path

import pytest

from PyamilySeq.order_analysis import (
    Contig,
    Feature,
    Genome,
    apply_fasta_lengths,
    infer_contig_orientations,
    load_family_map,
    load_fasta_lengths,
    parse_gff,
)


def test_repeated_gene_names_do_not_collapse_distinct_loci(tmp_path):
    gff = tmp_path / "sample.gff3"
    gff.write_text(
        "ctg\tRef\tCDS\t10\t20\t.\t+\t0\tID=cds1;Name=dup\n"
        "ctg\tRef\tCDS\t40\t50\t.\t-\t0\tID=cds2;Name=dup\n"
    )
    genome = parse_gff(gff)
    assert len(genome.genes) == 2
    assert [gene.feature_id for gene in genome.genes] == ["cds1", "cds2"]


def test_parent_attributes_resolve_through_transcript_chain(tmp_path):
    gff = tmp_path / "sample.gff3"
    gff.write_text(
        "ctg\tRef\tgene\t1\t100\t.\t+\t.\tID=gene1;locus_tag=LT1;product=parent product\n"
        "ctg\tRef\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "ctg\tRef\tCDS\t10\t90\t.\t+\t0\tID=cds1;Parent=tx1\n"
    )
    gene = parse_gff(gff).genes[0]
    assert gene.locus_tag == "LT1"
    assert gene.product == "parent product"
    assert gene.attributes["ID"] == "cds1"


def test_invalid_sequence_region_is_warned_and_falls_back_to_feature_end(tmp_path):
    gff = tmp_path / "sample.gff3"
    gff.write_text(
        "##sequence-region ctg 100 1\n"
        "ctg\tRef\tCDS\t10\t90\t.\t+\t0\tID=cds1\n"
    )
    genome = parse_gff(gff)
    assert genome.contigs["ctg"].length == 90
    assert genome.contigs["ctg"].length_source == "max_feature_end"
    assert any("invalid ##sequence-region" in warning for warning in genome.warnings)


def test_conflicting_scoped_family_assignments_fail(tmp_path):
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text("genome\tgene_id\tfamily_id\nA\tg1\tF1\nA\tg1\tF2\n")
    with pytest.raises(ValueError, match="Conflicting family assignments"):
        load_family_map(mapping)


def test_pyamilyseq_wide_gpa_csv_maps_space_separated_members(tmp_path):
    mapping = tmp_path / "gene_presence_absence.csv"
    mapping.write_text(
        "Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,"
        "Genome Fragment,Order within Fragment,Accessory Fragment,Accessory Order with Fragment,QC,"
        "Min group size nuc,Max group size nuc,Avg group size nuc,A,B\n"
        "family_1,,,2,3,1.5,,,,,,100,100,100,A_gene1  A_gene2,B_gene1\n"
    )
    family_map = load_family_map(mapping)
    assert family_map.lookup("A", ["A_gene1"]) == "family_1"
    assert family_map.lookup("A", ["A_gene2"]) == "family_1"
    assert family_map.lookup("B", ["B_gene1"]) == "family_1"


def test_ambiguous_fasta_match_is_not_selected(tmp_path):
    gff = tmp_path / "g.gff3"
    gff.write_text("ctg\tRef\tCDS\t1\t9\t.\t+\t0\tID=g1\n")
    fasta_dir = tmp_path / "fastas"
    fasta_dir.mkdir()
    (fasta_dir / "g.fa").write_text(">ctg\nAAAAAAAAAA\n")
    (fasta_dir / "g.fna").write_text(">ctg\nAAAAAAAAAAAA\n")
    genome = parse_gff(gff)
    warnings = apply_fasta_lengths([genome], fasta_dir)
    assert genome.contigs["ctg"].length == 9
    assert any("Multiple FASTA files match" in warning for warning in warnings)


def test_fasta_contigs_without_gff_features_are_included_in_assembly_stats(tmp_path):
    gff = tmp_path / "g.gff3"
    gff.write_text("ctg1\tRef\tCDS\t1\t9\t.\t+\t0\tID=g1\n")
    fasta_dir = tmp_path / "fastas"
    fasta_dir.mkdir()
    (fasta_dir / "g.fa").write_text(">ctg1\nAAAAAAAAAA\n>unannotated\nACGTACGTACGT\n")
    genome = parse_gff(gff)
    warnings = apply_fasta_lengths([genome], fasta_dir)
    assert genome.contigs["ctg1"].length == 10
    assert genome.contigs["unannotated"].length == 12
    assert genome.contigs["unannotated"].length_source == "fasta_no_gff_features"
    assert "included in assembly contig counts" in warnings[0]


def test_fasta_duplicate_identifiers_are_rejected(tmp_path):
    gff = tmp_path / "g.gff3"
    gff.write_text("ctg\tRef\tCDS\t1\t9\t.\t+\t0\tID=g1\n")
    fasta_dir = tmp_path / "fastas"
    fasta_dir.mkdir()
    (fasta_dir / "g.fa").write_text(">ctg description one\nAAAA\n>ctg description two\nTTTT\n")
    with pytest.raises(ValueError, match="Duplicate FASTA contig identifier"):
        load_fasta_lengths(parse_gff(gff), fasta_dir)


def _genome(name, contig_orders):
    genome = Genome(name=name, path=Path(f"{name}.gff"))
    for contig_name, families in contig_orders.items():
        contig = Contig(name=contig_name, length=len(families) * 100)
        for index, family in enumerate(families):
            feature = Feature(name, contig_name, index * 100 + 1, index * 100 + 90,
                              "+", "CDS", f"{name}_{family}", family=family)
            contig.genes.append(feature)
            genome.genes.append(feature)
        genome.contigs[contig_name] = contig
    return genome


def test_orientation_handles_reverse_order_and_reference_contig_ties():
    reference = _genome("ref", {"r1": ["A", "B", "C", "D", "E"], "r2": ["V", "W", "X", "Y", "Z"]})
    query = _genome("query", {"q1": ["E", "D", "C", "B", "A"]})
    _, rows = infer_contig_orientations([reference, query], "ref", min_anchors=5)
    q1 = next(row for row in rows if row["genome"] == "query")
    assert q1["orientation"] == "reverse"
    assert q1["reference_contig"] == "r1"

    tied_reference = _genome("ref2", {"r1": ["A", "B", "C", "D", "E"], "r2": ["A", "B", "C", "D", "E"]})
    tied_query = _genome("query2", {"q1": ["A", "B", "C", "D", "E"]})
    _, rows = infer_contig_orientations([tied_reference, tied_query], "ref2", min_anchors=5)
    q1 = next(row for row in rows if row["genome"] == "query2")
    assert q1["orientation"] == "undetermined"
    assert q1["orientation_status"] == "ambiguous_reference_contig"
    assert q1["reference_contig"] == ""


def test_orientation_tau_threshold_is_configurable():
    reference = _genome("ref", {"r1": ["A", "B", "C", "D", "E"]})
    query = _genome("query", {"q1": ["A", "C", "B", "E", "D"]})
    _, rows = infer_contig_orientations([reference, query], "ref", min_anchors=5, min_abs_tau=0.9)
    q1 = next(row for row in rows if row["genome"] == "query")
    assert q1["orientation"] == "undetermined"
    assert q1["orientation_status"] == "insufficient_or_ambiguous_anchors"
