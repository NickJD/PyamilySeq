from types import SimpleNamespace

from PyamilySeq.Group_Compare import compare_cd_hit_clusters
from PyamilySeq.Group_Splitter import calculate_new_rep_seq
from PyamilySeq.Group_Summary import categorise_percentage
from PyamilySeq.Seq_Finder import find_ids_in_csv
from PyamilySeq.utils import detect_sequence_type, feature_matches_gene_ident, read_fasta, translate_frame


def test_group_splitter_representative_minimizes_length_and_identity_deviation():
    entries = [
        {"header": "best", "length": 100, "percent_identity": 90},
        {"header": "worst", "length": 100, "percent_identity": 100},
        {"header": "other", "length": 100, "percent_identity": 80},
    ]
    assert calculate_new_rep_seq(entries)["header"] == "best"


def test_group_summary_categories_cover_zero_and_one_hundred_percent():
    assert categorise_percentage(0) == "0-20%"
    assert categorise_percentage(100) == "99-100%"


def test_group_compare_ignores_cluster_number_reordering(tmp_path):
    first = tmp_path / "first.clstr"
    second = tmp_path / "second.clstr"
    output = tmp_path / "comparison.tsv"
    first.write_text(
        ">Cluster 0\n0\t3nt, >a... *\n1\t3nt, >b... at 100.00%\n"
        ">Cluster 1\n0\t3nt, >c... *\n"
    )
    second.write_text(
        ">Cluster 0\n0\t3nt, >c... *\n"
        ">Cluster 1\n0\t3nt, >a... *\n1\t3nt, >b... at 100.00%\n"
    )
    compare_cd_hit_clusters(str(first), str(second), str(output))
    assert output.read_text().splitlines() == [
        "Sequence\tCluster ID (File 1)\tCluster ID (File 2)\tChange Type"
    ]


def test_sequence_utilities_use_canonical_ids_and_preserve_stops(tmp_path):
    fasta = tmp_path / "descriptions.fa"
    fasta.write_text(">G1|geneA description one\natgaaatag\n")
    assert read_fasta(str(fasta)) == {"G1|geneA": "atgaaatag"}
    assert translate_frame("atgaaatag") == "MK*"

    protein = tmp_path / "protein.fa"
    protein.write_text(">protein\nMKTWYV\n")
    assert detect_sequence_type(str(protein)) is False
    assert feature_matches_gene_ident("pseudogene", "gene") is False


def test_seq_finder_reports_identifier_then_matching_rows(tmp_path):
    matrix = tmp_path / "matrix.csv"
    matrix.write_text("family,G1,G2\nfam1,a,b\nfam2,c,a\n")
    found = find_ids_in_csv(["a", "c"], str(matrix))
    assert found == {"a": ["fam1", "fam2"], "c": ["fam2"]}
