import gzip

import pytest

from PyamilySeq.order_analysis import parse_assembly_report


HEADER = (
    "# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\t"
    "GenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name\n"
)


def _row(name, role, molecule, location, gb, ref, unit="Primary Assembly", length="100"):
    return f"{name}\t{role}\t{molecule}\t{location}\t{gb}\t=\t{ref}\t{unit}\t{length}\t-\n"


def test_parser_indexes_accessions_and_classifies_assembly_report_records(tmp_path):
    path = tmp_path / "assembly_report.txt"
    path.write_text(
        "# report comment\n" + HEADER
        + _row("chr1", "assembled-molecule", "1", "Chromosome", "CM000001.1", "NC_000001.1")
        + _row("chr1_part2", "assembled-molecule", "1", "Chromosome", "CM000002.1", "NC_000002.1")
        + _row("pl1", "assembled-molecule", "plasmid pA", "Plasmid", "CP000001.1", "NC_100001.1")
        + _row("pl2", "assembled-molecule", "plasmid pB", "Plasmid", "CP000002.1", "NC_100002.1")
        + _row("unloc", "unlocalized-scaffold", "1", "Chromosome", "", "NW_000001.1", unit="-", length="50")
        + _row("unplaced", "unplaced-scaffold", "-", "-", "NZ_000001.1", "", unit="-", length="40")
    )
    records = parse_assembly_report(path)
    assert records["CM000001.1"]["Sequence-Name"] == "chr1"
    assert records["NC_000002.1"]["Assigned-Molecule"] == "1"
    assert records["CM000001.1"]["replicon_type"] == "chromosome"
    assert records["CP000001.1"]["replicon_type"] == records["CP000002.1"]["replicon_type"] == "plasmid"
    assert records["NW_000001.1"]["classification_status"] == "unlocalized_scaffold"
    assert records["NZ_000001.1"]["classification_status"] == "unplaced_scaffold"
    assert records["NW_000001.1"]["Assembly-Unit"] == "-"
    assert records["NZ_000001.1"]["replicon_topology"] == "unknown"


def test_parser_preserves_accessions_exactly_and_does_not_fuzzy_match_sequence_name(tmp_path):
    path = tmp_path / "assembly_report.txt"
    path.write_text(HEADER + _row("sequenceA", "assembled-molecule", "1", "Chromosome",
                                  "CM000010.1", "NC_000010.1"))
    records = parse_assembly_report(path)
    assert "CM000010.1" in records
    assert "NC_000010.1" in records
    assert "sequenceA" not in records
    assert "CM000010" not in records


def test_parser_supports_gzip_and_missing_optional_values(tmp_path):
    path = tmp_path / "assembly_report.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(HEADER + _row("scaf", "unplaced-scaffold", "-", "-", "-", "NW_000010.1",
                                   unit="-", length="17"))
    record = parse_assembly_report(path)["NW_000010.1"]
    assert record["GenBank-Accn"] == "-"
    assert record["Assigned-Molecule"] == "-"
    assert record["Sequence-Length"] == "17"


def test_parser_rejects_duplicate_and_contradictory_accessions(tmp_path):
    path = tmp_path / "assembly_report.txt"
    path.write_text(HEADER +
                    _row("one", "assembled-molecule", "1", "Chromosome", "CM000020.1", "NC_000020.1") +
                    _row("two", "unplaced-scaffold", "-", "-", "CM000020.1", "NC_000021.1"))
    with pytest.raises(ValueError, match="Contradictory assembly-report record.*CM000020.1"):
        parse_assembly_report(path)


@pytest.mark.parametrize("content,match", [
    ("", "no '# Sequence-Name' column header"),
    ("# report only\n", "no '# Sequence-Name' column header"),
    (HEADER + "short\trow\n", "has 2 fields; expected 10"),
    (HEADER + _row("x", "assembled-molecule", "1", "Chromosome", "CM1", "NC1", length="x"),
     "invalid Sequence-Length"),
    (HEADER + _row("x", "assembled-molecule", "1", "Chromosome", "CM1", "NC1") * 2,
     "Duplicate assembly-report record"),
])
def test_parser_rejects_invalid_or_empty_reports(tmp_path, content, match):
    path = tmp_path / "assembly_report.txt"
    path.write_text(content)
    with pytest.raises(ValueError, match=match):
        parse_assembly_report(path)
