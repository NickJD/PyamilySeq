from PyamilySeq.order_analysis import resolve_replicon_metadata


def ncbi(role="assembled-molecule", molecule="1", location="Chromosome", unit="Primary Assembly", topology=None):
    value = {"Sequence-Role": role, "Assigned-Molecule": molecule,
             "Assigned-Molecule-Location/Type": location, "Assembly-Unit": unit,
             "replicon_type": "unknown" if "unplaced" in role or "unlocalized" in role else
                              ("plasmid" if "plasmid" in location.casefold() else "chromosome"),
             "classification_status": "unplaced_scaffold" if "unplaced" in role else
                                      "unlocalized_scaffold" if "unlocalized" in role else "classified"}
    if topology is not None:
        value["replicon_topology"] = topology
    return value


def test_ncbi_and_gff_agree_on_chromosome_identity_and_preserve_evidence():
    n = ncbi()
    result = resolve_replicon_metadata("ctg1", n, gff_region={"replicon_type": "chromosome", "replicon_id": "Primary Assembly:chromosome:1"})
    assert result["replicon_type"] == "chromosome"
    assert result["replicon_id"] == "Primary Assembly:chromosome:1"
    assert result["conflict_status"] == "none"
    assert result["metadata_evidence"]["assembly_report"] == n
    assert set(result["metadata_sources"]) == {"assembly_report", "gff_region"}


def test_ncbi_and_gff_agree_on_plasmid_identity():
    n = ncbi(molecule="plasmid pA", location="Plasmid")
    out = resolve_replicon_metadata("pA_ctg", n, gff_region={"replicon_type": "plasmid", "replicon_id": "Primary Assembly:plasmid:plasmid pA"})
    assert out["replicon_type"] == "plasmid"
    assert out["replicon_id"] == "Primary Assembly:plasmid:plasmid pA"
    assert out["conflict_status"] == "none"


def test_chromosome_fragments_resolve_to_same_replicon_identity():
    a = resolve_replicon_metadata("frag1", ncbi())
    b = resolve_replicon_metadata("frag2", ncbi())
    assert a["replicon_id"] == b["replicon_id"]


def test_separate_plasmids_have_distinct_replicon_ids():
    a = resolve_replicon_metadata("p1", ncbi(molecule="plasmid pA", location="Plasmid"))
    b = resolve_replicon_metadata("p2", ncbi(molecule="plasmid pB", location="Plasmid"))
    assert a["replicon_id"] != b["replicon_id"]


def test_unlocalized_scaffold_retains_assigned_chromosome_identity():
    out = resolve_replicon_metadata("scaf", ncbi(role="unlocalized-scaffold"))
    assert out["replicon_type"] == "chromosome"
    assert out["replicon_id"] == "Primary Assembly:chromosome:1"
    assert out["sequence_role"] == "unlocalized-scaffold"


def test_unplaced_scaffold_remains_unassigned():
    out = resolve_replicon_metadata("scaf", ncbi(role="unplaced-scaffold", molecule="-", location="-"))
    assert out["replicon_id"] == "unknown"
    assert out["classification_status"] == "unplaced_scaffold"


def test_explicit_circular_topology_is_preserved_but_report_alone_is_unknown():
    assert resolve_replicon_metadata("c", ncbi(topology="circular"))["replicon_topology"] == "circular"
    assert resolve_replicon_metadata("c", ncbi())["replicon_topology"] == "unknown"


def test_missing_topology_and_unknown_values_do_not_conflict():
    out = resolve_replicon_metadata("c", ncbi(), fasta_header={"replicon_type": "unknown", "topology": "unknown"},
                                    gff_region={"Is_circular": "", "replicon_type": "unknown"})
    assert out["replicon_topology"] == "unknown"
    assert out["conflict_status"] == "none"


def test_conflicting_replicon_types_are_reported_and_resolve_unknown():
    out = resolve_replicon_metadata("c", ncbi(), fasta_header="c plasmid")
    assert out["replicon_type"] == "unknown"
    assert out["conflict_status"] == "type"
    assert "fasta_header" in out["metadata_sources"]


def test_conflicting_topologies_are_reported_and_resolve_unknown():
    out = resolve_replicon_metadata("c", ncbi(topology="circular"), gff_region={"Is_circular": "false"})
    assert out["replicon_topology"] == "unknown"
    assert out["conflict_status"] == "topology"


def test_missing_source_metadata_is_supported():
    out = resolve_replicon_metadata("unknown_contig")
    assert out["contig_id"] == "unknown_contig"
    assert out["replicon_id"] == out["replicon_type"] == out["replicon_topology"] == "unknown"
    assert out["classification_status"] == "unknown"
    assert out["metadata_evidence"] == {"assembly_report": {}, "fasta_header": None, "gff_region": {}}


def test_gff_feature_id_is_not_mistaken_for_shared_replicon_identity():
    out = resolve_replicon_metadata("c", gff_region={"ID": "unique-region-123", "gbkey": "Src"})
    assert out["replicon_id"] == "unknown"


def test_multiple_gff_regions_retain_evidence_and_conflicts():
    regions = [
        {"ID": "regionA", "replicon_type": "chromosome", "replicon_id": "chr1", "Is_circular": "true"},
        {"ID": "regionB", "replicon_type": "chromosome", "replicon_id": "chr1", "Is_circular": "true"},
    ]
    result = resolve_replicon_metadata("ctg", gff_region=regions)
    assert result["replicon_type"] == "chromosome"
    assert result["replicon_id"] == "chr1"
    assert result["metadata_evidence"]["gff_region"] == regions
    conflict = resolve_replicon_metadata("ctg", gff_region=regions + [
        {"replicon_type": "plasmid", "replicon_id": "p1", "Is_circular": "false"}])
    assert conflict["replicon_type"] == "unknown"
    assert conflict["replicon_topology"] == "unknown"
    assert "type" in conflict["conflict_status"]

