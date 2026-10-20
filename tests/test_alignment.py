import pytest

from whiscy.modules.alignment import (
    AlignmentMappingError,
    map_pdb_to_reference,
    project_on_reference,
    read_phylseq,
    write_phylseq,
)

# 1PPE chain I
PDB_SEQ = "RVCPRILMECKKDSDCLAECVCLEHGYCG"
PDB_RESIDUES = {n: aa for n, aa in zip(range(1, len(PDB_SEQ) + 1), PDB_SEQ)}


def test_project_removes_reference_gap_columns():
    rows = [
        ("target", "-RV-CP"),
        ("hit1", "MRVACP"),
        ("hit2", "-KVA-P"),
    ]
    assert project_on_reference(rows) == [
        ("target", "RVCP"),
        ("hit1", "RVCP"),
        ("hit2", "KV-P"),
    ]


def test_project_rejects_unequal_lengths():
    with pytest.raises(AlignmentMappingError):
        project_on_reference([("a", "RVCP"), ("b", "RVC")])


def test_map_identical_sequence():
    mapping = map_pdb_to_reference(PDB_RESIDUES, PDB_SEQ)
    assert mapping == {n: n for n in PDB_RESIDUES}


def test_map_pdb_numbering_offset():
    # PDB numbering starting at 33, as in an N-terminally truncated construct
    residues = {n + 32: aa for n, aa in PDB_RESIDUES.items()}
    mapping = map_pdb_to_reference(residues, PDB_SEQ)
    assert mapping[33] == 1
    assert mapping[33 + len(PDB_SEQ) - 1] == len(PDB_SEQ)


def test_map_structure_covers_part_of_reference():
    # The reference is the full-length protein, the structure a fragment
    reference = "MSTNKQLAGHE" + PDB_SEQ + "WWKKPLLE"
    mapping = map_pdb_to_reference(PDB_RESIDUES, reference)
    assert mapping[1] == 12
    assert mapping[len(PDB_SEQ)] == 11 + len(PDB_SEQ)


def test_map_pdb_residue_numbering_gap():
    # Residues 10-12 are missing in the structure
    residues = {n: aa for n, aa in PDB_RESIDUES.items() if n not in (10, 11, 12)}
    mapping = map_pdb_to_reference(residues, PDB_SEQ)
    assert mapping == {n: n for n in residues}


def test_map_point_mutation_is_accepted():
    mutated = dict(PDB_RESIDUES)
    mutated[5] = "A"
    mapping = map_pdb_to_reference(mutated, PDB_SEQ)
    assert mapping[5] == 5


def test_map_unrelated_reference_is_rejected():
    with pytest.raises(AlignmentMappingError, match="does not match the PDB chain"):
        map_pdb_to_reference(PDB_RESIDUES, "MKTAYIAKQRQISFVKSHFSRQLEERLGL")


def test_gapped_reference_end_to_end():
    """Map the uploaded-alignment case that crashed consadjust (KeyError: '-')."""
    # After projection every PDB residue maps to a column where the reference
    # has a residue, never a gap
    reference = "MSTN--KQLAGHE" + PDB_SEQ[:10] + "---" + PDB_SEQ[10:] + "WW-KK"
    rows = [("target", reference), ("hit", reference.replace("-", "A"))]
    projected = project_on_reference(rows)
    ref_seq = projected[0][1]
    mapping = map_pdb_to_reference(PDB_RESIDUES, ref_seq)
    for pdb_res, col in mapping.items():
        assert ref_seq[col - 1] == PDB_RESIDUES[pdb_res]


def test_phylseq_round_trip(tmp_path):
    rows = [("MASTER", "RVCP"), ("gb|QAV5757", "RV-P")]
    path = tmp_path / "x.phylseq"
    write_phylseq(rows, path)
    assert path.read_text().splitlines()[0] == "2  4"
    assert [(n.strip(), s) for n, s in read_phylseq(path)] == rows
