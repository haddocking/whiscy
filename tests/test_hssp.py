from pathlib import Path

import pytest

from whiscy.modules import hssp, pdbutil
from whiscy.modules.alignment import (
    map_pdb_to_reference,
    project_on_reference,
    read_phylseq,
)

from . import GOLDEN_DATA_PATH

HSSP_FILE = Path(GOLDEN_DATA_PATH, "setup", "1ppe.hssp")
PDB_FILE = Path(GOLDEN_DATA_PATH, "setup", "1ppe.pdb")


@pytest.mark.parametrize("chain_id", ["E", "I"])
def test_hssp_alignment_maps_onto_pdb_chain(tmp_path, chain_id):
    """Map the HSSP alignment onto the PDB chain, including insertion codes."""
    # Chain E (trypsin) uses insertion codes (e.g. 184A): those residues are in
    # the HSSP alignment but not in the PDB sequence WHISCY reads, so the
    # MASTER row must come from the HSSP file to match the alignment columns.
    pdb_residues = pdbutil.get_pdb_sequence(
        str(PDB_FILE), chain_id, mapping_output=True
    )
    master_sequence = "".join(pdb_residues[n] for n in sorted(pdb_residues))
    phylip_file = Path(tmp_path, f"1ppe_{chain_id}.phylseq")

    hssp.hssp_file_to_phylip(str(HSSP_FILE), phylip_file, chain_id, master_sequence)

    rows = read_phylseq(phylip_file)
    assert len({len(seq) for _, seq in rows}) == 1
    assert int(phylip_file.read_text().split()[0]) == len(rows)

    projected = project_on_reference(rows)
    mapping = map_pdb_to_reference(pdb_residues, projected[0][1], min_identity=1.0)
    assert len(mapping) == len(pdb_residues)
