"""Util functions involving structure and PDB files"""

import os
import re

from Bio import AlignIO
from Bio.PDB.PDBList import PDBList
from Bio.PDB.PDBIO import Select
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.PDBData import protein_letters_3to1


class NotAlternative(Select):
    """Removes alternative AAs"""

    def accept_residue(self, residue):
        return is_aa(residue) and residue.id[2] == " "


_hydrogen = re.compile("[123 ]*H.*")


def is_hydrogen(atom):
    """Checks if atom is an hydrogen"""
    name = atom.get_id()
    return _hydrogen.match(name)


def download_pdb_structure(pdb_code, pdb_file_name, file_path="."):
    """Downloads a PDB structure from the Protein Data Bank"""
    pdbl = PDBList()
    file_name = pdbl.retrieve_pdb_file(
        pdb_code, file_format="pdb", pdir=file_path, overwrite=True
    )
    if os.path.exists(file_name):
        os.rename(file_name, pdb_file_name)
    else:
        raise Exception("Can not download structure: {0}".format(pdb_code))


def get_pdb_sequence(
    input_pdb_file: str,
    chain_id: str,
    mapping_output: bool = False,
    with_gaps: bool = False,
) -> str | dict:
    """Gets the PDB sequence in a dictionary"""
    mapping = {}
    pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = pdb_parser.get_structure(input_pdb_file, input_pdb_file)
    model = structure[0]
    chain = model[chain_id]
    residues = list(chain)
    for res in residues:
        # Remove alternative location residues
        if "CA" in res.child_dict and is_aa(res) and res.id[2] == " ":
            try:
                mapping[res.id[1]] = protein_letters_3to1[res.get_resname()]
            except KeyError:
                # Ignore non standard residues such as HIC, MSE, etc.
                pass

    if with_gaps:
        # Add missing gap residues by their residue number
        res_numbers = sorted(mapping.keys())
        start, end = res_numbers[0], res_numbers[-1]
        missing = sorted(set(range(start, end + 1)).difference(res_numbers))
        for m in missing:
            mapping[m] = "-"

    if mapping_output:
        return mapping
    else:
        return "".join([mapping[k] for k in sorted(mapping.keys())])


def map_protein_to_sequence_alignment(
    pdb_file, chain_id, sequence, phylip_file, output_file_name
):
    """Creates a dictionary .conv file mapping protein residue numeration to aligment"""
    mapping = get_pdb_sequence(pdb_file, chain_id, mapping_output=True)
    # Check if sequence is the same
    pdb_seq = "".join([mapping[k] for k in sorted(mapping.keys())])  # type: ignore
    if pdb_seq != sequence:
        raise SystemExit("ERROR: PDB sequence doest not match sequence alignment")

    # Account for gaps in phylipseq file
    # alignment = list(AlignIO.parse(phylip_file, format='phylip-sequential'))[0]
    # master_phylip = alignment[0].seq
    # if str(master_phylip.ungap('-')) != pdb_seq:
    #    raise SystemExit("ERROR: PDB sequence doest not match sequence alignment in phylip file")

    with open(output_file_name, "w") as output_handle:
        output_handle.write(
            "# Conversion table from {} and chain {} to sequence{}".format(
                pdb_file, chain_id, os.linesep
            )
        )
        seq_res_id = 1
        for pdb_res_id in sorted(mapping.keys()):  # type: ignore
            # Do not map gaps if any
            if mapping[pdb_res_id] != "-":
                output_handle.write(
                    "{0}     {1}{2}".format(pdb_res_id, seq_res_id, os.linesep)
                )
            seq_res_id += 1
