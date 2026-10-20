"""Util functions involving structure and PDB files"""

import os
import re

from Bio.PDB import PDBList
from Bio.PDB.PDBIO import Select
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa, three_to_one


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


def get_pdb_sequence(input_pdb_file, chain_id, mapping_output=False, with_gaps=False):
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
                mapping[res.id[1]] = three_to_one(res.get_resname())
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
