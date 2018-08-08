
"""Util functions involving structure and PDB files"""

import os
import re
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa, three_to_one
from Bio.PDB.PDBIO import Select


class NotAlternative(Select):
    """Removes alternative AAs"""
    def accept_residue(self, residue):
        return (is_aa(residue) and residue.id[2] == ' ')


_hydrogen = re.compile("[123 ]*H.*")
def is_hydrogen(atom):
    """Checks if atom is an hydrogen"""
    name = atom.get_id() 
    return _hydrogen.match(name)


def download_pdb_structure(pdb_code, pdb_file_name, file_path='.'):
    """Downloads a PDB structure from the Protein Data Bank"""
    pdbl = PDBList()
    file_name = pdbl.retrieve_pdb_file(pdb_code, file_format='pdb', pdir=file_path, overwrite=True)
    if os.path.exists(file_name):
        os.rename(file_name, pdb_file_name)
    else:
        logger.error("Can not download structure: {0}".format(pdb_code))
        raise SystemExit


def get_pdb_sequence(input_pdb_file, chain_id, mapping_output=False):
    """Gets the PDB sequence in a dictionary"""
    mapping = {}
    pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = pdb_parser.get_structure(input_pdb_file, input_pdb_file)
    model = structure[0]
    chain = model[chain_id]
    for res in chain:
        # Remove alternative location residues
        if "CA" in res.child_dict and is_aa(res) and res.id[2] == ' ':
            mapping[res.id[1]] = three_to_one(res.get_resname())
    if mapping_output:
        return mapping
    else:
        return ''.join([mapping[k] for k in sorted(mapping.keys())])


def map_protein_to_sequence_alignment(pdb_file, chain_id, sequence, output_file_name):
    """Creates a dictionary .conv file mapping protein residue numeration to aligment"""
    mapping = get_pdb_sequence(pdb_file, chain_id, mapping_output=True)
    # Check if sequence is the same
    pdb_seq = ''.join([mapping[k] for k in sorted(mapping.keys())])
    if pdb_seq != sequence:
        raise SystemExit("ERROR: PDB sequence doest not match sequence alignment")

    with open(output_file_name, 'w') as output_handle:
        for seq_res_id, pdb_res_id in enumerate(sorted(mapping.keys())):
            output_handle.write("{0}     {1}{2}".format(pdb_res_id, seq_res_id+1, os.linesep))
