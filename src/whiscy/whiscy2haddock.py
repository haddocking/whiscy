#!/usr/bin/env python3

"""Calculates the active/passive residues for HADDOCK

- Active residues: >=0.18 (WHISCY score)
- Passive residues: radius of 6.5A of active residues 
"""

__version__ = 1.0

import argparse
import os
import sys

import Bio.PDB

STANDARD_TYPES = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                  'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                  'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                  'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}


def parse_whiscy_scores(file_name):
    """Parses a WHISCY scores output file"""
    scores = {}
    with open(file_name) as input:
        for line in input:
            line = line.rstrip(os.linesep)
            fields = line.split()
            res = fields[1]
            score = float(fields[0])
            scores[res] = score
    return scores


if __name__ == "__main__":

    # Parse command line
    parser = argparse.ArgumentParser(prog='whiscy2haddock')
    parser.add_argument("input_pdb_file", help="Input PDB structure file", metavar="input_pdb_file")
    parser.add_argument("output_pdb_file", help="Output PDB file", metavar="output_pdb_file")
    parser.add_argument("whiscy_scores_file", help="WHISCY scores file", metavar="whiscy_scores_file")
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__))
    args = parser.parse_args()

    input_pdb_file_name = args.input_pdb_file
    output_pdb_file_name = args.output_pdb_file
    whiscy_scores_file_name = args.whiscy_scores_file

    # Parse scores from the WHISCY output file
    scores = parse_whiscy_scores(whiscy_scores_file_name)

    # Parse input PDB file
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(input_pdb_file_name, input_pdb_file_name)
    residues = structure.get_residues()

    # Calculate active residues
    active_residues = []
    for residue in residues:
        try:
            res_id = "{}{}".format(STANDARD_TYPES[residue.get_resname()], residue.get_id()[1])
            if scores[res_id] >= 0.18:
                active_residues.append(residue)
        except KeyError:
            pass

    # Calculate passive residues
    atoms  = Bio.PDB.Selection.unfold_entities(structure, 'A')
    neighbor_search = Bio.PDB.NeighborSearch(atoms)
    passive_residue_ids = []
    for active in active_residues:
        center = active['CA'].coord
        neighbors = neighbor_search.search(center, 6.5, level='R')
        # Remove itself
        neighbors.remove(active)
        for passive in neighbors:
            res_id = "{}{}".format(STANDARD_TYPES[passive.get_resname()], passive.get_id()[1])
            passive_residue_ids.append(res_id)

    # Create unique ids for passive residues
    passive_residue_ids = list(set(passive_residue_ids))

    # Create same ID list for active:
    active_residue_ids = []
    for residue in active_residues:
        res_id = "{}{}".format(STANDARD_TYPES[residue.get_resname()], residue.get_id()[1])
        active_residue_ids.append(res_id)
    active_residue_ids = list(set(active_residue_ids))

    # Finally, remove from passive list the active:
    passive_residue_ids = list(set(passive_residue_ids) - set(active_residue_ids))

    with open(input_pdb_file_name, 'r') as input:
        with open(output_pdb_file_name, 'w') as output:
            for line in input:
                line = line.rstrip(os.linesep)
                if line and line.startswith("ATOM  "):
                    res = STANDARD_TYPES[line[17:20].strip()]
                    res_num = line[22:26].strip()
                    res_id = "{}{}".format(res, res_num)
                    try:
                        if res_id in active_residue_ids:
                            bfactor = 100.0
                        elif res_id in passive_residue_ids:
                            bfactor = -100.0
                        else:
                            bfactor = 0.0
                    except KeyError:
                        bfactor = 0.0
                    line = line[:61] + "%6.2f" % bfactor + line[66:]
                output.write(line + os.linesep)

    with open('active.txt', 'w') as output:
        output.write(','.join([str(x) for x in sorted([int(id[1:]) for id in active_residue_ids])]))

    with open('passive.txt', 'w') as output:
        output.write(','.join([str(x) for x in sorted([int(id[1:]) for id in passive_residue_ids])]))
