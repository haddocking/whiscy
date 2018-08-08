#!/usr/bin/env python3

"""This script calculates the minimum distance between all residues in the PDB structure.

It assumes that the PDB file contains only one chain and there are no residues with
alternative positions. As this script should be executed after whiscy_setup, it is safe
to make such assumptions.
"""

import os
import argparse
import re
from Bio.PDB.PDBParser import PDBParser
# Logging
import logging
logging.basicConfig(format='%(name)s [%(levelname)s] %(message)s', level=logging.INFO)
logger = logging.getLogger("residue_distance")


def load_conversion_table(file_name):
    """Loads the data from a conversion dictionary file.

    Tipically, this file has extension .conv and its content is two columns.
    Left column is the id of the position of the residue in the original PDB 
    structure. Right Column is the id of the position of the residue in Phylseq
    aligment file.
    """
    conversion_table = {}
    with open(file_name, "rU") as handle:
        for line in handle:
            if line:
                fields = line.rstrip(os.linesep).split()
                try:
                    from_index, to_index = int(fields[0]), int(fields[1])
                    conversion_table[from_index] = to_index
                except ValueError:
                    pass
    return conversion_table


_hydrogen = re.compile("[123 ]*H.*")
def is_hydrogen(atom):
    """Checks if atom is an hydrogen"""
    name = atom.get_id() 
    return _hydrogen.match(name)


if __name__ == "__main__":

    # Parse command line
    parser = argparse.ArgumentParser(prog='residue_distance')
    parser.add_argument("pdb_file", help="PDB file", metavar="pdb_file")
    parser.add_argument("conv_file", help="Conversion table file", metavar="conv_file")
    parser.add_argument("output_file", help="Output file name", metavar="output_file")
    args = parser.parse_args()

    # Read conversion table file
    logger.info("Reading conversion table")
    conversion_table = load_conversion_table(args.conv_file)

    # Parse PDB structure
    logger.info("Reading PDB structure from {}".format(args.pdb_file))
    pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = pdb_parser.get_structure(args.pdb_file, args.pdb_file)
    model = structure[0]
    residues = [residue for residue in model.get_residues()]

    # Calculate minimum distance between residues
    keys = conversion_table.keys()
    with open(args.output_file, 'w') as output_handle:
        for i in range(len(residues)):
            res1 = residues[i]
            if res1.id[1] in keys:
                for j in range(i+1, len(residues)):
                    min_distance = float("inf")
                    res2 = residues[j]
                    if res2.id[1] in keys:
                        for a1 in res1:
                            if not is_hydrogen(a1):
                                for a2 in res2:
                                    if not is_hydrogen(a2):
                                        d = a1 - a2
                                        if d < min_distance:
                                            min_distance = d
                        output_handle.write("{0} {1} {2:.6f}{3}".format(res1.id[1], res2.id[1], 
                                                                        min_distance, os.linesep))
    
    logger.info("Residue distances written to {}".format(args.output_file))
