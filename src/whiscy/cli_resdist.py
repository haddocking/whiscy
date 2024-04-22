#!/usr/bin/env python3

"""This script calculates the minimum distance between all residues in the PDB structure.

It assumes that the PDB file contains only one chain and there are no residues with
alternative positions. As this script should be executed after whiscy_setup, it is safe
to make such assumptions.
"""

__version__ = 1.0

import argparse

# Logging
import logging
import os
import sys

from Bio.PDB.PDBParser import PDBParser

from whiscy.modules.pdbutil import is_hydrogen
from whiscy.modules.whiscy_data import load_conversion_table

logger = logging.getLogger("residue_distance")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(name)s [%(levelname)s] %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)


def main():

    # Parse command line
    parser = argparse.ArgumentParser(prog="residue_distance")
    parser.add_argument("pdb_file", help="PDB file", metavar="pdb_file")
    parser.add_argument("conv_file", help="Conversion table file", metavar="conv_file")
    parser.add_argument("output_file", help="Output file name", metavar="output_file")
    parser.add_argument(
        "--version", action="version", version="%(prog)s {}".format(__version__)
    )
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
    with open(args.output_file, "w") as output_handle:
        for i in range(len(residues)):
            res1 = residues[i]
            if res1.id[1] in keys:
                for j in range(i + 1, len(residues)):
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
                        output_handle.write(
                            "{0} {1} {2:.6f}{3}".format(
                                res1.id[1], res2.id[1], min_distance, os.linesep
                            )
                        )

    logger.info("Residue distances written to {}".format(args.output_file))


if __name__ == "__main__":
    main()
