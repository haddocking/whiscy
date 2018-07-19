#!/usr/bin/env python3

"""Adds the WHISCY scoring prediction to the B-factor column of a given PDB"""

import os
import sys
import argparse


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
    parser = argparse.ArgumentParser(prog='whiscy2bfactor')
    parser.add_argument("input_pdb_file", help="Input PDB structure file", metavar="input_pdb_file")
    parser.add_argument("output_pdb_file", help="Output PDB file", metavar="output_pdb_file")
    parser.add_argument("whiscy_scores_file", help="WHISCY scores file", metavar="whiscy_scores_file")
    args = parser.parse_args()

    input_pdb_file_name = args.input_pdb_file
    output_pdb_file_name = args.output_pdb_file
    whiscy_scores_file_name = args.whiscy_scores_file

    # Parse scores from the WHISCY output file
    scores = parse_whiscy_scores(whiscy_scores_file_name)

    with open(input_pdb_file_name, 'r') as input:
        with open(output_pdb_file_name, 'w') as output:
            for line in input:
                line = line.rstrip(os.linesep)
                if line and line.startswith("ATOM  "):
                    res = STANDARD_TYPES[line[17:20].strip()]
                    res_num = line[22:26].strip()
                    res_id = "{}{}".format(res, res_num)
                    try:
                        # Normalized Score = 100 * WHISCYSCORE + 50
                        # Values are truncated between 0 and 100
                        normalized_score = max(0.0, min(100.0, 100.0 * (scores[res_id] + 50.0)))
                    except KeyError:
                        normalized_score = 0.0

                    line = line[:61] + "%6.2f" % normalized_score + line[66:]
                output.write(line + os.linesep)

    print("{} PDB file with WHISCY scores in B-factor column has been created".format(output_pdb_file_name))
