#!/usr/bin/env python3

"""Python 3 implementation of the Whiscy C binary"""

import argparse
import os
from libwhiscy.pam_calc import pam_load_sequences, pam_calc_similarity
from libwhiscy.pam_data import code
from libwhiscy.quotes import get_one


def load_surface_list(file_name):
    """Loads the data from a surface list file.

    Tipically, this file has extension .sur and its content is a 
    simple column containing the index of the residue which has been
    marked as surface.
    """
    surface_list = []
    with open(file_name, "rU") as input:
        for line in input:
            if line:
                try:
                    surface_list.append(int(line))
                except ValueError:
                    pass
    return surface_list


def load_conversion_table(file_name):
    """Loads the data from a conversion dictionary file.

    Tipically, this file has extension .conv and its content is two columns.
    Left column is the id of the position of the residue in the original PDB 
    structure. Right Column is the id of the position of the residue in Phylseq
    aligment file.
    """
    conversion_table = {}
    with open(file_name, "rU") as input:
        for line in input:
            if line:
                fields = line.rstrip(os.linesep).split()
                try:
                    from_index, to_index = int(fields[0]), int(fields[1])
                    conversion_table[from_index] = to_index
                except ValueError:
                    pass
    return conversion_table


class Residue():
    """Represents a residue number and its score"""
    def __init__(self, nr, score):
        self.nr = nr
        self.score = score



if __name__ == "__main__":

    # Parse command line
    parser = argparse.ArgumentParser(prog='whiscy')
    parser.add_argument("surface_list", help="Surface list", metavar="surface_list")
    parser.add_argument("conversion_table", help="Conversion table", metavar="conversion_table")
    parser.add_argument("alignment_file", help="Alignment file", metavar="alignment_file")
    parser.add_argument("distance_file", help="Distance file", metavar="distance_file")
    parser.add_argument("-o", "--output", help="If sets, output prediction to this file", 
                        dest="output_file", metavar="output_file")
    args = parser.parse_args()

    print("Parsing surface list...")
    surface_list = load_surface_list(args.surface_list)

    print("Loading conversion table...")
    conversion_table = load_conversion_table(args.conversion_table)

    print("Converting...")
    converted_surface = []
    for n in range(len(surface_list)):
        if not surface_list[n] in conversion_table or conversion_table[surface_list[n]] < 1: 
            print("WARNING: Surface residue number {0} cannot be converted", surface_list[n])
            print("Continuing program...")
        else:
            converted_surface.append(conversion_table[surface_list[n]])
    
    if not len(converted_surface):
        raise SystemExit("ERROR: No surface residues")

    print("Initializing score calculation...")
    seqnr, seqlen, refseq, seq_distances, sequences, seqtodis = pam_load_sequences(args.alignment_file, args.distance_file)

    print("Calculating scores...")
    realsum = 0
    totlist = []
    for n in range(len(converted_surface)):
        r = Residue(converted_surface[n], -1000.)
        totlist.append(r)
        if r.nr > seqlen or r.nr < 1:
            raise SystemExit("ERROR: surface residue out of range")
        posnr, distances, scores = pam_calc_similarity(converted_surface[n]-1, seqnr, sequences, seq_distances)
        if posnr <= 0 or posnr > seqnr:
            continue
        realsum += 1
        r.score = scores[posnr - 1]

    if realsum == 0:
        raise SystemExit("ERROR: No sequence information for any surface residues")

    print("Subtracting average value ...")

    scoresum = 0.0
    for n in range(realsum):
        scoresum += totlist[n].score
  
    scoreaverage = scoresum / realsum
    for n in range(realsum):
        totlist[n].score -= scoreaverage

    print("Sorting scores...")
    sorted_totlist = sorted(totlist, key=lambda res: res.score, reverse=True)

    print("Writing scores...")
    prediction = ""
    for n in range(realsum):
        r = sorted_totlist[n]
        rid = "{0}{1}".format(refseq[r.nr-1], r.nr)
        prediction += "{:7.5f}  {:^5s}{}".format(r.score, rid, os.linesep)

    if args.output_file:
        with open(args.output_file, "w") as output:
            output.write(prediction)
        print("Prediction written to {0}".format(args.output_file))
    else:
        print(prediction)

    # Ending with a quote
    print(get_one())
