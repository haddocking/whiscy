#!/usr/bin/env python3

"""Python 3 implementation of the Whiscy C binary"""

__version__ = "1.0"

import argparse

# Logging
import logging
import os
import sys

from libwhiscy.pam_calc import pam_calc_similarity, pam_load_sequences
from libwhiscy.pam_data import code
from libwhiscy.quotes import get_one
from libwhiscy.whiscy_data import load_conversion_table, load_surface_list

logger = logging.getLogger("whiscy")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(name)s [%(levelname)s] %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)


class Residue:
    """Represents a residue number and its score"""

    def __init__(self, nr, score):
        self.nr = nr
        self.score = score


def main():

    # Parse command line
    parser = argparse.ArgumentParser(prog="whiscy")
    parser.add_argument("surface_list", help="Surface list", metavar="surface_list")
    parser.add_argument(
        "conversion_table", help="Conversion table", metavar="conversion_table"
    )
    parser.add_argument(
        "alignment_file", help="Alignment file", metavar="alignment_file"
    )
    parser.add_argument("distance_file", help="Distance file", metavar="distance_file")
    parser.add_argument(
        "-o",
        "--output",
        help="If set, output prediction to this file",
        dest="output_file",
        metavar="output_file",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s {}".format(__version__)
    )
    args = parser.parse_args()

    logger.info("Parsing surface list...")
    surface_list = load_surface_list(args.surface_list)

    logger.info("Loading conversion table...")
    conversion_table = load_conversion_table(args.conversion_table)

    logger.info("Converting...")
    converted_surface = []
    for n in range(len(surface_list)):
        if (
            not surface_list[n] in conversion_table
            or conversion_table[surface_list[n]] < 1
        ):
            logger.warning(
                "Surface residue number {0} cannot be converted".format(surface_list[n])
            )
            logger.info("Continuing program...")
        else:
            converted_surface.append(conversion_table[surface_list[n]])

    if not len(converted_surface):
        logger.error("No surface residues")
        raise SystemExit

    logger.info("Initializing score calculation...")
    try:
        seqnr, seqlen, refseq, seq_distances, sequences, seqtodis = pam_load_sequences(
            args.alignment_file, args.distance_file
        )
    except Exception as err:
        logger.error(str(err))
        raise SystemExit

    logger.info("Calculating scores...")
    realsum = 0
    totlist = []
    for n in range(len(converted_surface)):
        r = Residue(converted_surface[n], -1000.0)
        totlist.append(r)
        if r.nr > seqlen or r.nr < 1:
            logger.error("Surface residue out of range")
            raise SystemExit
        posnr, distances, scores = pam_calc_similarity(
            converted_surface[n] - 1, seqnr, sequences, seq_distances
        )
        if posnr <= 0 or posnr > seqnr:
            continue
        realsum += 1
        r.score = scores[posnr - 1]

    if realsum == 0:
        logger.error("No sequence information for any surface residues")
        raise SystemExit

    logger.info("Subtracting average value ...")

    scoresum = 0.0
    for n in range(realsum):
        scoresum += totlist[n].score

    scoreaverage = scoresum / realsum
    for n in range(realsum):
        totlist[n].score -= scoreaverage

    logger.info("Sorting scores...")
    sorted_totlist = sorted(totlist, key=lambda res: res.score, reverse=True)

    logger.info("Writing scores...")
    prediction = ""
    inv_conversion_table = {v: k for k, v in conversion_table.items()}
    for n in range(realsum):
        r = sorted_totlist[n]
        rid = "{0}{1}".format(refseq[r.nr - 1], inv_conversion_table[r.nr])
        prediction += "{:7.5f}  {:^5s}{}".format(r.score, rid, os.linesep)

    if args.output_file:
        with open(args.output_file, "w") as output:
            output.write(prediction)
        logger.info("Prediction written to {0}".format(args.output_file))
    else:
        print(prediction)

    # Ending with a quote
    print()
    print(get_one())


if __name__ == "__main__":
    main()
