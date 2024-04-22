import argparse
import logging
import math
import os
import sys
from pathlib import Path

import numpy as np

from whiscy.modules import PARAM_PATH
from whiscy.modules.whiscy_data import (
    load_cons_file,
    load_residue_weights,
    load_z_table,
)

logger = logging.getLogger("consadjust")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(name)s [%(levelname)s] %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)

RESIDUE_WEIGHT_FILE = Path(PARAM_PATH, "weight_ma.txt")
ZTABLE_FILE = Path(PARAM_PATH, "ztable.txt")


def main():

    # Parse command line
    parser = argparse.ArgumentParser(prog="consadjust")
    parser.add_argument("cons_file", help="Conservation file", metavar="cons_file")

    parser.add_argument(
        "-o",
        "--output",
        help="If set, output prediction to this file",
        dest="output_file",
        metavar="output_file",
    )
    args = parser.parse_args()

    if not os.path.exists(args.cons_file):
        logger.error("Conservation file {} does not exist".format(args.cons_file))
        raise SystemExit

    logger.info("Reading input files")

    resweight = load_residue_weights(RESIDUE_WEIGHT_FILE)

    residues = load_cons_file(args.cons_file)

    # Statistics
    sum_scores = 0.0
    square_sum_scores = 0.0
    for residue in residues:
        sum_scores += residue.score
        square_sum_scores += residue.score * residue.score
    consnr = len(residues)
    stddev = math.sqrt((square_sum_scores - sum_scores * sum_scores / consnr) / consnr)
    mean = sum_scores / consnr

    if stddev <= 0.0:
        logger.error(
            "All identical values in conservation file {}".format(args.cons_file)
        )
        raise SystemExit

    zcalc = [(res.score - mean) / stddev for res in residues]

    # Load Z-matrix
    z_values = np.array(load_z_table(ZTABLE_FILE))
    if len(z_values) != 25000:
        logger.error("Reading error in Z-table {}".format(args.cons_file))
        raise SystemExit

    for n in range(consnr):
        pscore = 0.0
        neg = False
        currz = zcalc[n]
        if currz < 0:
            currz *= -1
            neg = True

        index = np.argmax(z_values < currz)
        z_values[index]

        if index >= 24999:
            pscore = 0.5
        if index <= 0:
            pscore = 0.0
        else:
            pscore = index - (currz - z_values[index]) / (
                z_values[index - 1] - z_values[index]
            )

        newneg = False
        pscore0 = pscore
        if neg:
            pscore0 = 50000 - pscore

        padj0 = 1000 * pscore0
        if resweight[residues[n].code]:
            padj0 = pscore0 / resweight[residues[n].code]

        padj = padj0
        if padj > 25000:
            padj = 50000 - padj
            newneg = True

        if padj < 0:
            padj = 0

        pmin = int(math.floor(padj))
        pmax = int(math.ceil(padj))
        fac = padj - pmin
        znorm = fac * z_values[pmin] + (1 - fac) * z_values[pmax]

        residues[n].score = mean + znorm * (1.0 - 2.0 * newneg) * stddev

    logger.info("Subtracting average value...")
    scoresum = 0.0
    for n in range(consnr):
        scoresum += residues[n].score
    scoreaverage = scoresum / consnr
    for n in range(consnr):
        residues[n].score -= scoreaverage

    logger.info("Sorting scores...")
    residues_sorted = sorted(residues, key=lambda res: res.score, reverse=True)

    logger.info("Writing scores...")
    # If output to file
    if args.output_file:
        with open(args.output_file, "w") as output_handle:
            for res in residues_sorted:
                output_handle.write(
                    "{0:7.5f}  {1}{2}{3}".format(
                        res.score, res.code, res.nr, os.linesep
                    )
                )

    else:
        for res in residues_sorted:
            print("{0:7.5f}  {1}{2}".format(res.score, res.code, res.nr))


if __name__ == "__main__":
    main()
