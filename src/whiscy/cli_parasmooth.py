#!/usr/bin/env python3

__version__ = 1.0

import argparse

# Logging
import logging
import math
import os
import sys

logger = logging.getLogger("parasmooth")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(name)s [%(levelname)s] %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)
from libwhiscy.whiscy_data import (
    load_cons_file,
    load_residue_distance_matrix,
    load_smoothing_parameter_file,
)


def get_weight(par, val):
    """Calculates the weight in the par dictionary of the key val"""
    i = 0
    for k in sorted(par.keys()):
        if k > val:
            vh, wh = k, par[k]
            break
        i += 1

    if math.isclose(val, vh) or i == 0:
        return wh
    else:
        i -= 1
        key = list(par.keys())[i]
        vl = key
        wl = par[key]

        return wl + (wh - wl) * (val - vl) / (vh - vl)


def calculate_parasmooth(res_sur, res_lac, resdist, par):
    """Main function, calculates parasmooth values"""

    # Get the max distance
    maxdis = sorted(par.keys())[-1]

    for n in range(len(res_sur)):
        weight = 1.0
        score = res_sur[n].score
        for i in range(len(resdist)):
            if resdist[i].dis > maxdis:
                continue
            partner = -1

            if resdist[i].nr1 == res_sur[n].nr:
                partner = resdist[i].nr2
            elif resdist[i].nr2 == res_sur[n].nr:
                partner = resdist[i].nr1
            if partner == -1:
                continue

            found = False
            if resdist[i].dis < maxdis:
                for nn in range(len(res_sur)):
                    if n == nn:
                        continue
                    if res_sur[nn].nr == partner:
                        found = True
                        currweight = get_weight(par, resdist[i].dis)
                        currscore = currweight * res_sur[nn].score
                        weight += currweight
                        score += currscore
                        break
            if found:
                continue

            if resdist[i].dis < maxdis:
                for nn in range(len(res_lac)):
                    if res_lac[nn].nr == partner:
                        currweight = get_weight(par, resdist[i].dis)
                        currscore = currweight * res_lac[nn].score
                        weight += currweight
                        score += currscore
                        break

        res_sur[n].score = score / weight

    # Sort by scores in reverse order:
    res_sur = sorted(res_sur, key=lambda res: res.score, reverse=True)

    return res_sur


def main():

    # Parse command line
    parser = argparse.ArgumentParser(prog="parasmooth")
    parser.add_argument(
        "surface_cons_file",
        help="Surface conservation file",
        metavar="surface_cons_file",
    )
    parser.add_argument(
        "low_accessible_cons_file",
        help="Low accessible conservation file",
        metavar="low_accessible_cons_file",
    )
    parser.add_argument(
        "residue_distance_matrix",
        help="Residue distance matrix",
        metavar="residue_distance_matrix",
    )
    parser.add_argument(
        "smoothing_parameter_file",
        help="Smoothing parameter file",
        metavar="smoothing_parameter_file",
    )
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

    logger.info("Reading input files")

    try:
        # Read .acons file
        res_sur = load_cons_file(args.surface_cons_file)

        # Read .lcons file
        res_lac = load_cons_file(args.low_accessible_cons_file)

        # Read .rd file
        resdist = load_residue_distance_matrix(args.residue_distance_matrix)

        # Read .par file
        par = load_smoothing_parameter_file(args.smoothing_parameter_file)

        # Calculate parasmooth values
        logger.info("Calculating parameter smoothing")
        res_sur = calculate_parasmooth(res_sur, res_lac, resdist, par)
    except Exception as err:
        logger.error(str(err))
        raise SystemExit

    # If output to file
    if args.output_file:
        with open(args.output_file, "w") as handle:
            for n in range(len(res_sur)):
                handle.write(
                    "{:8.5f}   {}{}{}".format(
                        res_sur[n].score, res_sur[n].code, res_sur[n].nr, os.linesep
                    )
                )
        logger.info("Result written to {}".format(args.output_file))
    else:
        for n in range(len(res_sur)):
            print(
                "{:8.5f}   {}{}".format(
                    res_sur[n].score, res_sur[n].code, res_sur[n].nr
                )
            )


if __name__ == "__main__":
    main()
