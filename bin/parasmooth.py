#!/usr/bin/env python3

import math
import os
import argparse
# Logging
import logging
logging.basicConfig(format='%(name)s [%(levelname)s] %(message)s', level=logging.INFO)
logger = logging.getLogger("parasmooth")


class Residue:
    """Represents a residue"""
    def __init__(self, nr, code, score):
        self.nr = nr
        self.code = code
        self.score = score

    def __str__(self):
        return "{}.{}: {}".format(self.code, self.nr, self.score)


class Distance:
    """Represents the distance between two residues"""
    def __init__(self, nr1, nr2, dis):
        self.nr1 = nr1
        self.nr2 = nr2
        self.dis = dis

    def __str__(self):
        return "{} - {}: {}".format(self.nr1, self.nr2, self.dis)


def get_weight(par, val):
    """Calculates the weight in the par dictionary of the key val"""
    i = 0
    for k in sorted(par.keys()):
        if k > val:
            vh, wh = k, par[k]
            break
        i += 1

    if math.isclose(val, vh) or i==0:
        return wh 
    else:
        i -= 1
        key = list(par.keys())[i]
        vl = key
        wl = par[key]

        return wl + (wh - wl) * (val - vl)/(vh - vl)


def read_surface_cons_file(file_name):
    """Reads and parses a .acons file"""
    residues = []
    with open(file_name, 'rU') as handle:
        for line in handle:
            if line:
                fields = line.rstrip(os.linesep).split()
                try:
                    score = float(fields[0])
                    code = fields[1][0].upper()
                    nr = int(fields[1][1:])
                    residues.append(Residue(nr, code, score))
                except:
                    logger.error("Reading error in surface conservation file {}".format(file_name))
                    raise SystemExit
    return residues


def read_low_accessible_cons_file(file_name):
    """Reads and parses a .lcons file"""
    residues = []
    with open(file_name, 'rU') as handle:
        for line in handle:
            if line:
                fields = line.rstrip(os.linesep).split()
                try:
                    score = float(fields[0])
                    code = fields[1][0].upper()
                    nr = int(fields[1][1:])
                    residues.append(Residue(nr, code, score))
                except:
                    logger.error("Reading error in low-accessible conservation file {}".format(file_name))
                    raise SystemExit
    return residues


def read_residue_distance_matrix(file_name):
    """Reads and parses the residue distance matrix"""
    distances = []
    with open(file_name, 'rU') as handle:
        for line in handle:
            if line:
                fields = line.rstrip(os.linesep).split()
                try:
                    nr1 = int(fields[0])
                    nr2 = int(fields[1])
                    distance = float(fields[2])
                    distances.append(Distance(nr1, nr2, distance))
                except:
                    logger.error("Reading error in distance matrix file {}".format(file_name))
                    raise SystemExit
    return distances


def read_smoothing_parameter_file(file_name):
    """Reads and parses the smoothing parameter file"""
    par = {}
    with open(file_name, 'rU') as handle:
        for line in handle:
            if line:
                fields = line.rstrip(os.linesep).split()
                try:
                    v = float(fields[0])
                    w = float(fields[1])
                    par[v] = w
                except:
                    logger.error("Reading error in smoothing parameter file {}".format(file_name))
                    raise SystemExit
    return par


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


if __name__ == "__main__":

    # Parse command line
    parser = argparse.ArgumentParser(prog='parasmooth')
    parser.add_argument("surface_cons_file", help="Surface conservation file", metavar="surface_cons_file")
    parser.add_argument("low_accessible_cons_file", help="Low accessible conservation file", metavar="low_accessible_cons_file")
    parser.add_argument("residue_distance_matrix", help="Residue distance matrix", metavar="residue_distance_matrix")
    parser.add_argument("smoothing_parameter_file", help="Smoothing parameter file", metavar="smoothing_parameter_file")
    parser.add_argument("-o", "--output", help="If sets, output prediction to this file", 
                        dest="output_file", metavar="output_file")
    args = parser.parse_args()

    logger.info("Reading input files")

    # Read .acons file
    res_sur = read_surface_cons_file(args.surface_cons_file)

    # Read .lcons file
    res_lac = read_low_accessible_cons_file(args.low_accessible_cons_file)

    # Read .rd file
    resdist = read_residue_distance_matrix(args.residue_distance_matrix)

    # Read .par file
    par = read_smoothing_parameter_file(args.smoothing_parameter_file)

    # Calculate parasmooth values
    logger.info("Calculating parameter smoothing")
    res_sur = calculate_parasmooth(res_sur, res_lac, resdist, par)

    # If output to file
    if args.output_file:
        with open(args.output_file, 'w') as handle:
            for n in range(len(res_sur)):
                handle.write("{:8.5f}   {}{}{}".format(res_sur[n].score, 
                                                       res_sur[n].code, 
                                                       res_sur[n].nr,
                                                       os.linesep))
        logger.info("Result written to {}".format(args.output_file))
    else:
        for n in range(len(res_sur)):
            print("{:8.5f}   {}{}".format(res_sur[n].score, res_sur[n].code, res_sur[n].nr))
