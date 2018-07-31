#!/usr/bin/env python3

import math
import os
from collections import OrderedDict
import argparse


class Residue:
    def __init__(self, nr, code, score):
        self.nr = nr
        self.code = code
        self.score = score

    def __str__(self):
        return "{}.{}: {}".format(self.code, self.nr, self.score)


class Distance:
    def __init__(self, nr1, nr2, dis):
        self.nr1 = nr1
        self.nr2 = nr2
        self.dis = dis

    def __str__(self):
        return "{} - {}: {}".format(self.nr1, self.nr2, self.dis)


def get_weight(par, val):
    """Calculates the weight in the OrderedDict par of the key val"""
    i = 0
    for k, v in par.items():
        if v <= val:
            vh, wh = k, v
        else:
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
                    raise SystemExit("ERROR: Reading error in surface conservation file {}".format(file_name))
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
                    raise SystemExit("ERROR: Reading error in low-accessible conservation file {}".format(file_name))
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
                    raise SystemExit("ERROR: Reading error in distance matrix file {}".format(file_name))
    return distances


if __name__ == "__main__":

    # Parse command line
    parser = argparse.ArgumentParser(prog='param_smooth')
    parser.add_argument("surface_cons_file", help="Surface conservation file", metavar="surface_cons_file")
    parser.add_argument("low_accessible_cons_file", help="Low accessible conservation file", metavar="low_accessible_cons_file")
    parser.add_argument("residue_distance_matrix", help="Residue distance matrix", metavar="residue_distance_matrix")
    parser.add_argument("smoothing_parameter_file", help="Smoothing parameter file", metavar="smoothing_parameter_file")
    parser.add_argument("-o", "--output", help="If sets, output prediction to this file", 
                        dest="output_file", metavar="output_file")
    args = parser.parse_args()

    # Read .acons file:
    res_sur = read_surface_cons_file(args.surface_cons_file)

    # Read .lcons file:
    res_lac = read_low_accessible_cons_file(args.low_accessible_cons_file)

    # Read .rd file:
    resdist = read_residue_distance_matrix(args.residue_distance_matrix)
