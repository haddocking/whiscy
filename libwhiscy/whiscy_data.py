
"""Whiscy Data files"""

import os


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


_accepted_residue_codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
                           'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


def load_surface_list(file_name):
    """Loads the data from a surface list file.

    Tipically, this file has extension .sur and its content is a 
    simple column containing the index of the residue which has been
    marked as surface.
    """
    surface_list = []
    with open(file_name, "rU") as handle:
        for line in handle:
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


def load_surface_cons_file(file_name):
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
                    raise Exception("Reading error in surface conservation file {}".format(file_name))
    return residues


def load_low_accessible_cons_file(file_name):
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
                    raise Exception("Reading error in low-accessible conservation file {}".format(file_name))
    return residues


def load_residue_distance_matrix(file_name):
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
                    raise Exception("Reading error in distance matrix file {}".format(file_name))
    return distances


def load_smoothing_parameter_file(file_name):
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
                    raise Exception("Reading error in smoothing parameter file {}".format(file_name))
    return par


def load_residue_weights(file_name):
    """Loads a file containing for each line a residue in one letter code and a weight.

    Ignores lines starting with '#', only accepts a residue the first time it appears in the file
    and only if it belongs to standard list of residues.
    """
    resweight = {}
    with open(file_name, 'rU') as handle:
        for line in handle:
            line = line.rstrip(os.linesep)
            if line[0] != '#':
                try:
                    fields = line.split()
                    residue_id = fields[0].upper()
                    weight = float(fields[1])
                    if residue_id in _accepted_residue_codes and residue_id not in resweight:
                        resweight[residue_id] = weight
                except:
                    pass
    return resweight
    