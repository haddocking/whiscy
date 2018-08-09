
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
