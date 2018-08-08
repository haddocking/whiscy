
"""Whiscy Data files"""

import os


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
