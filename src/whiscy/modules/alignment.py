"""Map a multiple sequence alignment onto the residues of a PDB chain.

WHISCY scores alignment column ``i`` as residue ``i`` of the target protein,
so the alignment written to the .phylseq file must have exactly one column per
target residue, and the .conv file must map each PDB residue number to its
column. This module builds both from any alignment whose reference row is the
target protein (the uploaded alignment's first sequence, or the PDB sequence
added to a BLAST/MUSCLE alignment).
"""

import os

from Bio import Align
from Bio.Align import substitution_matrices

# Minimum fraction of the PDB residues that must be aligned to an identical
# residue of the reference sequence
MIN_IDENTITY = 0.9


class AlignmentMappingError(Exception):
    pass


def project_on_reference(rows, reference_index=0):
    """Keep only the alignment columns where the reference row has a residue.

    ``rows`` is a list of (name, aligned_sequence) tuples. Returns the rows
    with gap columns of the reference removed, so the reference row becomes
    its ungapped sequence and column ``i`` is residue ``i`` of the reference.
    """
    lengths = {len(seq) for _, seq in rows}
    if len(lengths) != 1:
        raise AlignmentMappingError(
            "Sequences in the alignment do not all have the same length"
        )
    reference = rows[reference_index][1]
    keep = [i for i, c in enumerate(reference) if c not in "-."]
    return [(name, "".join(seq[i] for i in keep)) for name, seq in rows]


def map_pdb_to_reference(pdb_residues, reference, min_identity=MIN_IDENTITY):
    """Map PDB residue numbers to 1-based positions in the reference sequence.

    ``pdb_residues`` is a dict {residue number: one-letter code}. The PDB
    sequence is aligned to the (ungapped) reference, which may be longer,
    e.g. the full-length protein when the structure covers a domain only.
    Raises AlignmentMappingError if the PDB chain does not match the
    reference closely enough.
    """
    numbers = sorted(pdb_residues)
    pdb_seq = "".join(pdb_residues[n] for n in numbers)
    reference = reference.upper()

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    # The reference may extend beyond the PDB chain on either side
    aligner.target_end_gap_score = 0
    aligner.query_end_gap_score = 0
    alignment = aligner.align(reference, pdb_seq)[0]

    mapping = {}
    identical = 0
    for (ref_start, ref_end), (pdb_start, pdb_end) in zip(*alignment.aligned):
        for ref_pos, pdb_pos in zip(range(ref_start, ref_end), range(pdb_start, pdb_end)):
            mapping[numbers[pdb_pos]] = ref_pos + 1
            if reference[ref_pos] == pdb_seq[pdb_pos]:
                identical += 1

    identity = identical / len(pdb_seq)
    if identity < min_identity:
        raise AlignmentMappingError(
            "The first sequence of the alignment does not match the PDB chain "
            f"({identity:.0%} of the PDB residues identical, at least "
            f"{min_identity:.0%} required). The first sequence must be the "
            "sequence of the protein in the PDB file"
        )
    return mapping


def read_phylseq(file_name):
    """Read a WHISCY .phylseq file: header, then 10-character name + sequence."""
    with open(file_name) as handle:
        lines = [line.rstrip("\r\n") for line in handle]
    return [(line[:10], line[10:]) for line in lines[1:] if line.strip()]


def write_phylseq(rows, file_name):
    with open(file_name, "w") as handle:
        handle.write(f"{len(rows)}  {len(rows[0][1])}{os.linesep}")
        for name, seq in rows:
            handle.write(f"{name[:10]:10s}{seq}{os.linesep}")


def write_conversion_table(mapping, pdb_file, chain_id, file_name):
    with open(file_name, "w") as handle:
        handle.write(
            f"# Conversion table from {pdb_file} and chain {chain_id} to sequence{os.linesep}"
        )
        for pdb_res_id in sorted(mapping):
            handle.write(f"{pdb_res_id}     {mapping[pdb_res_id]}{os.linesep}")
