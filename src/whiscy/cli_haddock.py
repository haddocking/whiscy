"""
Calculates the active/passive residues for HADDOCK

- Active residues: >=0.18 (WHISCY score), and the 10 best-scored residues
- Passive residues: surface residues within 6.5A of active residues. The
  surface is the WHISCY surface list (.sur) when given, otherwise the
  residues with a WHISCY score.

Writes active.txt, passive.txt and haddock-cutoff.pdb for the cutoff, and
active_top10.txt, passive_top10.txt and haddock-top10.pdb for the top 10.
"""

import argparse
import filecmp
import logging
import os
import shutil
import sys

import Bio.PDB

from whiscy import LOG_FILE, log_to_file

logger = logging.getLogger("whiscy_haddock")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(name)s [%(levelname)s] %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)

ACTIVE_CUTOFF = 0.18  # WHISCY score
PASSIVE_DISTANCE = 6.5  # A, from the CA of an active residue
TOP_N = 10  # active residues of the top-N selection

STANDARD_TYPES = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def parse_whiscy_scores(file_name):
    """Parses a WHISCY scores output file"""
    scores = {}
    with open(file_name) as input:
        for line in input:
            line = line.rstrip(os.linesep)
            fields = line.split()
            res = fields[1]
            score = float(fields[0])
            scores[res] = score
    return scores


def parse_surface_list(file_name):
    """Parses a WHISCY surface list (.sur): one residue number per line."""
    with open(file_name) as input:
        return {int(line.split()[0]) for line in input if line.strip()}


def residue_id(residue):
    """Return the WHISCY id of a residue (e.g. V75), or None if not standard."""
    code = STANDARD_TYPES.get(residue.get_resname())
    return code and "{}{}".format(code, residue.get_id()[1])


def haddock_residues(structure, active_residues, surface):
    """Return the ids of the active and passive residues for HADDOCK."""
    # Passive: surface residues with an atom within PASSIVE_DISTANCE of the CA
    # of an active residue, and not active themselves
    atoms = Bio.PDB.Selection.unfold_entities(structure, "A")
    neighbor_search = Bio.PDB.NeighborSearch(atoms)
    active_ids = {residue_id(r) for r in active_residues}
    passive_ids = set()
    for active in active_residues:
        center = active["CA"].coord
        for passive in neighbor_search.search(center, PASSIVE_DISTANCE, level="R"):
            if passive.get_id()[1] in surface:
                passive_ids.add(residue_id(passive))
    return active_ids, passive_ids - active_ids - {None}


def ids_to_numbers(residue_ids):
    """Return the sorted residue numbers of residue ids (e.g. V75 -> 75)."""
    return sorted(int(res_id[1:]) for res_id in residue_ids)


def write_residue_list(file_name, residue_ids):
    """Write the residue numbers comma-separated, as HADDOCK takes them."""
    with open(file_name, "w") as output:
        output.write(",".join(str(x) for x in ids_to_numbers(residue_ids)))


def write_haddock_pdb(input_pdb_file_name, output_pdb_file_name, active, passive):
    """Write the structure with B-factors 100 (active), -100 (passive) or 0."""
    with open(input_pdb_file_name, "r") as input:
        with open(output_pdb_file_name, "w") as output:
            for line in input:
                line = line.rstrip(os.linesep)
                if line and line.startswith("ATOM  "):
                    res = STANDARD_TYPES.get(line[17:20].strip())
                    res_id = "{}{}".format(res, line[22:26].strip())
                    if res_id in active:
                        bfactor = 100.0
                    elif res_id in passive:
                        bfactor = -100.0
                    else:
                        bfactor = 0.0
                    line = line[:61] + "%6.2f" % bfactor + line[66:]
                output.write(line + os.linesep)


def main():

    # Parse command line
    parser = argparse.ArgumentParser(prog="whiscy2haddock")
    parser.add_argument(
        "input_pdb_file", help="Input PDB structure file", metavar="input_pdb_file"
    )
    parser.add_argument(
        "whiscy_scores_file", help="WHISCY scores file", metavar="whiscy_scores_file"
    )
    parser.add_argument(
        "--surface_list",
        help="WHISCY surface list (.sur); passive residues are taken from it",
        default=None,
    )
    parser.add_argument(
        "--no-report",
        dest="report",
        action="store_false",
        help="Do not write report.html and state.mvsj (the prediction with a "
        "Mol* view)",
    )
    args = parser.parse_args()
    log_to_file(logger)

    input_pdb_file_name = args.input_pdb_file
    whiscy_scores_file_name = args.whiscy_scores_file

    # Parse scores from the WHISCY output file
    scores = parse_whiscy_scores(whiscy_scores_file_name)

    # Surface residues: only these can be passive
    if args.surface_list:
        surface = parse_surface_list(args.surface_list)
    else:
        surface = {int(res_id[1:]) for res_id in scores}

    # Parse input PDB file
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(input_pdb_file_name, input_pdb_file_name)
    scored = [r for r in structure.get_residues() if residue_id(r) in scores]

    # Active residues: by score cutoff, and the TOP_N best-scored residues
    active_residues = [r for r in scored if scores[residue_id(r)] >= ACTIVE_CUTOFF]
    top_residues = sorted(scored, key=lambda r: -scores[residue_id(r)])[:TOP_N]

    active_residue_ids, passive_residue_ids = haddock_residues(
        structure, active_residues, surface
    )
    top_active_ids, top_passive_ids = haddock_residues(structure, top_residues, surface)

    write_haddock_pdb(
        input_pdb_file_name,
        "haddock-cutoff.pdb",
        active_residue_ids,
        passive_residue_ids,
    )
    write_haddock_pdb(
        input_pdb_file_name,
        "haddock-top{}.pdb".format(TOP_N),
        top_active_ids,
        top_passive_ids,
    )
    write_residue_list("active.txt", active_residue_ids)
    write_residue_list("passive.txt", passive_residue_ids)
    write_residue_list("active_top{}.txt".format(TOP_N), top_active_ids)
    write_residue_list("passive_top{}.txt".format(TOP_N), top_passive_ids)

    if args.report:
        from whiscy.report import generate_report

        options = {
            "structure": input_pdb_file_name,
            "scores": whiscy_scores_file_name,
            "surface list": args.surface_list or "residues with a score",
            "active residues": "score >= {}".format(ACTIVE_CUTOFF),
            "top-{} active residues".format(TOP_N): "the {} best-scored".format(TOP_N),
            "passive residues": "surface, within {} A of an active CA".format(
                PASSIVE_DISTANCE
            ),
        }
        top = "top{}".format(TOP_N)
        files = [whiscy_scores_file_name]
        files += ["haddock-cutoff.pdb", "active.txt", "passive.txt"]
        files += ["haddock-{}.pdb".format(top)]
        files += ["active_{}.txt".format(top), "passive_{}.txt".format(top)]
        files += [LOG_FILE]
        # Link the input structure, copied next to the report if needed
        structure = os.path.basename(input_pdb_file_name)
        if not os.path.exists(structure):
            shutil.copyfile(input_pdb_file_name, structure)
        if filecmp.cmp(structure, input_pdb_file_name, shallow=False):
            files.insert(0, structure)
        selections = {
            "score \u2265 {}".format(ACTIVE_CUTOFF): (
                ids_to_numbers(active_residue_ids),
                ids_to_numbers(passive_residue_ids),
            ),
            "top {}".format(TOP_N): (
                ids_to_numbers(top_active_ids),
                ids_to_numbers(top_passive_ids),
            ),
        }
        fname = generate_report(
            "haddock-cutoff.pdb", files, scores, options, selections
        )
        logger.info("Report written to {}".format(os.path.relpath(fname)))
        logger.info(
            "To view it, run `python -m http.server` in this directory and open "
            "http://localhost:8000/report.html"
        )


if __name__ == "__main__":
    main()
