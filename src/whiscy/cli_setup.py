import argparse
import os
import shutil
import subprocess
import sys
import warnings

from Bio import AlignIO, BiopythonWarning, SeqIO
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser

warnings.simplefilter("ignore", BiopythonWarning)

# Logging
import logging

from whiscy.modules import (
    CUTOFF,
    HSSPCONV_BIN,
    PROTDIST_BIN,
    access,
    alignment,
    homologs,
    hssp,
    pdbutil,
)

logger = logging.getLogger("whiscy_setup")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.INFO)
formatter = logging.Formatter(
    "%(asctime)s %(module)s:%(lineno)d %(levelname)s - %(message)s"
)
ch.setFormatter(formatter)
logger.addHandler(ch)


def write_homologs(master_sequence, rows, output_file):
    """Write the structure sequence and its aligned homologs as FASTA."""
    with open(output_file, "w") as handle:
        handle.write(">MASTER{0}{1}{0}".format(os.linesep, master_sequence))
        for name, row in rows:
            handle.write(">{0}{1}{2}{1}".format(name, os.linesep, row))


def hsspconv(hssp_file, converted_hssp_file):
    """Converts a HSSP in version 3 to original format"""
    cmd = f"{HSSPCONV_BIN} < {hssp_file} > {converted_hssp_file}"
    try:
        subprocess.run(cmd, shell=True)
    except:
        subprocess.check_call(cmd, shell=True)


def calculate_protdist(phylip_file, protdist_output_file):
    """Calculates the protdist of the given MSA"""
    cmd = f"{PROTDIST_BIN} {phylip_file} {protdist_output_file} > /dev/null 2>&1"
    try:
        subprocess.run(cmd, shell=True)
    except:
        subprocess.check_call(cmd, shell=True)


def write_to_fasta(output_fasta_file, sequence):
    """Writes a sequence to a FASTA format file"""
    with open(output_fasta_file, "w") as output_handle:
        output_handle.write(">{0}{1}".format(output_fasta_file, os.linesep))
        n = 60
        seq = [sequence[i : i + n] for i in range(0, len(sequence), n)]
        for chunk in seq:
            output_handle.write("{0}{1}".format(chunk, os.linesep))


def main():

    # Parse command line
    parser = argparse.ArgumentParser(prog="whiscy_setup")
    parser.add_argument(
        "pdb_file_name",
        help="PDB file name (.pdb extension) or PDB code",
        metavar="pdb_file_name",
    )
    parser.add_argument("chain_id", help="Chain ID to be predicted", metavar="chain_id")
    parser.add_argument(
        "--hssp", help="HSSP ID code", dest="hssp", metavar="hssp", default=None
    )
    parser.add_argument(
        "--alignment",
        help="Alignment file",
        dest="alignment",
        metavar="alignment",
        default=None,
    )
    parser.add_argument(
        "--alignment_format",
        help="Alignment format",
        dest="alignment_format",
        metavar="alignment_format",
        default="fasta",
    )

    args = parser.parse_args()

    filename, file_extension = os.path.splitext(os.path.basename(args.pdb_file_name))
    with_pdb_code = file_extension == ""
    input_pdb_file = None

    if with_pdb_code:
        # PDB code has been specified instead of a PDB file name
        pdb_code = args.pdb_file_name
        input_pdb_file = "{0}.pdb".format(pdb_code)
        if not os.path.exists(input_pdb_file):
            try:
                pdbutil.download_pdb_structure(pdb_code, input_pdb_file)
            except Exception as err:
                logger.error(str(err))
                sys.exit(1)
        else:
            logger.warning(
                "PDB structure already exists ({0}), no need to download it again".format(
                    input_pdb_file
                )
            )
    else:
        pdb_code = filename
        input_pdb_file = args.pdb_file_name

    if not os.path.exists(input_pdb_file):
        logger.error("PDB structure file {0} not found".format(input_pdb_file))
        sys.exit(1)

    # Check if chain belongs to this PDB
    pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = pdb_parser.get_structure(filename, input_pdb_file)
    chain_ids = [chain.id for chain in structure.get_chains()]
    chain_id = args.chain_id.upper()
    if len(chain_id) > 1:
        logger.error("Wrong chain id {0}".format(chain_id))
        sys.exit(1)
    if chain_id not in chain_ids:
        logger.error(
            "Chain {0} provided not in available chains: {1}".format(
                chain_id, str(chain_ids)
            )
        )
        sys.exit(1)

    # Save only the given chain and discard residues with alternative positions
    io = PDBIO()
    current_pdb_file = "{0}_{1}.pdb".format(pdb_code, chain_id)
    for chain in structure.get_chains():
        if chain.id == chain_id:
            io.set_structure(chain)
            io.save(current_pdb_file, select=pdbutil.NotAlternative())
    logger.info(
        "PDB structure with chain {0} saved to {1}".format(chain_id, current_pdb_file)
    )

    # Calculate SASA:
    rsa_output_file = "{0}_{1}.rsa".format(pdb_code, chain_id)
    access.calculate_accessibility(current_pdb_file, rsa_output_file)
    logger.info("Atom accessibility calculated to {0}".format(rsa_output_file))

    # Calculate the different accessibility files according to the cutoffs:
    access.create_cutoff_files(rsa_output_file, pdb_code, chain_id, CUTOFF)
    logger.info("Surface and buried residues calculated")

    # Get structure sequence
    master_sequence = pdbutil.get_pdb_sequence(input_pdb_file, chain_id)
    write_to_fasta("{0}_{1}.fasta".format(filename, chain_id), master_sequence)

    if args.hssp:
        hssp_file = "{0}.hssp".format(args.hssp)
        phylip_file = "{0}_{1}.phylseq".format(args.hssp, chain_id)
    else:
        hssp_file = "{0}.hssp".format(pdb_code)
        phylip_file = "{0}_{1}.phylseq".format(pdb_code, chain_id)

    if with_pdb_code or args.hssp:
        # Get the HSSP alignment from FTP if pdb code specified
        if not os.path.exists(hssp_file):
            logger.info("Downloading HSSP alignment...")
            try:
                compressed_hssp_file = hssp.get_from_ftp(pdb_code)
                if not compressed_hssp_file:
                    compressed_hssp_file = hssp.get_from_url(pdb_code)

                if "hssp3" in compressed_hssp_file:
                    hssp_file = hssp_file.replace("hssp", "hssp3")
                hssp.decompress_bz2(compressed_hssp_file, hssp_file)
                logger.info("HSSP alignment stored to {0}".format(hssp_file))
            except Exception:
                logger.warning("HSSP file could not be downloaded")
        try:
            if "hssp3" in hssp_file:
                # HSSP downloaded file is in new HSSP3 format, need to be
                # converted back to original HSSP format using hsspconv
                converted_hssp_file = hssp_file.replace("hssp3", "hssp")
                hsspconv(hssp_file, converted_hssp_file)
                hssp_file = converted_hssp_file

            hssp.hssp_file_to_phylip(hssp_file, phylip_file, chain_id, master_sequence)
            logger.info("HSSP file converted to PHYLIP format")

        except Exception as e:
            logger.debug(e)

    if not os.path.exists(hssp_file):
        # if no alignment is provided, will do it automatically
        if not args.alignment:
            logger.info("HSSP file not found, searching for homologous sequences")
            # The homologs are aligned to the structure sequence, one column
            # per residue; a previous search result is reused if present
            homologs_file = "{0}_{1}_homologs.fasta".format(filename, chain_id)
            if os.path.exists(homologs_file):
                logger.info("Using the homologs in {0}".format(homologs_file))
                rows = [
                    (record.id, str(record.seq))
                    for record in SeqIO.parse(homologs_file, "fasta")
                ][1:]
            else:
                try:
                    rows = homologs.find_homologs(
                        str(master_sequence), os.environ.get("WHISCY_EBI_EMAIL")
                    )
                except homologs.HomologSearchError as err:
                    logger.error(str(err))
                    sys.exit(1)
                write_homologs(master_sequence, rows, homologs_file)
                logger.info("Homologs stored in {0}".format(homologs_file))

            output_phylseq_file = "{0}_{1}.phylseq".format(filename, chain_id)
            phylseq_rows = [("MASTER", str(master_sequence))]
            for n, (name, row) in enumerate(rows, 1):
                phylseq_rows.append(("{0}_{1}".format(n, name), row))
            alignment.write_phylseq(phylseq_rows, output_phylseq_file)
            logger.info(
                "{0} written with {1} homologs".format(output_phylseq_file, len(rows))
            )

        else:
            # Alignment is provided, need to convert it to phylipseq
            # which is the format WHISCY is expecting
            alignment_file_name = args.alignment
            alignment_format = args.alignment_format.lower()
            phylip_file = "{0}_{1}.phylseq".format(pdb_code, chain_id)
            if alignment_format == "phylip":
                # Nothing to convert, just rename file if needed
                if os.path.basename(alignment_file_name) != phylip_file:
                    shutil.copyfile(alignment_file_name, phylip_file)
            else:
                AlignIO.convert(
                    alignment_file_name,
                    alignment_format,
                    phylip_file,
                    "phylip-sequential",
                )

    if not os.path.exists(phylip_file):
        logger.error("PHYLIP sequence file {0} not found".format(phylip_file))
        sys.exit(1)

    # WHISCY scores alignment column i as residue i of the target: keep only
    # the columns where the first (target) sequence has a residue, and map the
    # PDB residues onto them
    try:
        rows = alignment.project_on_reference(alignment.read_phylseq(phylip_file))
        pdb_residues = pdbutil.get_pdb_sequence(
            current_pdb_file, chain_id, mapping_output=True
        )
        residue_mapping = alignment.map_pdb_to_reference(pdb_residues, rows[0][1])
    except alignment.AlignmentMappingError as err:
        logger.error(str(err))
        sys.exit(1)
    alignment.write_phylseq(rows, phylip_file)
    logger.info(
        "Alignment mapped: {0} of {1} PDB residues, {2} alignment columns".format(
            len(residue_mapping), len(pdb_residues), len(rows[0][1])
        )
    )

    # Calculate protdist
    protdist_output_file = "{0}_{1}.out".format(filename, chain_id)
    calculate_protdist(phylip_file, protdist_output_file)
    logger.info("Protdist calculated")

    # Generate conversion table file
    conv_output_file = "{0}_{1}.conv".format(filename, chain_id)
    alignment.write_conversion_table(
        residue_mapping, current_pdb_file, chain_id, conv_output_file
    )
    logger.info("Conversion table file generated")

    logger.info("Whiscy setup finished")


if __name__ == "__main__":
    main()
