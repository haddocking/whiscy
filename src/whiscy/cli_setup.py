#!/usr/bin/env python3

"""Whiscy predictor setup"""

__version__ = 1.0

import argparse
import json
import os
import shutil
import subprocess
import sys
import warnings

# Import SearchIO and suppress experimental warning
from Bio import AlignIO, BiopythonExperimentalWarning, BiopythonWarning, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from Bio.Blast import NCBIWWW
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser

warnings.simplefilter("ignore", BiopythonWarning)
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio import SearchIO

# Logging
import logging

from libwhiscy import (
    PROTDIST_BIN,
    access,
    hssp,
    pdbutil,
    CUTOFF,
    HSSPCONV_BIN,
    MUSCLE_BIN,
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


def muscle_msa(
    input_sequence_file: str, output_alignment_file: str
) -> MultipleSeqAlignment:
    """Calculates a MSA using MUSCLE's Biopython wrapper"""

    assert MUSCLE_BIN is not None, "MUSCLE_BIN not found in system variables"

    muscle_cline = MuscleCommandline(
        MUSCLE_BIN, input=input_sequence_file, out=output_alignment_file
    )
    _, _ = muscle_cline()
    msa = AlignIO.read(output_alignment_file, "fasta")
    return msa


def ncbi_blast(fasta_file: str, output_file: str) -> None:
    """Performs a remote BLAST against the NCBI server"""
    record = SeqIO.read(fasta_file, format="fasta")
    result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"))
    with open(output_file, "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()


def msa_to_phylseq(
    msa: MultipleSeqAlignment, master_sequence: str, output_file: str
) -> None:
    """Converts a MSA to a Phylip Seq file"""
    with open(output_file, "w") as output_handle:
        # Write header
        output_handle.write(f"{len(msa)}  {len(master_sequence)}{os.linesep}")
        # Write master sequence
        output_handle.write(f"MASTER    {master_sequence}{os.linesep}")
        # Write the rest of alignments
        for alignment in msa:
            output_handle.write(f"{alignment.id[:10]:10s}{alignment.seq}{os.linesep}")


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
    parser.add_argument(
        "--version", action="version", version="%(prog)s {}".format(__version__)
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
                raise SystemExit
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
        raise SystemExit

    # Check if chain belongs to this PDB
    pdb_parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = pdb_parser.get_structure(filename, input_pdb_file)
    chain_ids = [chain.id for chain in structure.get_chains()]
    chain_id = args.chain_id.upper()
    if len(chain_id) > 1:
        logger.error("Wrong chain id {0}".format(chain_id))
        raise SystemExit
    if chain_id not in chain_ids:
        logger.error(
            "Chain {0} provided not in available chains: {1}".format(
                chain_id, str(chain_ids)
            )
        )
        raise SystemExit

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
            except Exception as err:
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
        logger.info("HSSP file not found, fallback to generating MSA with blastp")
        # if no alignment is provided, will do it automatically
        if not args.alignment:
            # Run BLASTP if needed
            blast_output_file = "{0}_{1}_blast.xml".format(filename, chain_id)
            input_sequence_file = "{0}_{1}.fasta".format(filename, chain_id)
            if not os.path.exists(blast_output_file):
                logger.info("Running blastp via NCBI")
                ncbi_blast(input_sequence_file, blast_output_file)
                logger.info("Result stored in {0}".format(blast_output_file))

            assert os.path.exists(blast_output_file), "BLAST file not found"

            # Convert file to FASTA format
            blast_qresult = SearchIO.read(blast_output_file, "blast-xml")
            records = []
            for hit in blast_qresult:
                records.append(hit[0].hit)
            blast_fasta_file = "{0}_{1}_blast.fasta".format(filename, chain_id)
            SeqIO.write(records, blast_fasta_file, "fasta")

            # Multiple sequence alignment
            logger.info("Generating MSA using MUSCLE...")
            output_alignment_file = "{0}_{1}_msa.fasta".format(filename, chain_id)
            msa = muscle_msa(blast_fasta_file, output_alignment_file)
            logger.info("Done.")

            # Convert MSA to Phylipseq
            logger.info("Converting MSA file to Phylseq format...")
            output_phylseq_file = "{0}_{1}.phylseq".format(filename, chain_id)
            msa_to_phylseq(
                msa=msa,
                master_sequence=str(master_sequence),
                output_file=output_phylseq_file,
            )

            assert os.path.exists(output_phylseq_file), "PHYLIP file not found"
            logger.info(f"{output_phylseq_file} file written")

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
        raise SystemExit

    # Calculate protdist
    protdist_output_file = "{0}_{1}.out".format(filename, chain_id)
    calculate_protdist(phylip_file, protdist_output_file)
    logger.info("Protdist calculated")

    # Generate conversion table file
    conv_output_file = "{0}_{1}.conv".format(filename, chain_id)
    pdbutil.map_protein_to_sequence_alignment(
        current_pdb_file, chain_id, master_sequence, phylip_file, conv_output_file
    )
    logger.info("Conversion table file generated")

    logger.info("Whiscy setup finished")


if __name__ == "__main__":
    main()
