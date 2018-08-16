#!/usr/bin/env python3

"""Whiscy predictor setup"""

__version__ = 1.0

import argparse
import os
import json
import subprocess
import warnings
import shutil
from Bio import AlignIO, SeqIO
from Bio.Blast import NCBIWWW
from Bio.Align.Applications import MuscleCommandline
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
# Import SearchIO and suppress experimental warning
from Bio import BiopythonExperimentalWarning, BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO
from libwhiscy import hssp
from libwhiscy import access
from libwhiscy import pdbutil

# Logging
import logging
logging.basicConfig(format='%(name)s [%(levelname)s] %(message)s', level=logging.INFO)
logger = logging.getLogger("whiscy_setup")

#Some extra checks...
if 'WHISCY_PATH' not in os.environ:
    logger.warning("WHISCY_PATH variable not defined in the environment. Assuming the path of the script...")
    os.environ['WHISCY_PATH'] = os.path.abspath(os.path.dirname(__file__))+os.path.sep
if shutil.which("freesasa") is None:
    logger.critical("FreeSASA not in PATH. Please set environment correctly.")
    raise SystemExit

def load_config(config_file='etc/local.json'):
    """Load Whiscy configuration"""
    with open(os.path.join(os.environ['WHISCY_PATH'], config_file), 'r') as f:
        config = json.load(f)
        return config


def muscle_msa(config, input_sequence_file, output_alignment_file):
    """Calculates a MSA using MUSCLE's Biopython wrapper"""
    muscle_bin = config['ALIGN']['MUSCLE_BIN']
    muscle_cline = MuscleCommandline(muscle_bin, input=input_sequence_file, out=output_alignment_file)
    if not os.path.exists(muscle_bin):
        logger.critical("The path defined for the MUSCLE binary is not correct. Check the configuration file!")
        raise SystemExit
    stdout, stderr = muscle_cline()
    MultipleSeqAlignment = AlignIO.read(output_alignment_file, "fasta") 
    return MultipleSeqAlignment


def ncbi_blast(fasta_file, output_file):
    """Performs a remote BLAST against the NCBI server"""
    record = SeqIO.read(fasta_file, format="fasta")
    result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"))
    with open(output_file, "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle.close()


def msa_to_phylseq(msa, master_sequence, output_file):
    """Converts a MSA to a Phylip Seq file"""
    with open(output_file, 'w') as output_handle:
        # Write header
        output_handle.write("{0}  {1}{2}".format(len(msa), 
                                                 len(master_sequence),
                                                 os.linesep))

        # Write master sequence
        output_handle.write("MASTER    {0}{1}".format(master_sequence,
                                                      os.linesep))
        # Write the rest of alignments
        for alignment in msa:
            output_handle.write("{:10s}{}{}".format(alignment.id[:10],
                                                    alignment.seq,
                                                    os.linesep))


def calculate_protdist(phylip_file, protdist_output_file):
    """Calculates the protdist of the given MSA"""
    protdist_bin = os.path.join(os.path.dirname(os.path.realpath(__file__)), "bin", "protdist", "protdist")
    if not os.path.exists(protdist_bin):
        logger.critical("Protdist was not compiled. Please check the installation instructions")
        raise SystemExit
    cmd = "{0} {1} {2} > /dev/null 2>&1".format(protdist_bin, phylip_file, protdist_output_file)
    try:
        subprocess.run(cmd, shell=True)
    except:
        subprocess.check_call(cmd, shell=True)


def write_to_fasta(output_fasta_file, sequence):
    """Writes a sequence to a FASTA format file"""
    with open(output_fasta_file, 'w') as output_handle:
        output_handle.write(">{0}{1}".format(output_fasta_file, os.linesep))
        n = 60
        seq = [sequence[i:i+n] for i in range(0, len(sequence), n)]
        for chunk in seq:
            output_handle.write("{0}{1}".format(chunk, os.linesep))


if __name__ == "__main__":

    # Parse command line
    parser = argparse.ArgumentParser(prog='whiscy_setup')
    parser.add_argument("pdb_file_name", help="PDB file name (.pdb extension) or PDB code", metavar="pdb_file_name")
    parser.add_argument("chain_id", help="Chain ID to be predicted", metavar="chain_id")
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__))
    args = parser.parse_args()

    # Load configuration
    config = load_config()

    filename, file_extension = os.path.splitext(os.path.basename(args.pdb_file_name))
    with_pdb_code = (file_extension == '')
    input_pdb_file = None

    if with_pdb_code:
        # PDB code has been specified instead of a PDB file name
        pdb_code = args.pdb_file_name
        input_pdb_file = '{0}.pdb'.format(pdb_code)
        if not os.path.exists(input_pdb_file):
            try:
                pdbutil.download_pdb_structure(pdb_code, input_pdb_file)
            except Exception as err:
                logger.error(str(err))
                raise SystemExit
        else:
            logger.warning("PDB structure already exists ({0}), no need to download it again".format(input_pdb_file))
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
        logger.error("Chain {0} provided not in available chains: {1}".format(chain_id, str(chain_ids)))
        raise SystemExit

    # Save only the given chain and discard residues with alternative positions
    io = PDBIO()
    current_pdb_file = "{0}_{1}.pdb".format(pdb_code, chain_id)
    for chain in structure.get_chains():
        if chain.id == chain_id:
            io.set_structure(chain)
            io.save(current_pdb_file, select=pdbutil.NotAlternative())
    logger.info("PDB structure with chain {0} saved to {1}".format(chain_id, current_pdb_file))

    # Calculate SASA:
    rsa_output_file = "{0}_{1}.rsa".format(pdb_code, chain_id)
    access.calculate_accessibility(current_pdb_file, rsa_output_file)
    logger.info("Atom accessibility calculated to {0}".format(rsa_output_file))

    # Calculate the different accessibility files according to the cutoffs:
    cutoffs = config['CUTOFF']
    access.create_cutoff_files(rsa_output_file, pdb_code, chain_id, cutoffs)
    logger.info("Surface and buried residues calculated")

    # Get structure sequence
    master_sequence = pdbutil.get_pdb_sequence(input_pdb_file, chain_id)
    write_to_fasta("{0}_{1}.fasta".format(filename, chain_id), master_sequence)

    hssp_file = "{0}.hssp".format(pdb_code)
    phylip_file = "{0}_{1}.phylseq".format(pdb_code, chain_id)

    if with_pdb_code:
        # Get the HSSP alignment from FTP if pdb code specified
        if not os.path.exists(hssp_file):
            logger.info("Downloading HSSP alignment...")
            try:
                compressed_hssp_file = hssp.get_from_ftp(pdb_code)
                hssp.decompress_bz2(compressed_hssp_file, hssp_file)
                logger.info("HSSP alignment stored to {0}".format(hssp_file))
            except Exception as err:
                logger.error("HSSP file could not be generated")
                raise SystemExit("Error is: {0}".format(err))
        try:
            hssp.hssp_file_to_phylip(hssp_file, phylip_file, chain_id, master_sequence)
        except Exception as err:
            logger.error(str(err))
            raise SystemExit
        
        logger.info("HSSP file converted to PHYLIP format")

    if not os.path.exists(hssp_file):
        # Run BLASTP if needed
        blast_output_file = "{0}_{1}_blast.xml".format(filename, chain_id)
        input_sequence_file = "{0}_{1}.fasta".format(filename, chain_id)
        if not os.path.exists(blast_output_file):
            logger.info("Please wait while running BLASTP against NCBI servers...")
            ncbi_blast(input_sequence_file, blast_output_file)
            logger.info("Result stored in {0}".format(blast_output_file))
        else:
            logger.info("BLAST file found ({0}), nothing to do".format(blast_output_file))

        # Convert file to FASTA format
        blast_qresult = SearchIO.read(blast_output_file, 'blast-xml')
        records = []
        for hit in blast_qresult:
            records.append(hit[0].hit)
        blast_fasta_file = "{0}_{1}_blast.fasta".format(filename, chain_id)
        SeqIO.write(records, blast_fasta_file, "fasta")

        # Multiple sequence alignment
        logger.info("MSA using MUSCLE...")
        output_alignment_file = "{0}_{1}_msa.fasta".format(filename, chain_id)
        msa = muscle_msa(config, blast_fasta_file, output_alignment_file)
        logger.info("Done.")

        # Convert MSA to Phylipseq
        logger.info("Converting MSA file to Phylseq format...")
        output_phylseq_file = "{0}_{1}.phylseq".format(filename, chain_id)
        msa_to_phylseq(msa, master_sequence, output_phylseq_file)
        logger.info("{} file written".format(output_phylseq_file))

    if not os.path.exists(phylip_file):
        logger.error("PHYLIP sequence file {0} not found".format(phylip_file))
        raise SystemExit

    # Calculate protdist
    protdist_output_file = "{0}_{1}.out".format(filename, chain_id)
    calculate_protdist(phylip_file, protdist_output_file)
    logger.info("Protdist calculated")

    # Generate conversion table file
    conv_output_file = "{0}_{1}.conv".format(filename, chain_id)
    pdbutil.map_protein_to_sequence_alignment(current_pdb_file, chain_id, master_sequence, conv_output_file)
    logger.info("Conversion table file generated")

    logger.info("Whiscy setup finished")
