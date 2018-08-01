#!/usr/bin/env python3

"""Whiscy predictor setup"""

import argparse
import os
import json
import Bio
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO, SeqIO
from libwhiscy import hssp


def load_config(config_file='etc/local.json'):
    """Load Whiscy configuration"""
    with open(config_file, 'r') as f:
        config = json.load(f)
        return config


def muscle_msa(config, input_sequence_file, output_alignment_file):
    """Calculates a MSA using MUSCLE's Biopython wrapper"""
    muscle_bin = config['ALIGN']['MUSCLE_BIN']
    muscle_cline = MuscleCommandline(muscle_bin, input=input_sequence_file, out=output_alignment_file)
    stdout, stderr = muscle_cline()
    MultipleSeqAlignment = AlignIO.read(output_alignment_file, "fasta") 
    return MultipleSeqAlignment


if __name__ == "__main__":

    # Parse command line
    parser = argparse.ArgumentParser(prog='whiscy_setup')
    parser.add_argument("pdb_file_name", help="PDB file name (.pdb extension) or PDB code", metavar="pdb_file_name")
    parser.add_argument("chain_id", help="Chain ID to be predicted", metavar="chain_id")
    args = parser.parse_args()

    # Load configuration
    config = load_config()

    filename, file_extension = os.path.splitext(os.path.basename(args.pdb_file_name))
    input_pdb_file = args.pdb_file_name

    # PDB code has been specified instead of a PDB file name
    if file_extension == '':
        pdb_code = args.pdb_file_name
        pdbl = PDBList()
        file_name = pdbl.retrieve_pdb_file(pdb_code, file_format='pdb', pdir='.', overwrite=True)
        if os.path.exists(file_name):
            input_pdb_file = '{0}.pdb'.format(pdb_code)
            os.rename(file_name, input_pdb_file)
        else:
            raise SystemExit("Error downloading structure: {0}".format(pdb_code))

    # Check if chain belongs to this PDB
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    structure = parser.get_structure(filename, input_pdb_file)
    chain_ids = [chain.id for chain in structure.get_chains()]
    chain_id = args.chain_id.upper()
    if len(chain_id) > 1:
        raise SystemExit("Wrong chain id {0}".format(chain_id))
    if chain_id not in chain_ids:
        raise SystemExit("Chain {0} provided not in available chains: {1}".format(chain_id, str(chain_ids)))

    # Get structure sequence
    with open(input_pdb_file, "rU") as handle:
        sequences = []
        for record in SeqIO.parse(handle, "pdb-seqres"):
            if record.id[-1] == chain_id:
                sequences.append(record)

        with open("{0}_{1}.fasta".format(filename, chain_id), "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")

    # Get the HSSP alignment from FTP
    print("Downloading HSSP alignment...")
    hssp_file = hssp.get_from_ftp(pdb_code)
    print("HSSP alignment stored to {}".format(hssp_file))
    
    input_sequence_file = "{0}_{1}.fasta".format(filename, chain_id)
    output_alignment_file = "{0}_{1}_msa.fasta".format(filename, chain_id)
    msa = muscle_msa(config, input_sequence_file, output_alignment_file)
