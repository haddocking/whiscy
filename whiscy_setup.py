#!/usr/bin/env python3

"""Whiscy predictor setup"""

import argparse
import os
import json
import Bio
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.Blast import NCBIWWW
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO, SeqIO
from libwhiscy import hssp
import warnings
# Import SearchIO and suppress experimental warning
from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)
    from Bio import SearchIO


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
        output_handle.write("{}  {}{}".format(len(msa), 
                                              len(master_sequence),
                                              os.linesep))

        # Write master sequence
        output_handle.write("MASTER    {}{}".format(master_sequence,
                                                    os.linesep))
        # Write the rest of alignments
        for alignment in msa:
            output_handle.write("{:10s}{}{}".format(alignment.id[:10],
                                                    alignment.seq,
                                                    os.linesep))


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
        input_pdb_file = '{0}.pdb'.format(pdb_code)
        if not os.path.exists(input_pdb_file):
            pdbl = PDBList()
            file_name = pdbl.retrieve_pdb_file(pdb_code, file_format='pdb', pdir='.', overwrite=True)
            if os.path.exists(file_name):
                os.rename(file_name, input_pdb_file)
            else:
                raise SystemExit("Error downloading structure: {0}".format(pdb_code))
        else:
            print("PDB structure already exists ({}), no need to download it again".format(input_pdb_file))

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

    # Get the HSSP alignment from FTP if pdb code specified
    hssp_file = "{}.hssp".format(pdb_code)
    if pdb_code and not os.path.exists(hssp_file):
        print("Downloading HSSP alignment...")
        compressed_hssp_file = hssp.get_from_ftp(pdb_code)
        hssp.decompress_bz2(compressed_hssp_file, hssp_file)
        print("HSSP alignment stored to {}".format(hssp_file))
    
    # Run BLASTP if needed
    blast_output_file = "{0}_{1}_blast.xml".format(filename, chain_id)
    input_sequence_file = "{0}_{1}.fasta".format(filename, chain_id)
    if not os.path.exists(blast_output_file):
        print("Please wait while running BLASTP against NCBI servers...")
        ncbi_blast(input_sequence_file, blast_output_file)
        print("Result stored in {}".format(blast_output_file))
    else:
        print("BLAST file found ({}), nothing to do".format(blast_output_file))

    # Convert file to FASTA format
    blast_qresult = SearchIO.read(blast_output_file, 'blast-xml')
    records = []
    for hit in blast_qresult:
        records.append(hit[0].hit)
    blast_fasta_file = "{0}_{1}_blast.fasta".format(filename, chain_id)
    SeqIO.write(records, blast_fasta_file, "fasta")

    # Multiple sequence alignment
    print("MSA using MUSCLE...")
    output_alignment_file = "{0}_{1}_msa.fasta".format(filename, chain_id)
    msa = muscle_msa(config, blast_fasta_file, output_alignment_file)
    print("Done.")

    # Convert MSA to Phylipseq
    print("Converting MSA file to Phylseq format...")
    output_phylseq_file = "{0}_{1}.phylseq".format(filename, chain_id)
    master_sequence = SeqIO.read(input_sequence_file, format="fasta")
    msa_to_phylseq(msa, master_sequence.seq, output_phylseq_file)
    print("{} file written".format(output_phylseq_file))
