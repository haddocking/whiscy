from ftplib import FTP
import os
import bz2


def get_from_ftp(pdb_code, path_to_store='.',
                 ftp_server='ftp.cmbi.ru.nl', ftp_path='/pub/molbio/data/hssp/'):
    """Downloads using FTP protocol an HSSP alignment for the given pdb_code"""
    # Start connection
    ftp = FTP(ftp_server)
    # Anonymous login
    ftp.login()
    # Move to path where HSSP alignments are stored
    ftp.cwd(ftp_path)
    # File name format
    file_name = '{}.hssp.bz2'.format(pdb_code.lower())
    # Retrieve file
    path_to_file = os.path.join(path_to_store, file_name)
    ftp.retrbinary("RETR " + file_name, open(path_to_file, 'wb').write)
    # Close connection
    ftp.close()

    return path_to_file


def decompress_bz2(file_name_input, file_name_output):
    """Decompresses file_name_input in BZ2 format into file_name_output""" 
    with open(file_name_output, 'wb') as new_file, bz2.BZ2File(file_name_input, 'rb') as file:
        for data in iter(lambda : file.read(100 * 1024), b''):
            new_file.write(data)


def _parse_hssp_proteins(line_buffer):
    """Parses the '## PROTEINS section' from HSSP alignment file"""
    proteins = {}
    for line in line_buffer:
        if line.startswith("  NR.") or line.startswith("##"):
            continue
        # Only get the id and name of the protein in the alignment
        fields = (line[:20]).split(':')
        seq_id = int(fields[0]) - 1
        name = fields[1].strip()[:10]
        proteins[seq_id] = name
    return proteins


def _parse_hssp_alignments(line_buffer, chain_id, num_alignments):
    """Parses the '## ALIGNMENTS section' from HSSP alignment file"""
    alignments = [[] for i in range(num_alignments)]
    first_alignment = 0
    last_alignment = 0
    for line in line_buffer:
        if line.startswith(" SeqNo") or line[12] == '!':
            continue
        if line.startswith("## ALIGNMENTS"):
            fields = (line[13:]).split('-')
            # We are now parsing alignments from first to last specified
            # in the ALINGMENTS header
            first_alignment = int(fields[0]) - 1
            last_alignment = int(fields[1]) - 1
        else:
            if line[12] == chain_id:
                for i, s in enumerate(line[51:]):
                    # We will convert spaces or dots to -
                    if s == '.' or s == ' ':
                        s = '-'
                    else:
                        # Just in case to avoid insertion characters in minor case
                        s = s.upper()
                    alignments[first_alignment + i].append(s)
    alignments = [''.join(s) for s in alignments]
    return alignments


def hssp_file_to_phylip(hssp_file_name, phylip_file_name, chain_id, master_sequence):
    """Parses an HSSP file and returns a list of the sequences"""
    # We're only instered in the lenght of the sequence of our given chain_id, 
    # SEQLENGHT header gives us the sum of all.
    seqlength = len(master_sequence)
    num_alignments = 0
    num_chains = 0
    parsing = False
    parsing_alignment = False
    line_buffer = []
    parsing_proteins = False
    prot_line_buffer = []
    with open(hssp_file_name, "rU") as handle:
        for line in handle:
            line = line.rstrip(os.linesep)
            if line.startswith('NCHAIN'):
                num_chains = int(line.split()[1])
            if line.startswith('NALIGN'):
                num_alignments = int(line.split()[1])
            
            parsing = (seqlength != 0 and num_chains != 0 and num_alignments != 0)
            
            if parsing:
                if line.startswith("## ALIGNMENTS"):
                    parsing_alignment = True

                if line.startswith("## PROTEINS"):
                    parsing_proteins = True

                if line.startswith("##") and "ALIGNMENTS" not in line:
                    parsing_alignment = False

                if line.startswith("##") and "PROTEINS" not in line:
                    parsing_proteins = False

                if parsing_alignment:
                    line_buffer.append(line)

                if parsing_proteins:
                    prot_line_buffer.append(line)

        proteins = _parse_hssp_proteins(prot_line_buffer)
        alignments = _parse_hssp_alignments(line_buffer, chain_id.upper(), num_alignments)

        non_valid = [k for k in proteins.keys() if alignments[k].count('-') == seqlength]
        with open(phylip_file_name, 'w') as output_handle:
            # Write header
            output_handle.write("{}  {}{}".format(len(proteins) - len(non_valid), seqlength, os.linesep))
            # Write master sequence
            output_handle.write("MASTER    {}{}".format(master_sequence, os.linesep))
            # Write the rest of non null alignments
            for k in sorted(proteins.keys()):
                if k not in non_valid:
                    output_handle.write("{:10s}{}{}".format(proteins[k], alignments[k], os.linesep))
