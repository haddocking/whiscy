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
