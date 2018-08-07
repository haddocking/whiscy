import subprocess


def calculate_accessibility(pdb_file_name, output_file_name):
    """Calculates the SASA using freesasa.

    Uses the command line interface and not the Python bindings to be able to get 
    a RSA NACCESS-format like file.
    """
    cmd = "freesasa {} -n 20 --format=rsa --radii=naccess -o {}".format(pdb_file_name, output_file_name)
    subprocess.run(cmd, shell=True)
