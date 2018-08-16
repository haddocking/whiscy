import subprocess
import os


def calculate_accessibility(pdb_file_name, output_file_name):
    """Calculates the SASA using freesasa.

    Uses the command line interface and not the Python bindings to be able to get 
    a RSA NACCESS-format like file.
    """
    cmd = "freesasa {} -n 20 --format=rsa --radii=naccess -o {}".format(pdb_file_name, output_file_name)
    try:
        subprocess.run(cmd, shell=True)
    except:
        subprocess.check_call(cmd, shell=True)


class ResidueSASA():
    def __init__(self, chain_id, name, number, tot_rel, sd_rel, bk_rel):
        self.chain_id = chain_id
        self.name = name
        self.number = number
        self.tot_rel = tot_rel
        self.sd_rel = sd_rel
        self.bk_rel = bk_rel


def parse_rsa_file(rsa_file_name):
    """Parses a .rsa NACCESS (or freesasa) file and gets the relative SASAs"""
    residue_sasas = []
    with open(rsa_file_name, "rU") as input_handle:
        for line in input_handle:
            if line.startswith("RES"):
                if line[13] == ' ':
                    # Avoid alternative positions
                    name = line[4:7]
                    chain_id = line[8]
                    number = int(line[9:13])
                    # Difference between NACCESS and freesasa output
                    try:
                        tot_rel = float(line[23:29])
                    except ValueError:
                        tot_rel = -99.9
                    try:
                        sd_rel = float(line[36:42])
                    except ValueError:
                        sd_rel = -99.9
                    try:
                        bk_rel = float(line[49:55])
                    except ValueError:
                        bk_rel = -99.9

                    residue_sasas.append(ResidueSASA(chain_id, name, number,
                                                     tot_rel, sd_rel, bk_rel))

    return residue_sasas


def create_cutoff_files(rsa_file_name, pdb_code, chain_id, cutoffs, path='.'):
    """Creates three output files depending on cutoffs of accessibility.

    - pdb_code.sur: rsa_file_name filtered by surface sa_pred_cutoff
    - pdb_code.suract: rsa_file_name filtered by surface sa_act_cutoff
    - pdb_code.lac: rsa_file_name filtered by buried sa_pred_cutoff
    """
    residue_sasas = parse_rsa_file(rsa_file_name)

    # Create .sur file
    output_file_name = os.path.join(path, "{0}_{1}.sur".format(pdb_code, chain_id))
    cutoff = cutoffs['sa_pred_cutoff']
    with open(output_file_name, 'w') as output_handle:
        for res_sasa in residue_sasas:
            if res_sasa.chain_id == chain_id:
                if res_sasa.tot_rel >= cutoff or res_sasa.sd_rel >= cutoff or res_sasa.bk_rel >= cutoff:
                    output_handle.write("{}{}".format(res_sasa.number, os.linesep))

    # Create .suract file
    output_file_name = os.path.join(path, "{0}_{1}.suract".format(pdb_code, chain_id))
    cutoff = cutoffs['sa_act_cutoff']
    with open(output_file_name, 'w') as output_handle:
        for res_sasa in residue_sasas:
            if res_sasa.chain_id == chain_id:
                if res_sasa.tot_rel >= cutoff or res_sasa.sd_rel >= cutoff or res_sasa.bk_rel >= cutoff:
                    output_handle.write("{}{}".format(res_sasa.number, os.linesep))

    # Create .lac file
    output_file_name = os.path.join(path, "{0}_{1}.lac".format(pdb_code, chain_id))
    cutoff = cutoffs['sa_pred_cutoff']
    with open(output_file_name, 'w') as output_handle:
        for res_sasa in residue_sasas:
            if res_sasa.chain_id == chain_id:
                if (res_sasa.tot_rel < cutoff and res_sasa.tot_rel > 0) and \
                    res_sasa.sd_rel < cutoff and res_sasa.bk_rel < cutoff:
                    output_handle.write("{}{}".format(res_sasa.number, os.linesep))
