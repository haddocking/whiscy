#!/usr/bin/env python3

import math
import os
import argparse
from libwhiscy.whiscy_data import load_residue_weights
# Logging
import logging
logging.basicConfig(format='%(name)s [%(levelname)s] %(message)s', level=logging.INFO)
logger = logging.getLogger("consadjust")


if __name__ == "__main__":

    # Parse command line
    parser = argparse.ArgumentParser(prog='parasmooth')
    parser.add_argument("cons_file", help="Conservation file", metavar="cons_file")
    parser.add_argument("residue_weight_file", help="Residue weight file", metavar="residue_weight_file")
    parser.add_argument("z_table_file", help="Z-table file", metavar="z_table_file")
    parser.add_argument("-o", "--output", help="If set, output prediction to this file", 
                        dest="output_file", metavar="output_file")
    args = parser.parse_args()

    if not os.path.exists(args.cons_file):
        logger.error("Conservation file {} does not exist".format(args.cons_file))
        raise SystemExit

    if not os.path.exists(args.residue_weight_file):
        logger.error("Weight file {} does not exist".format(args.residue_weight_file))
        raise SystemExit    

    logger.info("Reading input files")

    resweight = load_residue_weights(args.residue_weight_file)

    # If output to file
    if args.output_file:
        with open(args.output_file, 'w') as output_handle:
            output_handle.write()

    else:
        pass