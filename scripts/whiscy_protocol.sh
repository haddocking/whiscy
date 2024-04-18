#!/bin/bash

protein=$1

# Default protocol
whiscy.py $protein.sur $protein.conv $protein.phylseq $protein.out -o $protein.cons

# Interface propensities
whiscy.py $protein.lac $protein.conv $protein.phylseq $protein.out -o $protein.lcons
${WHISCY_PATH}/bin/consadjust.py $protein.cons ${WHISCY_PATH}/param/weight_ma.txt ${WHISCY_PATH}/param/ztable.txt -o $protein.acons

# Surface smoothing
${WHISCY_PATH}/bin/residue_distance.py $protein.pdb $protein.conv $protein.rd
${WHISCY_PATH}/bin/parasmooth.py $protein.acons $protein.lcons $protein.rd ${WHISCY_PATH}/param/parasmooth.par -o $protein.pscons
