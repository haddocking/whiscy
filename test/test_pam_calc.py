import os
from libwhiscy.pam_calc import pam_load_sequences, pam_calc_similarity


golden_data_path = os.path.normpath(os.path.dirname(os.path.realpath(__file__))) + '/golden_data/'


def test_pam_load_sequences():
    alignment_file = os.path.join(golden_data_path, 'pam', '2SNIE.phylseq')
    distance_file = os.path.join(golden_data_path, 'pam', '2SNIE.out')

    reference_sequence = ("AQSVPYGVSQIKAPALHSQGYTGSNVKVAVIDSGIDSSHPDLKVAGGASMVPSETNPFQDNNSHGTHVAGTVAAL"
                          "NNSIGVLGVAPSASLYAVKVLGADGSGQYSWIINGIEWAIANNMDVINMSLGGPSGSAALKAAVDKAVASGVVVVAAAGNEGTSGSSSTVGYPGKY"
                          "PSVIAVGAVDSSNQRASFSSVGPELDVMAPGVSIQSTLPGNKYGAYNGTSMASPHVAGAAALILSKHPNWTNTQVRSSLENTTTKLGDSFYYGKGLINVQAAAQ")

    seqnr, seqlen, refseq, distances, sequences, seqtodis = pam_load_sequences(alignment_file, distance_file)

    assert reference_sequence == refseq


def test_pam_calc_similarity():
    alignment_file = os.path.join(golden_data_path, 'pam', '2SNIE.phylseq')
    distance_file = os.path.join(golden_data_path, 'pam', '2SNIE.out')

    seqnr, seqlen, refseq, seq_distances, sequences, seqtodis = pam_load_sequences(alignment_file, distance_file)

    posnr, distances, scores = pam_calc_similarity(2, seqnr, sequences, seq_distances)

    assert posnr == 267
