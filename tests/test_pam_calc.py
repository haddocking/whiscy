import os

import pytest

from libwhiscy.pam_calc import pam_calc_similarity, pam_load_sequences

from . import (GOLDEN_DATA_PATH, PAM_EXPECTED_DISTANCES, PAM_EXPECTED_SCORES,
               PAM_EXPECTED_SEQTODIS)


@pytest.fixture
def alignment_file():
    yield os.path.join(GOLDEN_DATA_PATH, "pam", "2SNIE.phylseq")


@pytest.fixture
def distance_file():
    yield os.path.join(GOLDEN_DATA_PATH, "pam", "2SNIE.out")


@pytest.fixture
def reference_sequence():
    yield (
        "AQSVPYGVSQIKAPALHSQGYTGSNVKVAVIDSGIDSSHPDLKVAGGASMVPSETNPFQDNNSHGTHVAGTVAAL"
        "NNSIGVLGVAPSASLYAVKVLGADGSGQYSWIINGIEWAIANNMDVINMSLGGPSGSAALKAAVDKAVASGVVVVAAAGNEGTSGSSSTVGYPGKY"
        "PSVIAVGAVDSSNQRASFSSVGPELDVMAPGVSIQSTLPGNKYGAYNGTSMASPHVAGAAALILSKHPNWTNTQVRSSLENTTTKLGDSFYYGKGLINVQAAAQ"
    )


def test_pam_load_sequences(reference_sequence, alignment_file, distance_file):

    seqnr, seqlen, refseq, _, _, seqtodis = pam_load_sequences(
        alignment_file, distance_file
    )

    assert seqnr == 629
    assert seqlen == 275
    assert refseq == reference_sequence

    assert len(seqtodis) == len(PAM_EXPECTED_SEQTODIS)

    for i, j in zip(seqtodis, PAM_EXPECTED_SEQTODIS):
        assert i == j


def test_pam_calc_similarity(reference_sequence, alignment_file, distance_file):

    seqnr, seqlen, refseq, seq_distances, sequences, seqtodis = pam_load_sequences(
        alignment_file, distance_file
    )

    assert seqnr == 629
    assert seqlen == 275
    assert refseq == reference_sequence

    assert len(seqtodis) == len(PAM_EXPECTED_SEQTODIS)

    for i, j in zip(seqtodis, PAM_EXPECTED_SEQTODIS):
        assert i == j

    posnr, distances, scores = pam_calc_similarity(2, seqnr, sequences, seq_distances)

    assert posnr == 267

    for d in distances:
        assert d == pytest.approx(PAM_EXPECTED_DISTANCES.pop(0), abs=0.01)

    for s in scores:
        assert s == pytest.approx(PAM_EXPECTED_SCORES.pop(0), abs=0.01)
