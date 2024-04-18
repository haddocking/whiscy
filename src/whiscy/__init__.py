# Define the third-party executables
import os

# Define the path to the third-party executables
FREESASA_BIN = os.environ.get("FREESASA_BIN")
MUSCLE_BIN = os.environ.get("MUSCLE_BIN")
HSSPCONV_BIN = os.environ.get("HSSPCONV_BIN")
PROTDIST_BIN = os.environ.get("PROTDIST_BIN")

# Make sure the path to the third-party executables is set
if not FREESASA_BIN:
    raise Exception("FREESASA_BIN environment variable is not set")

if not MUSCLE_BIN:
    raise Exception("MUSCLE_BIN environment variable is not set")

if not HSSPCONV_BIN:
    raise Exception("HSSPCONV_BIN environment variable is not set")

if not PROTDIST_BIN:
    raise Exception("PROTDIST_BIN environment variable is not set")

CUTOFF = {
    "sa_pred_cutoff": 15.0,
    "sa_act_cutoff": 40.0,
    "air_cutoff": 0.18,
    "air_dist_cutoff": 6.5,
}

AIR = {
    "air_pro_percentage": 10.0,
    "air_wm_pro_or": 98.52,
    "air_wm_whis_or": 0.370515,
    "air_wm_pro_and": 55.42,
    "air_wm_whis_and": 0.106667,
}
