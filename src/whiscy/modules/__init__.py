import os
from pathlib import Path

MUSCLE_BIN = os.environ.get("MUSCLE_BIN")
FREESASA_BIN = os.environ.get("FREESASA_BIN")
HSSPCONV_BIN = os.environ.get("HSSPCONV_BIN")
PROTDIST_BIN = os.environ.get("PROTDIST_BIN")

# Make sure none of them are None
if MUSCLE_BIN is None:
    raise ValueError("MUSCLE_BIN not found in system variables")

if FREESASA_BIN is None:
    raise ValueError("FREESASA_BIN not found in system variables")

if HSSPCONV_BIN is None:
    raise ValueError("HSSPCONV_BIN not found in system variables")

if PROTDIST_BIN is None:
    raise ValueError("PROTDIST_BIN not found in system variables")

PARAM_PATH = Path(Path(__file__).parent.parent, "param")


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
