"""WHISCY: protein-protein interface prediction from conservation."""

from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

try:
    __version__ = version("whiscy")
except PackageNotFoundError:
    __version__ = "0+unknown"

PARAM_PATH = Path(Path(__file__).parent, "param")
