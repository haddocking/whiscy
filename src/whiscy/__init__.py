"""WHISCY: protein-protein interface prediction from conservation."""

from importlib.metadata import version
from pathlib import Path

__version__ = version("whiscy")

PARAM_PATH = Path(Path(__file__).parent, "param")
