"""WHISCY: protein-protein interface prediction from conservation."""

import logging
import shlex
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

try:
    __version__ = version("whiscy")
except PackageNotFoundError:
    __version__ = "0+unknown"

PARAM_PATH = Path(Path(__file__).parent, "param")
LOG_FILE = "whiscy.log"


def log_to_file(logger):
    """Append the messages of logger to whiscy.log, after the version and command."""
    handler = logging.FileHandler(LOG_FILE)
    handler.setLevel(logging.INFO)
    handler.setFormatter(logger.handlers[0].formatter)
    logger.addHandler(handler)
    command = shlex.join([Path(sys.argv[0]).name] + sys.argv[1:])
    logger.info("WHISCY {}: {}".format(__version__, command))
