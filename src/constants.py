"""This file is used to save all possible project wide constants.

Includes source folder, the project path, etc.

Example:
    Import statement at top of script::

        from src.constants import PROJECT_PATH, FIGURE_PATH, GWS_DIR
"""

# Place all your constants here
import os
import pathlib

# Note: constants should be UPPER_CASE
constants_path = pathlib.Path(os.path.realpath(__file__))
SRC_PATH = pathlib.Path(os.path.dirname(constants_path))
PROJECT_PATH = pathlib.Path(os.path.dirname(SRC_PATH))
REPORT_PATH = pathlib.Path(os.path.join(PROJECT_PATH, "report"))
FIGURE_PATH = pathlib.Path(os.path.join(REPORT_PATH, "figures"))


OCEAN_PATH = PROJECT_PATH / "ocean"
OCEAN_DATA_PATH = OCEAN_PATH / "DATA"
OCEAN_OUTPUT_PATH = OCEAN_PATH / "RUN" / "output"
ATMOS_PATH = PROJECT_PATH / "atmos"
ATMOS_DATA_PATH = ATMOS_PATH / "DATA"
ATMOS_TMP_PATH = ATMOS_PATH / "tmp"

# Data directory on GWS
GWS_DIR = pathlib.Path("/gws/nopw/j04/ai4er")

# Report WIDTH
REPORT_WIDTH = 398.3386  # in pixels
