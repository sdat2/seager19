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
CONFIG_PATH = os.path.join(SRC_PATH, "configs")
CONFIG_NAME = "config"

# PATHS to the models
OCEAN_PATH = PROJECT_PATH / "ocean"
OCEAN_DATA_PATH = OCEAN_PATH / "DATA"
OCEAN_RUN_PATH = OCEAN_PATH / "RUN"
OCEAN_SRC_PATH = OCEAN_PATH / "SRC"
OCEAN_TS_PATH = OCEAN_PATH / "ts_nc"
OCEAN_OUTPUT_PATH = OCEAN_PATH / "output"
ATMOS_PATH = PROJECT_PATH / "atmos"
ATMOS_DATA_PATH = ATMOS_PATH / "DATA"
ATMOS_TMP_PATH = ATMOS_PATH / "tmp"
GIF_PATH = PROJECT_PATH / "gifs"
DATA_PATH = PROJECT_PATH / "data"

# Figure data.
FIGURE_DATA_NAME: str = "Seager_etal_NCC-2019_datasetdatafiles.nc"
FIGURE_DATA_PATH = DATA_PATH / "Seager_etal_NCC-2019_datasetdatafiles.nc"

# Data directory on GWS
GWS_DIR = pathlib.Path("/gws/nopw/j04/ai4er/users/sdat2")

# Report WIDTH
REPORT_WIDTH: float = 398.3386  # in pixels

# DATE FORMAT for plotting titles
DATE_TITLE_FORMAT: str = "%Y.%m.%d"
