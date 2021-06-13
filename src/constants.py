"""This file is used to save all possible project wide constants.

Includes source folder, the project path, etc.

Example:
    Import statement at top of script::

        from src.constants import PROJECT_PATH, FIGURE_PATH, GWS_DIR

"""

# import os/pathlib to manipulate file names.
import os
import pathlib

# Note: constants should be UPPER_CASE
constants_path = pathlib.Path(os.path.realpath(__file__))
SRC_PATH = pathlib.Path(os.path.dirname(constants_path))
PROJECT_PATH = pathlib.Path(os.path.dirname(SRC_PATH))
REPORT_PATH = pathlib.Path(os.path.join(PROJECT_PATH, "report"))
FIGURE_PATH = pathlib.Path(os.path.join(REPORT_PATH, "figures"))
CONFIG_PATH = os.path.join(SRC_PATH, "configs")
CONFIG_NAME: str = "config"

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

# General data from e.g. paper or cmip etc.
DATA_PATH = SRC_PATH / "data"

# Test DIREC
TEST_DIREC = SRC_PATH / "test" / "test_direc"

# Figure data.
FIGURE_DATA_NAME: str = "Seager_etal_NCC-2019_datasetdatafiles.nc"
FIGURE_DATA_PATH: pathlib.Path = DATA_PATH / FIGURE_DATA_NAME

# NOAA DATA
NOAA_DATA_NAME: str = "NOAA_NCDC_ERSST_v3b_SST.nc"
NOAA_DATA_PATH: pathlib.Path = DATA_PATH / NOAA_DATA_NAME

# Nino3.4 test results from the noaa data.
# get_test_nino_data
NINO3_4_TEST_CODE: str = "8j698fap5iq2v9y/"
NINO3_4_TEST_NAME: str = "noaa_nino3_4.nc"
NINO3_4_TEST_PATH: pathlib.Path = DATA_PATH / NINO3_4_TEST_NAME

# Data directory on GWS
GWS_DIR = pathlib.Path("/gws/nopw/j04/ai4er/users/sdat2")

# pylint: disable=using-constant-test
if False:  # os.path.exists(GWS_DIR):
    LOG_PATH = GWS_DIR / "logs"
else:
    LOG_PATH = PROJECT_PATH / "logs"

if not os.path.exists(str(LOG_PATH)):
    os.mkdir(str(LOG_PATH))

# Report WIDTH
REPORT_WIDTH: float = 398.3386  # in pixels

# DATE FORMAT for plotting titles
DATE_TITLE_FORMAT: str = "%Y.%m.%d"
