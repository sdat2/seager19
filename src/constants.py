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

# PATHS to the models
OCEAN_PATH = PROJECT_PATH / "ocean"
OCEAN_DATA_PATH = OCEAN_PATH / "DATA"
OCEAN_RUN_PATH = OCEAN_PATH / "RUN"
OCEAN_SRC_PATH = OCEAN_PATH / "SRC"
OCEAN_OUTPUT_PATH = OCEAN_PATH / "output"
ATMOS_PATH = PROJECT_PATH / "atmos"
ATMOS_DATA_PATH = ATMOS_PATH / "DATA"
ATMOS_TMP_PATH = ATMOS_PATH / "tmp"

# Data directory on GWS
GWS_DIR = pathlib.Path("/gws/nopw/j04/ai4er/users/sdat2")

# mv wandb /gws/nopw/j04/ai4er/users/sdat2/
# mv OLD /gws/nopw/j04/ai4er/users/sdat2/
# mv tmp /gws/nopw/j04/ai4er/users/sdat2/
# mv nemo_ba130o_1h_20040101-20040131_grid-combined.nc /gws/nopw/j04/ai4er/users/sdat2/
# mv netcdf /gws/nopw/j04/ai4er/users/sdat2/
# mv gtc-biodiversity /gws/nopw/j04/ai4er/users/sdat2/
# mv surge /gws/nopw/j04/ai4er/users/sdat2/
# mv etc /gws/nopw/j04/ai4er/users/sdat2/
# mv Pytorch-UNet /gws/nopw/j04/ai4er/users/sdat2/
# mv sosie /gws/nopw/j04/ai4er/users/sdat2/
# mv milesial_Pytorch-UNet_master /gws/nopw/j04/ai4er/users/sdat2/
# mv feb.sh /gws/nopw/j04/ai4er/users/sdat2/
# mv import_katrina.sh /gws/nopw/j04/ai4er/users/sdat2/
# mv read_xarray.py /gws/nopw/j04/ai4er/users/sdat2/
# mv retry_download.sh /gws/nopw/j04/ai4er/users/sdat2/

# Report WIDTH
REPORT_WIDTH = 398.3386  # in pixels

# DATE FORMAT for plotting titles
DATE_TITLE_FORMAT = "%Y.%m.%d      %H:%M"
