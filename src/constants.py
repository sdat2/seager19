"""This file is used to save all possible project wide constants.

Includes source folder, the project path, etc.

Example:
    Import statement at top of script::

        from src.constants import PROJECT_PATH, FIGURE_PATH, GWS_DIR

"""

# import os/pathlib to manipulate file names.
import os
import pathlib
from omegaconf import DictConfig


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
OCEAN_OUTPUT_PATH = OCEAN_PATH / "output"
ATMOS_PATH = PROJECT_PATH / "atmos"
ATMOS_DATA_PATH = ATMOS_PATH / "DATA"
ATMOS_TMP_PATH = ATMOS_PATH / "tmp"
GIF_PATH = PROJECT_PATH / "gifs"

# General data from e.g. paper or cmip etc.
DATA_PATH = SRC_PATH / "data"
CMIP_TS_PATH = DATA_PATH / "ts_nc"

# Wandb-summary file download:
ORIG_WANDB_DATA = DATA_PATH / "results.csv"
NEW_WANDB_DATA = DATA_PATH / "results22Jun.csv"

# Test DIREC
TEST_DIREC = SRC_PATH / "test" / "test_direc"

# Original paper figure data.
FIGURE_DATA_NAME: str = "Seager_etal_NCC-2019_datasetdatafiles.nc"
FIGURE_DATA_PATH: pathlib.Path = DATA_PATH / FIGURE_DATA_NAME

# region selection dictionary
SEL_DICT = {
    # [west, east] [south, north] boundaries
    # X in degrees east, Y in degrees north.
    # colors compatible with matplotlib.
    "pac": {"X": (100, 290), "Y": (-30, 30), "color": "#411900"},
    "nino1+2": {"X": (270, 280), "Y": (-10, 0), "color": "#0652ff"},
    "nino3": {"X": (210, 270), "Y": (-5, 5), "color": "#01a049"},
    "nino3.4": {"X": (190, 240), "Y": (-5, 5), "color": "#380282"},
    "nino4": {"X": (160, 210), "Y": (-5, 5), "color": "#fe01b1"},
}

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

OLD_LOGS = GWS_DIR / "logs"
CD_LOGS = GWS_DIR / "sensitivity" / "cd_logs"
NEW_LOGS = PROJECT_PATH / "logs"
K_LOGS = GWS_DIR / "sensitivity" / "k_days_logs"
EPS_LOGS = GWS_DIR / "sensitivity"/ "eps_days_logs"


# pylint: disable=using-constant-test
# if False:  # os.path.exists(GWS_DIR):
#    LOG_PATH = GWS_DIR / "logs"
# else:
LOG_PATH = PROJECT_PATH / "logs"
FIN_LOG_PATH = EPS_LOGS
# LOG_PATH = PROJECT_PATH / "k_days_logs"

if not os.path.exists(str(LOG_PATH)):
    os.mkdir(str(LOG_PATH))

# Report WIDTH
REPORT_WIDTH: float = 398.3386  # in pixels

# DATE FORMAT for plotting titles
DATE_TITLE_FORMAT: str = "%Y.%m.%d"


def run_path(cfg: DictConfig, unit_test: bool = False) -> str:
    """
    Returns run path to store data in.

    Args:
        cfg (DictConfig): The config struct.
        unit_test (bool, optional): Whether this is a unit test. Defaults to False.

    Returns:
        str: The path to the relevant directory that exists.
    """
    if not unit_test:
        run_dir = os.path.join(LOG_PATH, cfg.name)
    else:
        run_dir = str(TEST_DIREC)
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)
    return run_dir
