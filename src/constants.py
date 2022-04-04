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
SENS_RANGES: str = "sens_ranges"
SENS_SETTINGS: str = "sens_settings"

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
MASK = OCEAN_DATA_PATH / "om_mask.nc"

# General data from e.g. paper or cmip etc.
DATA_PATH = SRC_PATH / "data"
NC_PATH = DATA_PATH / "nc"
CMIP_TS_PATH = DATA_PATH / "ts_nc"
CMIP6_TS_PATH = DATA_PATH / "nc80"
CMIP6_CLIM60_PATH = DATA_PATH / "nc_mean"

# NINO34 trend from fig5e
NINO34_TRENDS = DATA_PATH / "nino34-trends.csv"
NINO34_TRENDS_CMIP5 = DATA_PATH / "nino34-trends-cmip5.csv"

# Wandb-summary file download:
ORIG_WANDB_DATA = DATA_PATH / "results.csv"
NEW_WANDB_DATA = DATA_PATH / "results22Jun.csv"

# Multi model mean u and v fields
UV_PATH = DATA_PATH / "hist-winds-cmip5"
U_HIST = UV_PATH / "ua.nc"
V_HIST = UV_PATH / "va.nc"

# Multi model mean surface fields
MMM_V23_PATH = DATA_PATH / "mmm-v2.3-full-rep"
MMM_V23_HIST = MMM_V23_PATH / "cmip5-mmm-v2.3-historical.nc"
MMM_V23_RCP85 = MMM_V23_PATH / "cmip5-mmm-v2.3-rcp85.nc"

# Test DIREC
TEST_DIREC = SRC_PATH / "test" / "test_direc"

# Original paper figure data.
FIGURE_DATA_NAME: str = "Seager_etal_NCC-2019_datasetdatafiles.nc"
FIGURE_DATA_PATH: pathlib.Path = DATA_PATH / FIGURE_DATA_NAME

# NOAA DATA
NOAA_DATA_NAME: str = "NOAA_NCDC_ERSST_v3b_SST.nc"
NOAA_DATA_PATH: pathlib.Path = DATA_PATH / NOAA_DATA_NAME

# PSL DATA:
PSL_INDICES_NAME: str = "PSL_INDICES.nc"
PSL_INDICES_PATH: pathlib.Path = DATA_PATH / PSL_INDICES_NAME

# ERASSTv5
ERSSTV5_NAME: str = "ERSSTv5.sst.nc"
ERSSTV5_PATH: pathlib.Path = DATA_PATH / ERSSTV5_NAME

# Nino3.4 test results from the noaa data.
# get_test_nino_data
NINO3_4_TEST_CODE: str = "8j698fap5iq2v9y/"
NINO3_4_TEST_NAME: str = "noaa_nino3_4.nc"
NINO3_4_TEST_PATH: pathlib.Path = DATA_PATH / NINO3_4_TEST_NAME

# Model names:
MODEL_NAMES = {
    "E": "ECMWF",
    "F": "ECMWF-orig",
    "B": "CMIP5-39m",
    "C": "CMIP5",
    "6": "CMIP6",
    "S": "CMIP6",
    "D": "CMIP5-orig",
    "H": "HadGEM2",
    "f": "fixed",
    "e": "fixed78",
    "g": "fixed82",
    "W": "WHOI",
    "M": "MERRA",
    "I": "ISCCP",
}
VAR_DICT = {0: "ts", 1: "clt", 2: "sfcWind", 3: "rh", 4: "pr", 5: "ps", 6: "tau"}
# backwards compatibility: we want the new data to be stored without a atm.mem,
# but being able to process the old data where atm.mem was used.
#

ENSEMBLE_CSV = DATA_PATH / "ensemble_variable_members.csv"
# MINIMAL_ENSEMBLE_CSV = 

def atmos_input_file_path(
    var: str = "ts", model: str = "E", ending: str = "clim60"
) -> str:
    return str(
        ATMOS_DATA_PATH / str(var + "-" + MODEL_NAMES[model] + "-" + ending + ".nc")
    )


drop_var_d: dict = {"nc_clt": [], "nc_hur": [], "nc_pr": [], "nc_ts": []}
# https://www.dropbox.com/sh/pzp2s534m1i3081/AABsVz0HvpQTtXxlOXUS4eIla?dl=1
# https://www.dropbox.com/s/o82yp69pkpz50ze/nc_clt.zip?dl=0
# names of folders to download.
# Data directory on GWS
GWS_DIR = pathlib.Path("/gws/nopw/j04/ai4er/users/sdat2")
ARCHIVE_DIR = GWS_DIR / "rep"

# Subdirectories in GWS
OLD_LOGS = GWS_DIR / "logs"
CD_LOGS = GWS_DIR / "sensitivity" / "cd_logs"
NEW_LOGS = PROJECT_PATH / "logs"
K_LOGS = GWS_DIR / "sensitivity" / "k_days_logs"
EPS_LOGS = GWS_DIR / "sensitivity" / "eps_days_logs"
EPS_FRAC_LOGS = GWS_DIR / "sensitivity" / "eps_frac"
UC_LOGS = GWS_DIR / "uc_logs"

# pylint: disable=using-constant-test
# if False:  # os.path.exists(GWS_DIR):
#    LOG_PATH = GWS_DIR / "logs"
# else:
LOG_PATH = PROJECT_PATH / "logs"
FIN_LOG_PATH = EPS_LOGS  #  K_LOGS  # EPS_LOGS
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


# region selection dictionary
r"""
    Nino1-4 definitions are taken from:

    Trenberth, Kevin & National Center for Atmospheric Research Staff (Eds).
    Last modified 21 Jan 2020. "The Climate Data Guide: Nino SST Indices
    (Nino 1+2, 3, 3.4, 4; ONI and TNI)."
    Retrieved from
    https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni.

    Nino5 and Nino6 definitions taken from:

    @article{wang1999western,
    title={Western Pacific interannual variability
    associated with the El Ni{\~n}o-Southern Oscillation},
    author={Wang, Chunzai and Weisberg, Robert H and Virmani, Jyotika I},
    journal={Journal of Geophysical Research: Oceans},
    volume={104},
    number={C3},
    pages={5131--5149},
    year={1999},
    publisher={Wiley Online Library}
    }

    Nino X Index computation: (a) Compute area averaged total SST from Niño X
    region; (b) Compute monthly climatology (e.g., 1950-1979) for area averaged
    total SST from Niño X region, and subtract climatology from area averaged
    total SST time series to obtain anomalies; (c) Smooth the anomalies with a
    5-month running mean; (d) Normalize the smoothed values by its standard
    deviation over the climatological period.

    TNI computation: (a) Compute area averaged total SST from Niño 1+2 region;
    (b) Compute area averaged total SST from Niño 4 region; (c) Compute monthly
    climatologies (e.g., 1950-1979) for area averaged total SST from Niño 1+2
    region, and Niño 4 region, and subtract climatologies from area averaged
    total SST time series to obtain anomalies; (d) Normalize each time series of
    anomalies by their respective standard deviations over the climatological
    period; (e) Define the raw TNI as Niño 1+2 normalized anomalies minus Niño 4
    normalized anomalies; (f) Smooth the raw TNI with a 5-month running mean; (g)
    Normalize the smoothed TNI by its standard deviation over the climatological
    period.

    https://psl.noaa.gov/data/climateindices/list/

    https://psl.noaa.gov/data/correlation/tni.data

    # [west, east] [south, north] boundaries
    # X in degrees east, Y in degrees north.
    # colors compatible with matplotlib.
"""


SEL_DICT = {
    "pac": {"X": (100, 290), "Y": (-30, 30), "color": "#411900"},
    "nino1": {"X": (270, 280), "Y": (-10, -5), "color": "black"},
    "nino2": {"X": (270, 280), "Y": (-5, 0), "color": "black"},
    "nino1+2": {"X": (270, 280), "Y": (-10, 0), "color": "#0652ff"},
    "nino3": {"X": (210, 270), "Y": (-5, 5), "color": "#01a049"},
    "nino3.4": {"X": (190, 240), "Y": (-5, 5), "color": "#380282"},
    "nino4": {"X": (160, 210), "Y": (-5, 5), "color": "#fe01b1"},
    "nino5": {"X": (120, 140), "Y": (-5, 5), "color": "black"},
    "nino6": {"X": (140, 160), "Y": (8, 16), "color": "black"},
}
