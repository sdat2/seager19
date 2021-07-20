"""Download the data from dropbox links.

Example:
    Import statement::

        from src.data_loading import get_data

"""
import os
import shutil
import requests
import zipfile
from tqdm import tqdm
from src.utils import timeit
from src.constants import (
    OCEAN_PATH,
    ATMOS_PATH,
    DATA_PATH,
    FIGURE_DATA_NAME,
    NOAA_DATA_NAME,
    NINO3_4_TEST_NAME,
    NINO3_4_TEST_CODE,
)


@timeit
def get_and_unzip(direc: str, url: str, name: str) -> None:
    """Get the data and unzip it.

    Args:
        direc (str): directory to put the data in.
        url (str): url of the zip file.
        name (str): name of file.

    """

    write_path = os.path.join(direc, name)

    @timeit
    def get_zip() -> None:
        req = requests.get(url, stream=True)
        with open(write_path, "wb") as file:
            for chunk in tqdm(req.iter_content(chunk_size=128)):
                file.write(chunk)

    @timeit
    def un_zip() -> None:
        with zipfile.ZipFile(write_path, "r") as zip_ref:
            zip_ref.extractall(direc)

    @timeit
    def clean_up() -> None:
        os.remove(write_path)
        mac_path = os.path.join(direc, "__MACOSX")
        if os.path.exists(mac_path):
            shutil.rmtree(mac_path)

    get_zip()
    un_zip()
    clean_up()


PREFIX = "https://www.dropbox.com/s/"
SUFFIX = "?raw=1"


@timeit
def get_data() -> None:
    """Download the relevant dataset from a Dropbox link and extract it."""

    lol = [
        [
            str(ATMOS_PATH),
            [
                [PREFIX + "j7x3bjfnb8fdw3b/tmp.zip" + SUFFIX, "tmp.zip"],
                [
                    PREFIX + "r82fv6g8sdwykzj/DATA.zip" + SUFFIX,
                    "DATA.zip",
                ],
            ],
        ],
        [
            str(OCEAN_PATH),
            [
                [
                    PREFIX + "luoo402lzbsyi5r/DATA.zip" + SUFFIX,
                    "DATA.zip",
                ],
                [
                    PREFIX + "q2qxeucdsb9958j/output.zip" + SUFFIX,
                    "output.zip",
                ],
            ],
        ],
    ]
    _get_data(lol)


def get_figure_data(force_refresh: bool = False) -> None:
    """Downloads figure nc."""

    code = "r1vp6ny8wovyq2a/"
    name = FIGURE_DATA_NAME + ".zip"

    lol = [
        [
            str(DATA_PATH),
            [
                [
                    PREFIX + code + name + SUFFIX,
                    name,
                ],
            ],
        ],
    ]

    _get_data(lol, force_refresh=force_refresh)


def get_noaa_data() -> None:
    """Downloads noaa sst nc."""

    code = "f7sdvmvsziaxxj9/"
    name = NOAA_DATA_NAME + ".zip"

    lol = [
        [
            str(DATA_PATH),
            [
                [
                    PREFIX + code + name + SUFFIX,
                    name,
                ],
            ],
        ],
    ]

    _get_data(lol)


def get_test_nino_data() -> None:
    """Downloads noaa sst nc."""

    code = NINO3_4_TEST_CODE
    name = NINO3_4_TEST_NAME + ".zip"

    lol = [
        [
            str(DATA_PATH),
            [
                [
                    PREFIX + code + name + SUFFIX,
                    name,
                ],
            ],
        ],
    ]

    _get_data(lol)


def get_original_models() -> None:
    """Downloads original models.

    # https://www.dropbox.com/s/u5ufak4kv8peajd/atmos-model.zip?dl=0
    # https://www.dropbox.com/s/s6u80p1s1e5qom4/ocean-model.zip?dl=0
    """

    for code, name in [
        ["u5ufak4kv8peajd/", "atmos-model.zip"],
        ["s6u80p1s1e5qom4/", "ocean-model.zip"],
    ]:

        lol = [
            [
                str(DATA_PATH),
                [
                    [
                        PREFIX + code + name + SUFFIX,
                        name,
                    ],
                ],
            ],
        ]

        _get_data(lol)


def get_member_data(force_refresh: bool = False) -> None:
    """Downloads ensemble members."""
    lol = [
        [
            str(DATA_PATH),
            [
                [
                    PREFIX + "udui5x9c3q7y2ca/ts_nc.zip" + SUFFIX,
                    "ts_nc.zip",
                ],
            ],
        ],
    ]

    _get_data(lol, force_refresh=force_refresh)


def get_mmm(force_refresh: bool = False) -> None:
    """Downloads multimodel mean."""
    # https://www.dropbox.com/s/r2k6qwaijq5a8qv/mmm-v2.3-full-rep.zip?dl=0
    lol = [
        [
            str(DATA_PATH),
            [
                [
                    PREFIX + "r2k6qwaijq5a8qv/mmm-v2.3-full-rep.zip" + SUFFIX,
                    "mmm-v2.3-full-rep.zip",
                ],
            ],
        ],
    ]

    _get_data(lol, force_refresh=force_refresh)


# https://www.dropbox.com/s/18oca0kaft5tdy0/cmip6-mmm.zip?dl=0


def get_uv(force_refresh: bool = False) -> None:
    """Downloads uv mean."""
    # https://www.dropbox.com/s/9zuxdyhapwr1eat/hist-winds-cmip5.zip?dl=0
    lol = [
        [
            str(DATA_PATH),
            [
                [
                    PREFIX + "9zuxdyhapwr1eat/hist-winds-cmip5.zip" + SUFFIX,
                    "hist-winds-cmip5.zip",
                ],
            ],
        ],
    ]

    _get_data(lol, force_refresh=force_refresh)


def get_cmip6_mmm(force_refresh: bool = False) -> None:
    """Downloads cmip6 means."""
    lol = [
        [
            str(DATA_PATH),
            [
                [
                    PREFIX + "18oca0kaft5tdy0/cmip6-mmm.zip" + SUFFIX,
                    "cmip6-mmm.zip",
                ],
            ],
        ],
    ]

    _get_data(lol, force_refresh=force_refresh)


# https://www.dropbox.com/s/4fsdixia0dl9gt4/nc80.zip?dl=0
def get_cmip6_mmm_clim60(force_refresh: bool = False) -> None:
    """Downloads cmip6 means."""
    lol = [
        [
            str(DATA_PATH),
            [
                [
                    PREFIX + "c22xgus4wk6jzlv/nc_mean.zip" + SUFFIX,
                    "nc_mean.zip",
                ],
            ],
        ],
    ]

    _get_data(lol, force_refresh=force_refresh)


def get_ts_ensemble(force_refresh: bool = False) -> None:
    """Downloads cmip6 ensemble for ts."""
    lol = [
        [
            str(DATA_PATH),
            [
                [
                    PREFIX + "4fsdixia0dl9gt4/nc80.zip" + SUFFIX,
                    "nc80.zip",
                ],
            ],
        ],
    ]

    _get_data(lol, force_refresh=force_refresh)


def _get_data(lol: list, force_refresh: bool = False) -> None:
    """Gets the data using lol."""
    for item in lol:
        direc = item[0]
        if not os.path.exists(direc):
            os.mkdir(direc)
        for pair in item[1]:
            url = pair[0]
            name = pair[1]
            full_direc = str(os.path.join(direc, os.path.splitext(pair[1])[0]))
            if not os.path.exists(full_direc):
                print("Dowloading " + full_direc)
                get_and_unzip(direc, url, name)
                print(full_direc + " created.")
            else:
                if not force_refresh:
                    print(full_direc + " already exists, not going to redownload.")
                else:
                    print(full_direc + " already exists, but about to redownload.")
                    get_and_unzip(direc, url, name)


if __name__ == "__main__":
    # python3 src/data_loading/download.py
    # get_ts_ensemble()
    # get_data()
    # get_figure_data()
    # get_member_data()
    # get_original_models()
    # get_mmm()
    # get_uv()
    # get_cmip6_mmm()
    get_cmip6_mmm_clim60()
