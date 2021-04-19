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
from src.constants import OCEAN_PATH, ATMOS_PATH, PROJECT_PATH


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


def get_figure_data() -> None:
    """Downloads figure nc."""

    code = "r1vp6ny8wovyq2a/"
    name = "Seager_etal_NCC-2019_datasetdatafiles.nc.zip"

    lol = [
        [
            str(PROJECT_PATH),
            [
                [
                    PREFIX + code + name + SUFFIX,
                    name,
                ],
            ],
        ],
    ]

    _get_data(lol)


def get_member_data() -> None:
    """Downloads ensemble members."""
    lol = [
        [
            str(OCEAN_PATH),
            [
                [
                    PREFIX + "udui5x9c3q7y2ca/ts_nc.zip" + SUFFIX,
                    "ts_nc.zip",
                ],
            ],
        ],
    ]

    _get_data(lol)


def _get_data(lol: list) -> None:
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
                #  or "ocean/DATA" in full_direc:
                print("Dowloading " + full_direc)
                get_and_unzip(direc, url, name)
                print(full_direc + " created.")
            else:
                print(full_direc + " already exists, not going to redownload.")


# str(OCEAN_PATH)
# "https://www.dropbox.com/s/udui5x9c3q7y2ca/ts_nc.zip?raw=1"
# "tc_nc"

if __name__ == "__main__":
    # python3 src/data_loading/download.py
    get_data()
    get_figure_data()
    get_member_data()
