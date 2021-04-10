"""Download the data from dropbox links

# import shutil
# from typing import List

from src.data_loading import get_data

"""
import os
import requests
import zipfile
from tqdm import tqdm
from src.utils import timeit
from src.constants import (
    OCEAN_DATA_PATH,
    OCEAN_OUTPUT_PATH,
    ATMOS_TMP_PATH,
    ATMOS_DATA_PATH,
    DROP_ATMOS_DATA,
    DROP_ATMOS_TMP,
    DROP_OCEAN_OUTPUT,
    DROP_OCEAN_DATA,
)


@timeit
def get_and_unzip(direc: str, url: str, name: str) -> None:
    """Get the data and unzip it.

    Args:
        direc (str): directory to put the data in
        url (str): url of the zip file
        name (str): name of file
    """

    @timeit
    def get_zip() -> None:
        write_path = os.path.join(direc, name)
        req = requests.get(url, stream=True)
        with open(write_path, "wb") as file:
            for chunk in tqdm(req.iter_content(chunk_size=128)):
                file.write(chunk)

    @timeit
    def un_zip() -> None:
        write_path = os.path.join(direc, name)
        with zipfile.ZipFile(write_path, "r") as zip_ref:
            zip_ref.extractall(direc)

    get_zip()
    un_zip()


@timeit
def get_data(atmos: bool = True, ocean: bool = True):
    """Download the relevant dataset from a Dropbox link and extract it."""

    if atmos and ocean:
        print(OCEAN_DATA_PATH, OCEAN_OUTPUT_PATH, ATMOS_TMP_PATH, ATMOS_DATA_PATH)
        print(DROP_ATMOS_DATA, DROP_ATMOS_TMP, DROP_OCEAN_OUTPUT, DROP_OCEAN_DATA)

    lol = [
        [
            "atmos",
            [
                ["https://www.dropbox.com/s/j7x3bjfnb8fdw3b/tmp.zip?raw=1", "tmp.zip"],
                [
                    "https://www.dropbox.com/s/r82fv6g8sdwykzj/DATA.zip?raw=1",
                    "DATA.zip",
                ],
            ],
        ],
        [
            "ocean",
            [
                [
                    "https://www.dropbox.com/s/luoo402lzbsyi5r/DATA.zip?raw=1",
                    "DATA.zip",
                ],
                [
                    "https://www.dropbox.com/s/q2qxeucdsb9958j/output.zip?raw=1",
                    "output.zip",
                ],
            ],
        ],
    ]

    for item in lol:
        direc = item[0]
        if not os.path.exists(direc):
            os.mkdir(direc)
        for pair in item[1]:
            url = pair[0]
            name = pair[1]
            get_and_unzip(direc, url, name)


if __name__ == "__main__":
    get_data()
