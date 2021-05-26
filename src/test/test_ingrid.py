"""Test the ingrid replacement python functions.

Example:
    Test using::
        pytest src/test/test_ingrid.py

"""
import os
import shutil
import xarray as xr
import recursive_diff
from src.data_loading.ingrid import linear_qflx_replacement
from src.constants import OCEAN_DATA_PATH, TEST_DIREC, OCEAN_OUTPUT_PATH
from src.data_loading.download import get_data
from src.models.model_setup import ModelSetup
from src.configs.load_config import load_config
from src.utils import delete_folder_contents


def test_ingrid() -> None:
    """Test the qflx replacement function."""

    # get_data if it does not exists
    get_data()
    cfg = load_config()
    print(cfg.name)

    delete_folder_contents(str(TEST_DIREC))
    setup = ModelSetup(str(TEST_DIREC))

    shutil.copy(str(OCEAN_OUTPUT_PATH / "om_diag.nc"), str(setup.ocean_output_path))

    # make a qflx-test file.
    linear_qflx_replacement(setup, output_file_name="qflx-test.nc")

    print("files in paths:\t\t\n", str(os.listdir(str(OCEAN_DATA_PATH))))

    # load different qflx files.
    qflx_test = xr.open_dataarray(
        str(os.path.join(setup.ocean_data_path, "qflx-test.nc")), decode_times=False
    )
    qflx_old = xr.open_dataarray(
        str(os.path.join(setup.ocean_data_path, "qflx-0.nc")), decode_times=False
    )

    # look at stuff
    print(qflx_test)
    print(qflx_old)

    # xr.testing.assert_equal(qflx_test, qflx_old)
    xr.testing.assert_allclose(qflx_test, qflx_old, atol=1e-2)

    for x in recursive_diff.recursive_diff(qflx_test, qflx_old, abs_tol=1e-2):
        print(x)

        if "-999.0" not in str(x):

            # otherwise it would catch errors about what
            # the missing value attribute is.
            assert False

    delete_folder_contents(str(TEST_DIREC))
