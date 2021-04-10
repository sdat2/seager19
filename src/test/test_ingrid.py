"""Test the ingrid replacement python functions."""
import os
import xarray as xr
import recursive_diff
from src.data_loading.ingrid import linear_qflx_replacement
from src.constants import OCEAN_DATA_PATH
from src.data_loading.download import get_data


def test_ingrid() -> None:
    """Test the qflx replacement function."""

    # get_data if it does not exists
    get_data()

    # make a qflx-test file.
    linear_qflx_replacement(output_file_name="qflx-test.nc")

    print("files in paths:\t\t\n", str(os.listdir(str(OCEAN_DATA_PATH))))

    # load different qflx files.
    qflx_test = xr.open_dataarray(
        str(OCEAN_DATA_PATH / "qflx-test.nc"), decode_times=False
    )
    qflx_old = xr.open_dataarray(str(OCEAN_DATA_PATH / "qflx-0.nc"), decode_times=False)

    # look at stuff
    print(qflx_test)
    print(qflx_old)

    # xr.testing.assert_equal(qflx_test, qflx_old)
    xr.testing.assert_allclose(qflx_test, qflx_old, atol=1e-2)

    for x in recursive_diff.recursive_diff(qflx_test, qflx_old, abs_tol=1e-2):
        print(x)

        if "-999.0" not in str(x):

            # otherwise it would catch errors about what the managed is.
            assert False
