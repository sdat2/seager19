"""Test the ingrid replacement python functions.

from src.data_loading import get_data
"""
import xarray as xr
import recursive_diff
from src.data_loading.ingrid import linear_qflx_replacement
from src.constants import OCEAN_DATA_PATH
from src.data_loading import get_data


def test_ingrid() -> None:
    """Test the qflx replacement function."""
    
    get_data()

    linear_qflx_replacement(output_file_name="qflx-test.nc")
    qflx_test = xr.open_dataarray(OCEAN_DATA_PATH / "qflx-test.nc", decode_times=False)
    qflx_old = xr.open_dataarray(OCEAN_DATA_PATH / "qflx-backup.nc", decode_times=False)

    print(qflx_test)

    print(qflx_old)

    # xr.testing.assert_equal(qflx_test, qflx_old)
    xr.testing.assert_allclose(qflx_test, qflx_old, atol=1e-2)

    for x in recursive_diff.recursive_diff(qflx_test, qflx_old, abs_tol=1e-2):
        print(x)
        assert False
