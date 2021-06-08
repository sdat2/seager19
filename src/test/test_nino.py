"""Test the nino functions."""
import xarray as xr
import recursive_diff
from src.data_loading.download import get_test_nino_data
from src.constants import NINO3_4_TEST_PATH
from src.metrics import calculate_nino3_4_from_noaa


def test_nino() -> None:
    """Check if the nino function produces the same data as normal."""
    get_test_nino_data()
    recursive_diff.recursive_diff(
        calculate_nino3_4_from_noaa(),
        xr.open_dataarray(str(NINO3_4_TEST_PATH)),
    )
