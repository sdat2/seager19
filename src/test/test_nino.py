"""Test the nino functions."""
import xarray as xr
import recursive_diff
from src.data_loading.download import get_test_nino_data, get_noaa_data
from src.constants import NINO3_4_TEST_PATH
from src.metrics import (
    calculate_nino3_4_from_noaa,
    replace_nino3_4_from_noaa,
)
from src.xr_utils import get_trend


def test_nino() -> None:
    """Check if the nino function produces the same data as normal."""
    get_noaa_data()
    get_test_nino_data()
    metric, clim = calculate_nino3_4_from_noaa()
    for x in recursive_diff.recursive_diff(
        metric,
        xr.open_dataarray(str(NINO3_4_TEST_PATH)),
    ):
        print(x)
        assert False
    print(clim)
    replace_nino3_4_from_noaa()
    get_trend(metric)
