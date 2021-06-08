"""Test the nino functions."""
import xarray as xr
import recursive_diff
from src.data_loading.download import get_test_nino_data
from src.constants import NINO3_4_TEST_PATH
from src.metrics import calculate_nino3_4_from_noaa


def test_nino() -> None:
    """Check if the nino function produces the same data as normal."""
    get_test_nino_data()
    metric, mean_state = calculate_nino3_4_from_noaa()
    del metric.attrs["mean_state"]
    del metric.attrs["reg"]
    for x in recursive_diff.recursive_diff(
        metric,
        xr.open_dataarray(str(NINO3_4_TEST_PATH)),
    ):
        print(x)
        assert False
    print(mean_state)
