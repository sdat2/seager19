"""Test xr_utils."""
import xarray as xr
import recursive_diff
from src.xr_utils import fix_calendar, open_dataset, sel, cut_and_taper, can_coords
from src.data_loading.download import get_data
from src.constants import OCEAN_DATA_PATH


def test_fix_calendar() -> None:
    """Test `src.xr_utils.fix_calendar`."""
    get_data()
    ds_new: xr.Dataset = fix_calendar(open_dataset(OCEAN_DATA_PATH / "qflx.nc"))
    da_new: xr.DataArray = fix_calendar(open_dataset(OCEAN_DATA_PATH / "qflx.nc").qflx)
    print(ds_new)
    print(da_new)
    sel(da_new)
    sel(da_new, reg="nino3.4")

    recursive_diff.recursive_diff(
        cut_and_taper(da_new.isel(Z=0, T=0, variable=0)),
        cut_and_taper(can_coords(da_new.isel(Z=0, T=0)).isel(variable=0)),
    )
