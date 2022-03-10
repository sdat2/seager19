"""Module to process CMIP6 ensemble in gws."""
import os
import xarray as xr
from src.constants import (
    cmip6_ensemble_var,
    CMIP6_ENSEMBLE_TRENDS,
    CMIP6_ENSEMBLE_CLIMATOLOGIES,
    CMIP6_ENSEMBLE_MEANS,
)
from src.utils import timeit
from src.xr_utils import can_coords, get_trend, get_clim


@timeit
def process_var(var: str = "ts") -> None:
    """Process an ensemble of variables into key netcdfs."""
    direc = cmip6_ensemble_var(var)
    member_list = [x[len(var) + 1 : -6] for x in os.listdir(direc)]
    file_list = [os.path.join(direc, x) for x in os.listdir(direc)]

    # Add member dimension so that it can be loaded in

    for i, member in enumerate(member_list):
        ds = xr.load_dataset(file_list[i])
        if "member" not in ds.dims:
            print(i, member)
            ds = ds.expand_dims({"member": [member_list[i]]})
            ds.to_netcdf(file_list[i], mode="w")
            del ds

    # load the data into dask.
    da = can_coords(xr.open_mfdataset(file_list))[var].sel(T=slice("1958", "2017"))

    # calculate key fields from period
    climatology = get_clim(da).compute()
    mean = da.mean("T").compute()
    trend = get_trend(da).compute()

    # remove .80 from member id for consistency if it's there.
    trend = trend.assign_coords(
        {"member": [x[:-3] for x in trend.member.values if x.endswith(".80")]}
    )
    mean = mean.assign_coords(
        {"member": [x[:-3] for x in mean.member.values if x.endswith(".80")]}
    )
    climatology = climatology.assign_coords(
        {"member": [x[:-3] for x in climatology.member.values if x.endswith(".80")]}
    )

    # save key fields
    climatology.to_netcdf(os.path.join(CMIP6_ENSEMBLE_CLIMATOLOGIES, var + ".nc"))
    mean.to_netcdf(os.path.join(CMIP6_ENSEMBLE_MEANS, var + ".nc"))
    trend.to_netcdf(os.path.join(CMIP6_ENSEMBLE_TRENDS, var + ".nc"))


if __name__ == "__main__":
    # python src/data_loading/ensemble_processing.py
    for var_temp in ["hur", "pr", "sfcWind"]:
        process_var(var_temp)
