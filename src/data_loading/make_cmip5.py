"""Make CMIP5 variables to put into the atmosphere model."""
import numpy as np
import xarray as xr
from src.constants import MMM_V23_HIST, MMM_V23_RCP85, ATMOS_DATA_PATH, U_HIST, V_HIST
from src.xr_utils import open_dataset, open_dataarray, can_coords
from src.utils import timeit
from src.data_loading.download import get_uv, get_mmm
from src.visualisation.comp import return_figure_ds


@timeit
def qair2rh(qair: xr.DataArray, temp: xr.DataArray, pres: xr.DataArray) -> xr.DataArray:
    """
    get the relative humdity from the specific humidity.

    Args:
        qair (xr.DataArray): The specific humidity (dimensionless).
        temp (xr.DataArray): The temperature (kelvin).
        pres (xr.DataArray): The pressure (pascal).

    Returns:
        xr.DataArray: The relative humidity.
    """
    t_0c = 273.15
    es = 6.112 * np.exp((17.76 * (temp - t_0c)) / (temp - t_0c + 243.5))
    e = qair * pres / (0.378 * qair + 0.622)
    rh = e / es
    # rh[rh > 100] = 1
    # rh[rh < 0] = 0
    return rh


@timeit
def q2rh(qair: xr.DataArray, temp: xr.DataArray, pres: xr.DataArray) -> xr.DataArray:
    """
    get the relative humdity from the specific humidity.

    Args:
        qair (xr.DataArray): The specific humidity (dimensionless).
        temp (xr.DataArray): The temperature (kelvin)
        pres (xr.DataArray): The pressure (pascal).

    Returns:
        xr.DataArray: The relative humidity.
    """
    return 0.263 * pres * qair / (np.exp(17.67 * (temp - 273.16) / (temp - 29.65)))


@timeit
def make_rh() -> None:
    """
    Make the cmip5 relative humidity mean for 1956-2016.
    """
    get_mmm()
    ps = xr.concat(
        [open_dataset(MMM_V23_HIST)["ps"], open_dataset(MMM_V23_RCP85)["ps"]], "T"
    )
    ts = xr.concat(
        [open_dataset(MMM_V23_HIST)["ts"], open_dataset(MMM_V23_RCP85)["ts"]], "T"
    )
    qa = xr.concat(
        [open_dataset(MMM_V23_HIST)["huss"], open_dataset(MMM_V23_RCP85)["huss"]], "T"
    )
    rh = qair2rh(qa, ts, ps)
    rh = q2rh(qa, ts, ps)
    rh_e = open_dataarray(ATMOS_DATA_PATH / "rh-ECMWF-clim60.nc")
    rh_mean = (
        rh.sel(T=slice("1956", "2016"))
        .mean("T")
        .rename("rh")
        .interp_like(rh_e)
        .ffill("X")
        .bfill("X")
    )
    rh_e_ll = xr.open_dataarray(ATMOS_DATA_PATH / "rh-ECMWF-clim60.nc")
    rh_new = rh_e_ll.copy()
    rh_new[:, :] = rh_mean[:, :]
    rh_new.to_netcdf(ATMOS_DATA_PATH / "rh-CMIP5-clim60.nc")


@timeit
def make_sfcwind() -> None:
    """
    Make the cmip5 surface wind mean.
    """
    get_uv()
    us = (
        open_dataarray(U_HIST)
        .bfill("plev")
        .isel(plev=0, variable=0)
        .drop("variable")
        .drop("plev")
    )
    vs = (
        open_dataarray(V_HIST)
        .bfill("plev")
        .isel(plev=0, variable=0)
        .drop("variable")
        .drop("plev")
    )
    windsp = np.sqrt(np.square(us) + np.square(vs))
    swd_e = (
        open_dataarray(ATMOS_DATA_PATH / "sfcWind-ECMWF-clim.nc")
        .isel(variable=0)
        .drop("variable")
    )
    swd_e_ll = xr.open_dataarray(ATMOS_DATA_PATH / "sfcWind-ECMWF-clim.nc")
    wsp_mean = (
        windsp.mean("T").interp_like(swd_e).bfill("X").ffill("X").bfill("Y").ffill("Y")
    )
    wsp_new = swd_e_ll.copy()
    wsp_new[:, :] = wsp_mean[:, :]
    wsp_new.to_netcdf(ATMOS_DATA_PATH / "sfcWind-CMIP5-clim.nc")


if __name__ == "__main__":
    # python3 src/data_loading/make_cmip5.py
    print(return_figure_ds("3"))
    print(return_figure_ds("5c"))
    print(return_figure_ds("1d"))
    print(return_figure_ds("1e"))
    print(return_figure_ds("1f"))
    print(return_figure_ds("5c"))
    wsp = can_coords(return_figure_ds("5c")["wnspClim"])
    print(wsp)
    rh = can_coords(return_figure_ds("5c")["rh"])
    print(rh)
    print(xr.open_dataarray(ATMOS_DATA_PATH / "sfcWind-ECMWF-clim.nc"))
    print(xr.open_dataarray(ATMOS_DATA_PATH / "rh-ECMWF-clim60.nc"))
    # make_rh()
    # make_sfcwind()

# 1 - forced ocean
# d - rising CO2, observed winds
# e - rising CO2, fixed winds
# f - fixed CO2, rising winds
# 2 - forced atmosphere
# c - no heating over land, ECMWF forcing
# d - heating over land, ECMWF forcing
# 3 - coupled atmosphere ocean model
#     heating over land, ECMWF inputs
# a - sst change
# b - prcp, utrend, vtrend change
# 4 - Trend in the thermocline
# a - ORAS4 model
# b - forced with ORAS4 winds
# c - forced with atmosphere-ocean
# 5 - coupled models
# a - CM-ECMWF world
# b - CM-ECMWF C-RH - change the relative humidity to CMIP5
# c - CM-ECMWF C-RH W - chang the
# d - CM-CMIP5 world qfluxed towards the CMIP5 mmm
# f - unknown thermocline thing.
# g - ECMWF graph
# h - CMIP5 multimodel mean
