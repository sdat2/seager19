"""Make input variable fields to put into the atmosphere and ocean model."""
import os
import numpy as np
import xarray as xr
from src.constants import (
    MMM_V23_HIST,
    MMM_V23_RCP85,
    ATMOS_DATA_PATH,
    U_HIST,
    V_HIST,
    CMIP6_CLIM60_PATH,
    MODEL_NAMES,
    cmip6_mmm_mean,
    atmos_input_file_path,
    cmip6_mmm_trend,
)
from src.xr_utils import open_dataset, open_dataarray, can_coords
from src.utils import timeit
from src.data_loading.download import get_uv, get_mmm, get_figure_data
from src.visualisation.comp_v_seager19 import return_figure_ds


@timeit
def qair2rh(qair: xr.DataArray, temp: xr.DataArray, pres: xr.DataArray) -> xr.DataArray:
    """
    Get the relative humdity from the specific humidity.

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


@timeit
def get_sfcwind() -> None:
    """Get the CMIP5 surface wind from the figure data."""
    get_figure_data()
    wsp = can_coords(return_figure_ds("5d")["wnspClim"])
    swd_e60 = (
        open_dataarray(ATMOS_DATA_PATH / "sfcWind-ECMWF-clim60.nc")
        .isel(variable=0)
        .drop("variable")
    )
    swd_e = (
        open_dataarray(ATMOS_DATA_PATH / "sfcWind-ECMWF-clim.nc")
        .isel(variable=0)
        .drop("variable")
    )
    swd_e_ll = xr.open_dataarray(ATMOS_DATA_PATH / "sfcWind-ECMWF-clim.nc")
    swd_e_ll60 = xr.open_dataarray(ATMOS_DATA_PATH / "sfcWind-ECMWF-clim60.nc")
    wsp_clim = wsp.interp_like(swd_e).bfill("X").ffill("X").bfill("Y").ffill("Y")
    wsp_new = swd_e_ll.copy()
    wsp_new[:, :] = wsp_clim[:, :]
    wsp_new.to_netcdf(ATMOS_DATA_PATH / "sfcWind-CMIP5-clim.nc")
    wsp_clim60 = wsp.interp_like(swd_e60).bfill("X").ffill("X").bfill("Y").ffill("Y")
    wsp_new60 = swd_e_ll60.copy()
    wsp_new60[:, :] = wsp_clim60[:, :]
    wsp_new60.to_netcdf(ATMOS_DATA_PATH / "sfcWind-CMIP5-clim60.nc")


@timeit
def get_sfcwind_6() -> None:
    """Get the CMIP6 surface wind and put it in the right files."""
    # wsp = can_coords(return_figure_ds("5d")["wnspClim"])
    wsp = can_coords(xr.open_dataarray(os.path.join(CMIP6_CLIM60_PATH, "wsp.nc")))
    print(wsp.max())
    print(wsp.min())
    wsp = wsp.where(wsp < 20).fillna(20)
    # wsp = wsp.where(wsp > 2).fillna(2)
    print(wsp.max())
    # print(wsp.min())

    swd_e60 = (
        open_dataarray(ATMOS_DATA_PATH / "sfcWind-ECMWF-clim60.nc")
        .isel(variable=0)
        .drop("variable")
    )
    swd_e = (
        open_dataarray(ATMOS_DATA_PATH / "sfcWind-ECMWF-clim.nc")
        .isel(variable=0)
        .drop("variable")
    )
    swd_e_ll = xr.open_dataarray(ATMOS_DATA_PATH / "sfcWind-ECMWF-clim.nc")
    swd_e_ll60 = xr.open_dataarray(ATMOS_DATA_PATH / "sfcWind-ECMWF-clim60.nc")
    wsp_clim = wsp.interp_like(swd_e).bfill("X").ffill("X").bfill("Y").ffill("Y")
    wsp_new = swd_e_ll.copy()
    wsp_new[:, :] = wsp_clim[:, :]
    wsp_new.to_netcdf(ATMOS_DATA_PATH / "sfcWind-CMIP6-clim.nc")
    wsp_clim60 = wsp.interp_like(swd_e60).bfill("X").ffill("X").bfill("Y").ffill("Y")
    wsp_new60 = swd_e_ll60.copy()
    wsp_new60[:, :] = wsp_clim60[:, :]
    wsp_new60.to_netcdf(ATMOS_DATA_PATH / "sfcWind-CMIP6-clim60.nc")


@timeit
def get_rh() -> None:
    """Get the CMIP5 relative humidity from the figure data."""
    get_figure_data()
    rh = can_coords(return_figure_ds("5c")["rh"])
    rh_e = open_dataarray(ATMOS_DATA_PATH / "rh-ECMWF-clim60.nc")
    rh_clim = rh.interp_like(rh_e).bfill("X").ffill("X").bfill("Y").ffill("Y")
    rh_e_ll = xr.open_dataarray(ATMOS_DATA_PATH / "rh-ECMWF-clim60.nc")
    rh_new = rh_e_ll.copy()
    rh_new[:, :] = rh_clim[:, :]
    rh_new.to_netcdf(ATMOS_DATA_PATH / "rh-CMIP5-clim60.nc")


@timeit
def get_rh_6() -> None:
    """Get the CMIP6 relative humidity from the figure data."""
    rh = can_coords(xr.open_dataarray(os.path.join(CMIP6_CLIM60_PATH, "hur.nc")))
    rh = rh.where(rh < 100).fillna(100)
    rh = rh.where(rh > 20).fillna(20)
    print(rh.min())
    print(rh.max())
    rh_e = open_dataarray(ATMOS_DATA_PATH / "rh-ECMWF-clim60.nc")
    rh_clim = rh.interp_like(rh_e).bfill("X").ffill("X").bfill("Y").ffill("Y")
    rh_e_ll = xr.open_dataarray(ATMOS_DATA_PATH / "rh-ECMWF-clim60.nc")
    rh_new = rh_e_ll.copy()
    rh_new[:, :] = rh_clim[:, :]
    rh_new.to_netcdf(ATMOS_DATA_PATH / "rh-CMIP6-clim60.nc")


coord_rename_dict = {
    "clim60": {"Y": "lat", "X": "lon"},
    # "clim": {"Y": "X", "X": "X"},
    # "trend": {"Y": "X", "X": "X"},
}
sel_dict = {
    "clim60": {"Y": slice(-60, 60)},
    "clim": {
        "Y": slice(-90, 90),
    },
    "trend": {
        "Y": slice(-90, 90),
    },
}
cmip6_mmm_func_d = {
    "clim60": cmip6_mmm_mean,
    "clim": cmip6_mmm_mean,
    "trend": cmip6_mmm_trend,
}
var_rename_dict = {"rh": "hurs", "ps": "ps"}  # "psl"}
var_regrid_list = ["ps"]


def kgms_mmday() -> float:
    """Bencini, Roberto. (2016). Re: How do I convert ERA Interim precipitation estimates from kg/m2/s to mm/day?. Retrieved from: https://www.researchgate.net/post/How-do-I-convert-ERA-Interim-precipitation-estimates-from-kg-m2-s-to-mm-day/56d71a54b0366dd61e0de055/citation/download."""
    return 86400


def generate(var: str, model: str = "S", ending: str = "clim60"):
    """
    Generate an input variable.

    Args:
        var (str): Variable string. e.g. "ts".
        model (str, optional): Model character. Defaults to "S".
        ending (str, optional): Input file ending. Defaults to "clim60".
    """
    if var in var_rename_dict:
        old_var = var_rename_dict[var]
        cmip6_mean = xr.open_dataarray(cmip6_mmm_func_d[ending](old_var))
        cmip6_mean = cmip6_mean.rename(var)
    else:
        cmip6_mean = xr.open_dataarray(cmip6_mmm_func_d[ending](var))
    ecmwf_mean = xr.open_dataarray(atmos_input_file_path(var=var, ending=ending))
    cmip6_mean = cmip6_mean.sel(**sel_dict[ending])
    print(var, model, ending)
    print(ecmwf_mean.attrs["units"], "\t\t\t", cmip6_mean.attrs["units"])

    # if (
    #    var == "ts" and ending == "clim"
    # ):  # ecmwf_mean.attrs["units"] == "degree_Celsius":
    # cmip6_mean = cmip6_mean - 273.15
    if var == "ps":
        cmip6_mean = cmip6_mean / 100
        cmip6_mean.attrs["units"] = ecmwf_mean.attrs["units"]
    #if var == "pr":
    #    #cmip6_mean = cmip6_mean  # * kgms_mmday()
    #    #cmip6_mean.attrs["units"] = ecmwf_mean.attrs["units"]
    if ending in coord_rename_dict:
        cmip6_mean = cmip6_mean.rename(coord_rename_dict[ending])
    ecmwf_dims = ecmwf_mean.dims
    if var in var_regrid_list:
        # this will break if clim60 requires regridding.
        cmip6_mean = (
            cmip6_mean.interp_like(ecmwf_mean)
            .bfill(ecmwf_dims[0])
            .ffill(ecmwf_dims[0])
            .bfill(ecmwf_dims[1])
            .ffill(ecmwf_dims[1])
        )
    assert cmip6_mean.dims == ecmwf_mean.dims
    print(
        ecmwf_mean.mean(ecmwf_dims[1]).mean(ecmwf_dims[0]).values,
        "\t\t\t",
        cmip6_mean.mean(ecmwf_dims[1]).mean(ecmwf_dims[0]).values,
    )
    print("========================================================")
    new_mean = ecmwf_mean.copy()
    new_mean[:, :359] = cmip6_mean[:, :359]
    new_mean.attrs["center"] = MODEL_NAMES[model] + " multi model mean"
    new_mean.to_netcdf(atmos_input_file_path(var=var, ending=ending, model=model))


def generate_all() -> None:
    """
    Generate all cmip6 entries.
    """
    ending_d = {
        "pr": ["clim", "trend"],
        "rh": ["clim60"],
        "sfcWind": ["clim", "clim60"],
        "ts": ["clim", "clim60", "trend"],
        "ps": ["clim"],
        "clt": ["clim60"],
    }
    for var in ending_d:
        for ending in ending_d[var]:
            generate(var, ending=ending)


if __name__ == "__main__":
    # python3 src/data_loading/make_inputs.py
    # print(return_figure_ds("5f"))
    # make_rh()
    # make_sfcwind()
    # get_figure_data()
    # get_rh()  # get it from the figure data
    # get_sfcwind()  # get it from the figure data
    # get_rh_6()
    # get_sfcwind_6()
    generate_all()


# pylint: disable=pointless-string-statement
"""
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
"""
