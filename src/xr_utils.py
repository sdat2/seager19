"""Utilities around opening and processing netcdfs from this project."""
import pathlib
from typing import Union
import xarray as xr


def fix_calendar(
    xr_in: Union[xr.Dataset, xr.DataArray], timevar: str = "T"
) -> Union[xr.Dataset, xr.DataArray]:
    """Fix and decode the calendar.

    Args:
        xr_in (Union[xr.Dataset, xr.DataArray]): the xarray object input
        timevar (str, optional): The time variable name. Defaults to "T".

    Returns:
        Union[xr.Dataset, xr.DataArray]: same type and xr_in with fixed calendar.

    """
    # make all into dataset
    if isinstance(xr_in, xr.DataArray):
        ds = xr_in.to_dataset()
    else:
        ds = xr_in

    # add 360_day attribute
    if "calendar" not in ds[timevar].attrs:
        ds[timevar].attrs["calendar"] = "360_day"
    else:
        if ds[timevar].attrs["calendar"] == "360":
            ds[timevar].attrs["calendar"] = "360_day"

    # decode
    ds = xr.decode_cf(ds)

    # transform back into original type
    if isinstance(xr_in, xr.DataArray):
        xr_out = ds.to_array()
    else:
        xr_out = ds

    return xr_out


def can_coords(
    xr_obj: Union[xr.Dataset, xr.DataArray]
) -> Union[xr.Dataset, xr.DataArray]:
    """
    Transform an object into having the canonical coordinates if possible.

    Fail hard if impossible.

    TODO: Should fail hard at the moment, need to setup.

    Args:
        xr_obj (Union[xr.Dataset, xr.DataArray]): The dataset or datarray to
            canonicalise.

    Returns:
        Union[xr.Dataset, xr.DataArray]: The dataset that has been canoncilised.
            Function will raise an assertion error otherwise.
    """
    assert isinstance(xr_obj, Union[xr.DataArray, xr.Dataset])

    assert 1 < 2

    return xr_obj


def sel(
    xr_obj: Union[xr.Dataset, xr.DataArray], reg="pac"
) -> Union[xr.Dataset, xr.DataArray]:
    """
    Select a region of the dataset or datarray.

    Assumes
    reg options: "pac', 'nino3.4', "glob"
    https://www.ncdc.noaa.gov/teleconnections/enso/indicators/sst/

    Args:
        xr_obj (Union[xr.Dataset, xr.DataArray]): The xarray object.
            Needs to have canonical coordinates.
        reg (str, optional): The keyword region to select. Defaults to 'pac'.

    Returns:
        Union[xr.Dataset, xr.DataArray]: The downsized xarray object.

    Example:
        Effect example::

                if reg == "pac":
                    return xr_obj.sel(X=slice(100, 290), Y=slice(-30, 30))
                elif reg == "nino3.4":
                    return xr_obj.sel(X=slice(190, 240), Y=slice(-5, 5))
                elif reg == "nino4":
                    return xr_obj.sel(X=slice(160, 210), Y=slice(-5, 5))
                elif reg == "nino3":
                    return xr_obj.sel(X=slice(210, 270), Y=slice(-5, 5))
                elif reg == "nino1+2":
                    return xr_obj.sel(X=slice(270, 280), Y=slice(-10, 0))

    """

    sel_d = {
        "pac": {"X": (100, 200), "Y": (-30, 30)},
        "nino1+2": {"X": (270, 280), "Y": (-10, 0)},
        "nino3": {"X": (210, 270), "Y": (-5, 5)},
        "nino3.4": {"X": (190, 240), "Y": (-5, 5)},
        "nino4": {"X": (160, 210), "Y": (-5, 5)},
    }

    assert reg in sel_d

    sel_c = sel_d[reg]

    return xr_obj.sel(
        X=slice(sel_c["X"][0], sel_c["X"][1]), Y=slice(sel_c["Y"][0], sel_c["Y"][1])
    )


def open_dataset(path: Union[str, pathlib.Path]) -> xr.Dataset:
    """
    Open a dataset and format it.

    Args:
        path (Union[str, pathlib.Path]): the path to the netcdf dataset file.

    Returns:
        xr.Dataset: The formatted dataset.
    """
    return fix_calendar(xr.open_dataset(str(path), decode_times=False))


def open_dataarray(path: Union[str, pathlib.Path]) -> xr.DataArray:
    """
    Open a dataarray and format it.

    Args:
        path (Union[str, pathlib.Path]): the path to the netcdf datarray file.

    Returns:
        xr.DataArray: The formatted datarray.
    """
    return fix_calendar(xr.open_dataarray(str(path), decode_times=False))
