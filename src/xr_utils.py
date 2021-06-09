"""Utilities around opening and processing netcdfs from this project."""
import numpy as np
import pathlib
from typing import Union, Tuple
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

    # likely will not fail.

    t_list = ["T_0" + str(x) for x in range(5)]
    t_list.append("T")
    t_list.append("time")
    t_list.append(timevar)

    for t_dim in t_list:

        if t_dim in ds.dims:
            # add 360_day attribute
            if "calendar" not in ds[t_dim].attrs:
                ds[t_dim].attrs["calendar"] = "360_day"
            else:
                if ds[t_dim].attrs["calendar"] == "360":
                    ds[t_dim].attrs["calendar"] = "360_day"

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

    Args:
        xr_obj (Union[xr.Dataset, xr.DataArray]): The dataset or datarray to
            canonicalise.

    Returns:
        Union[xr.Dataset, xr.DataArray]: The dataset that has been canoncilised.
            Function will raise an assertion error otherwise.
    """
    assert isinstance(xr_obj, (xr.DataArray, xr.Dataset))

    def only1(l):
        # check that there is only one true value in a list.
        # answer taken from:
        # https://stackoverflow.com/questions/16801322/
        # how-can-i-check-that-a-list-has-one-and-only-one-truthy-value
        true_found = False
        for v in l:
            if v:
                # a True was found!
                if true_found:
                    # found too many True's
                    return False
                else:
                    # found the first True
                    true_found = True
        # found zero or one True value
        return true_found

    def upgr(
        xr_ob: Union[xr.DataArray, xr.Dataset], dstr: str, dimtup: Tuple[str]
    ) -> Union[xr.DataArray, xr.Dataset]:

        ext_pos = {"X": "lon", "Y": "lat", "L": "Z", "T": "time"}

        def check_and_rep(
            xr_ob1: Union[xr.DataArray, xr.Dataset], var: str, dstr1: str
        ) -> Union[xr.DataArray, xr.Dataset]:
            d_l = [var]
            for i in range(0, 5):
                d_l.append(var + "_0" + str(i))
            d_l.append(var + "u")
            d_l.append(var + "v")
            d_l.append(ext_pos[var])
            d_l.append(var.lower())
            d_l.append(ext_pos[var].lower())

            # check that the dimension is within the possiblities.
            assert dstr1 in d_l
            # check that there is only one dimension down that axis.
            assert only1([dimstr2 in d_l for dimstr2 in dimtup])
            # (otherwise things might get messed up)
            # rename that dimension to the default.
            if var != "L":
                return xr_ob1.rename({dstr1: var})
            else:
                return xr_ob1.rename({dstr1: "Z"})

        if "X" in dstr or "x" in dstr or "lon" in dstr:
            xr_ob = check_and_rep(xr_ob, "X", dstr)
        elif "Y" in dstr or "y" in dstr or "lat" in dstr:
            xr_ob = check_and_rep(xr_ob, "Y", dstr)
        elif "L" in dstr or "z" in dstr or "Z" in dstr:
            xr_ob = check_and_rep(xr_ob, "L", dstr)
        elif "T" in dstr or "t" in dstr or "time" in dstr:
            xr_ob = check_and_rep(xr_ob, "T", dstr)
        else:
            # what is this dimension, break this.
            print("warning, not changing: ", dstr)
            # assert False

        return xr_ob

    dims = xr_obj.dims

    assert isinstance(dims, (tuple, xr.core.utils.Frozen))

    for dim in dims:
        assert isinstance(dim, str)
        xr_obj = upgr(xr_obj, dim, dims)

    return xr_obj


def sel(
    xr_obj: Union[xr.Dataset, xr.DataArray], reg="pac"
) -> Union[xr.Dataset, xr.DataArray]:
    """
    Select a region of the dataset or datarray.

    Assumes

    reg options: "pac", "nino1+2", "nino3", "nino3.4", "nino3", "glob"
    https://www.ncdc.noaa.gov/teleconnections/enso/indicators/sst/

    From Figure 1:
    Distribution of 60-year trends in the NINO3.4 SST index
    (SST averaged over 5° S−5° N and 170° W−120° W) for end dates from 2008–2017.

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

    # X in degrees east, Y in degrees north.

    sel_d = {
        # [west, east] [south, north] boundaries
        "pac": {"X": (100, 290), "Y": (-30, 30)},
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


def open_dataset(
    path: Union[str, pathlib.Path], use_can_coords: bool = False
) -> xr.Dataset:
    """
    Open a dataset and formats it.

    Will only work if there is only one set of each coordinate at the moment.

    TODO: make it able to open and decode each sort of object
        if there are multiple time axes.

    Args:
        path (Union[str, pathlib.Path]): the path to the
            netcdf dataset file.
        can_coords (bool): whether or not to try and make
            the coordinate into the canonical names.

    Returns:
        xr.Dataset: The formatted dataset. Will have time variables decoded.
    """
    if not use_can_coords:
        return fix_calendar(xr.open_dataset(str(path), decode_times=False))
    else:
        can_coords(fix_calendar(xr.open_dataset(str(path), decode_times=False)))


def open_dataarray(path: Union[str, pathlib.Path]) -> xr.DataArray:
    """
    Open an xarray dataarray and format it.

    Will automatically try to make the dataset Coordinates
     into the canonical coordinate names (using can_coords).
    Will also decode the time axis.

    Args:
        path (Union[str, pathlib.Path]): the path to the netcdf datarray file.

    Returns:
        xr.DataArray: The formatted datarray.
    """
    return fix_calendar(can_coords(xr.open_dataarray(str(path), decode_times=False)))


def cut_and_taper(
    da: xr.DataArray,
    y_var: str = "Y",
    x_var: str = "X",
) -> xr.DataArray:
    """
    Cut and taper a field by latitude.

    Since the atmosphere model dynamics are only applicable
    in the tropics, the computed wind stress anomaly is only
    applied to the ocean model between 20° S and 20° N, and
    is linearly tapered to zero at 25° S and 25° N.

    Currently only copes if the array is two dimensional.

    Args:
        da (xr.DataArray): The datarray.
        y_var (str, optional): The name of the Y coordinate. Defaults to "Y".
        x_var (str, optional): The name of the X coordinate. Defaults to "X".

    Returns:
        xr.DataArray: The datarray with the function applied.

    Example:
        Should achieve::

            if da.Y > 25 or da.Y < -25:
                da = 0.0
            elif 20 <= da.Y <= 25:
                da = da - (0.2* (da.Y- 20))) * da
            else -20 >= da.Y >= -25:
                da = da - (0.2* (-da.Y - 20)) * da

        Usage::

            from src.xr_utils import open_dataset, cut_and_taper
            from src.constants import OCEAN_DATA_PATH
            da_new: xr.DataArray = open_dataset(OCEAN_DATA_PATH / "qflx.nc").qflx
            cut_and_taper(da_new.isel(Z=0, T=0, variable=0))

    """
    # make sure that they are in the correct order.
    da = da.transpose(y_var, x_var)

    @np.vectorize
    def test_vec(x: float, y: float):
        if y > 25 or y < -25:
            x = 0.0
        elif 20 <= y <= 25:
            x = x - (0.2 * (y - 20)) * x
        elif -20 >= y >= -25:
            x = x - (0.2 * (-y - 20)) * x
        return x

    for x in range(len(da.coords[x_var].values)):
        da[:, x] = test_vec(da[:, x], da.coords[y_var])

    return da


def spatial_mean(da: xr.DataArray) -> xr.DataArray:
    # pylint: disable=anomalous-backslash-in-string
    """
    Average a datarray over "X" and "Y" coordinates.

    Spatially weighted.

    Originally from:
    https://ncar.github.io/PySpark4Climate/tutorials/Oceanic-Ni%C3%B1o-Index/
    (although their version is wrong as it assumes numpy input is degrees)

    https://numpy.org/doc/stable/reference/generated/numpy.cos.html
    https://numpy.org/doc/stable/reference/generated/numpy.radians.html


    .. math::

        \\[
        \\bar{T}_{\\text {month }}=\\frac{\\sum_{j=1}^{n L a t}
        \\cos \\left(\\text { lat }_{j}\\right)
        \\bar{T}_{\\text {lat }, j}}{\\sum_{j=1}^{\\text {nLat }}
        \\cos \\left(\\text { lat }_{j}\\right)}
        \\]

    Args:
        da (xr.DataArray): da to average.

    Returns:
        xr.DataArray: avarage of da.
    """
    # Find mean temperature for each latitude
    mean_sst_lat = da.mean(dim="X")

    # Find Weighted mean of those values
    # https://numpy.org/doc/stable/reference/generated/numpy.cos.html
    # https://numpy.org/doc/stable/reference/generated/numpy.radians.html
    num = (np.cos(np.radians(da.Y)) * mean_sst_lat).sum(dim="Y")
    denom = np.sum(np.cos(np.radians(da.Y)))

    # Find mean global temperature
    mean_temp = num / denom

    return mean_temp
