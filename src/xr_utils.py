"""Utilities around opening and processing netcdfs from this project."""
import numpy as np
import pathlib
from typing import Union, Tuple, Optional, Literal
import xarray as xr
from uncertainties import ufloat
from src.plot_utils import add_units
from src.constants import SEL_DICT, MASK


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

    # decode times into normal datetime
    ds = xr.decode_cf(ds)

    # transform back into original type
    if isinstance(xr_in, xr.DataArray):
        xr_out = ds.to_array()
    else:
        xr_out = ds

    return xr_out


def _mon_increase(
    xr_obj: Union[xr.Dataset, xr.DataArray]
) -> Union[xr.Dataset, xr.DataArray]:
    """Make sure that an xarray axes has monotonically increasing values"""

    def positive_monotonic(A):
        return all(A[i] <= A[i + 1] for i in range(len(A) - 1))

    def negative_monotonic(A):
        return all(A[i] >= A[i + 1] for i in range(len(A) - 1))

    for var in ["X", "Y"]:
        if negative_monotonic(xr_obj.coords[var].values):
            xr_obj = xr_obj.reindex(**{var: xr_obj.coords[var][::-1]})
        elif not positive_monotonic(xr_obj.coords[var].values):
            assert False

    return xr_obj


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

    xr_obj = _mon_increase(xr_obj)

    return xr_obj


def sel(
    xr_obj: Union[xr.Dataset, xr.DataArray], reg: str = "pac"
) -> Union[xr.Dataset, xr.DataArray]:
    """
    Select a region of the dataset or datarray.

    Assumes

    reg options: "pac", "nino1+2", "nino3", "nino3.4", "nino3"
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

    assert reg in SEL_DICT

    sel_c = SEL_DICT[reg]

    return xr_obj.sel(
        X=slice(sel_c["X"][0], sel_c["X"][1]), Y=slice(sel_c["Y"][0], sel_c["Y"][1])
    )


def rem_var(da1: xr.DataArray) -> xr.DataArray:
    """remove the 'variable' from a dataarray"""
    if "variable" in da1.dims:
        da1 = da1.isel(variable=0).drop("variable")
    return da1


def clip(da: xr.DataArray, pac: bool = True, mask_land: bool = True) -> xr.DataArray:
    """
    Clip a datarray to the pacific using sel, and mask the land.


    Args:
        da (xr.DataArray): The datarray to pass in.
        pac (bool, optional): Whether to focus on pacific. Defaults to True.
        mask_land (bool, optional): Whether to nan out land. Defaults to True.

    Returns:
        xr.DataArray: da with those operations applied to it.
    """
    attrs = da.attrs
    mask = open_dataset(MASK).mask
    mask = rem_var(mask)
    da = fix_calendar(da.rename("unknown"))
    da = add_units(rem_var(da))
    if pac and not mask_land:
        da = sel(da)
    elif pac and mask_land:
        try:
            da = sel(da).where(sel(mask) != 0.0)
            # pylint: disable=bare-except
        except:
            print(da, type(da), mask)
            da = sel(da)
    elif not pac and mask_land:
        da = da.where(mask != 0.0)

    da.attrs = attrs
    return da


def open_dataset(
    path: Union[str, pathlib.Path], use_can_coords: bool = False
) -> xr.Dataset:
    """
    Open a dataset and formats it.

    Will only work if there is only one set of each coordinate at the moment.

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
        return can_coords(fix_calendar(xr.open_dataset(str(path), decode_times=False)))


def open_dataarray(path: Union[str, pathlib.Path]) -> xr.DataArray:
    """
    Open an xarray dataarray and format it.

    Will automatically try to make the dataset Coordinates
     into the canonical coordinate names (using can_coords).
    Will also decode the time axis.

    TODO: add option for opening of datarrayys that just ensures they open, rather than changing their atrributes.

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
    # TODO: More global checking of order would improve relaibility.
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

    The average should behave like:

    .. math::
        :nowrap:

        \\begin{equation}
            \\bar{T}_{\\text {lat }}=\\frac{1}{n \\text { Lon }}
            \\sum_{i=1}^{n \\text{Lon}} T_{\\text \\text{lon}, i}
        \\end{equation}

        \\begin{equation}
            \\bar{T}_{\\text {month }}=\\frac{\\sum_{j=1}^{n L a t}
            \\cos \\left(\\text { lat }_{j}\\right)
            \\bar{T}_{\\text {lat }, j}}{\\sum_{j=1}^{\\text
            {n \\text{Lat} }}
            \\cos \\left(\\text { lat }_{j}\\right)}
        \\end{equation}

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


def get_trend(
    da: xr.DataArray,
    min_clim_f: bool = False,
    output: Literal["slope", "rise"] = "rise",
    t_var: str = "T",
    make_hatch_mask: bool = False,
    keep_ds: bool = False,
    uncertainty: bool = False,
) -> Union[float, ufloat, xr.DataArray, Tuple[xr.DataArray, xr.DataArray]]:
    """
    Returns either the linear trend rise, or the linear trend slope,
    possibly with the array to hatch out where the trend is not significant.

    Uses `xr.polyfit` order 1 to do everything.

    Args:
        da (xr.DataArray): the timeseries.
        min_clim_f (bool, optional): whether to calculate and remove the climateology.
            Defaults to false.
        output (Literal[, optional): What to return. Defaults to "rise".
        t_var (str, optional): The time variable name. Defaults to "T".
            Could be changed to another variable that you want to fit along.
        make_hatch_mask (bool, optional): Whether or not to also return a DataArray
            of boolean values to indicate where is not significant.
            Defaults to False. Will only work if you're passing in an xarray object.
        uncertainty (bool, optional): Whether to return a ufloat object
            if doing linear regression on a single timeseries. Defaults to false.

    Returns:
        Union[float, ufloat, xr.DataArray, Tuple[xr.DataArray, xr.DataArray]]:
            The rise/slope over the time period, possibly with the hatch array if
            that opition is selected for a grid.
    """

    def length_time(da):
        if t_var == "T":
            return int((da.coords[t_var][-1] - da.coords[t_var][0]).values)
        else:
            return (da.coords[t_var][-1] - da.coords[t_var][0]).values

    if min_clim_f:
        da = min_clim(da)

    def get_float(inp: Union[np.ndarray, list, float]):
        try:
            if hasattr(inp, "__iter__"):
                inp = inp[0]
        # pylint: disable=bare-except
        except:
            print(type(inp))
        return float(inp)

    if "X" in da.dims or "Y" in da.dims or "member" in da.dims or keep_ds:

        fit_da = da.polyfit(t_var, 1, cov=make_hatch_mask)

        if make_hatch_mask:
            frac_error = np.abs(
                np.sqrt(fit_da.polyfit_covariance.isel(cov_i=0, cov_j=0))
                / fit_da.polyfit_coefficients.isel(degree=0)
            )
            hatch_mask = frac_error >= 1.0

        slope = fit_da.polyfit_coefficients.sel(degree=1).drop("degree")
    else:
        if uncertainty:
            # print("uncertainty running")
            fit_da = da.polyfit(t_var, 1, cov=True)
            error = np.sqrt(fit_da.polyfit_covariance.isel(cov_i=0, cov_j=0)).values
            error = get_float(error)
            slope = fit_da.polyfit_coefficients.values
            slope = get_float(slope)
            slope = ufloat(slope, error)
        else:
            slope = da.polyfit(t_var, 1).polyfit_coefficients.values
            slope = get_float(slope)

    if output == "rise":

        run = length_time(da)
        rise = slope * run

        if isinstance(rise, xr.DataArray):
            for pr_name in ["units", "long_name"]:
                if pr_name in da.attrs:
                    rise.attrs[pr_name] = da.attrs[pr_name]
            rise = add_units(rise).rename("rise")

        # print("run", run, "slope", slope, "rise = slope * run", rise)

        if make_hatch_mask and not isinstance(rise, float, ufloat):
            return rise, hatch_mask
        else:
            return rise

    elif output == "slope":
        if make_hatch_mask and not isinstance(slope, float):
            return slope, hatch_mask
        else:
            return slope


def get_clim(xr_da: xr.DataArray) -> xr.DataArray:
    """
    Get the climateology of an xr.DataArray.

    Args:
        xr_da (xr.DataArray): The input datarray.
            Assumes that the time coordinate is canonical "T".

    Returns:
        xr.DataArray: The climatology for the time period.
    """
    time_coord = xr_da.coords["T"].values
    time_str = time_coord[0].__str__()[0:4] + " to " + time_coord[-1].__str__()[0:4]

    if "long_name" in xr_da.attrs:
        init_long_name = xr_da.attrs["long_name"]
    else:
        init_long_name = ""
    if "units" in xr_da.attrs:
        init_units = xr_da.attrs["units"]
    else:
        init_units = ""

    climatology = xr_da.groupby("T.month").mean("T")
    climatology.attrs["units"] = init_units
    climatology.attrs["long_name"] = (
        "Climateology for " + time_str + " " + init_long_name.lower()
    )
    return climatology


def min_clim(xr_da: xr.DataArray, clim: Optional[xr.DataArray] = None) -> xr.DataArray:
    """
    Take away the climatology from an xr.DataArray.

    Args:
        xr_da (xr.DataArray): The xarray input. Canonical coords.
        clim (Optional[xr.DataArray], optional): The climateology.
            Defaults to None, which will remake climatology.

    Returns:
        xr.DataArray: The anomaly.
    """
    if clim is None:
        clim = get_clim(xr_da)

    anom = xr_da.groupby("T.month") - clim

    for pr_name in ["units", "long_name"]:
        if pr_name in xr_da.attrs:
            anom.attrs[pr_name] = xr_da.attrs[pr_name]

    return anom
