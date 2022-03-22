"""Regrid the data


Scripts initially adapted from:

https://github.com/pangeo-data/xESMF/blob/master/xesmf/util.py
"""
from typing import Union, Tuple
import numpy as np
import xarray as xr
import xesmf as xe
import warnings
from src.xr_utils import can_coords


def _grid_1d(
    start_b: float, end_b: float, step: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    1D grid centers and bounds

    Args:
        start_b (float): start position. Bounds, not centers.
        end_b (float): end position. Bounds, not centers.
        step (float): step size, i.e. grid resolution

    Returns:
        Tuple[np.ndarray, np.ndarray]: centers with 1D numpy array
            bounds with 1D numpy array, with one more element than centers
    """

    bounds = np.arange(start_b, end_b + step, step)
    centers = (bounds[:-1] + bounds[1:]) / 2

    return centers, bounds


def grid_1d(lon0_b, lon1_b, d_lon, lat0_b, lat1_b, d_lat):
    """
    1D rectilinear grid centers and bounds.

    Adapted to enforce 0 to 360 longitude range.

    Args:
        lon0_b (float): Lower longitude bound.
        lon1_b (float): Longitude bounds.
        d_lon (float): Longitude step size, i.e. grid resolution.
        lat0_b (float): Latitude lower bound.
        lat1_b (float): Upper latitude bounds
        d_lat (float): Latitude step size, i.e. grid resolution.

    Returns:
        xr.Dataset: Dataset with coordinate values.
    """

    lon_1d, lon_b_1d = [x % 360 for x in _grid_1d(lon0_b, lon1_b, d_lon)]
    lat_1d, lat_b_1d = _grid_1d(lat0_b, lat1_b, d_lat)

    ds = xr.Dataset(
        coords={
            "lon": (["x"], lon_1d, {"standard_name": "longitude"}),
            "lat": (["y"], lat_1d, {"standard_name": "latitude"}),
            "lon_b": (["x_b"], lon_b_1d),
            "lat_b": (["y_b"], lat_b_1d),
        }
    )

    return ds


def _regridding_ds_1d(with_bounds=False):
    """
    Global 1D rectilinear grid centers, and bounds if required.

    The data inputs centre on xyz.

    Args:
        with_bounds (bool, optional): Whether to include bounds,
            which is required for the conservative regridding method.
            Defaults to False.

    Returns:
        xr.Dataset: Coordinate values
    """
    if with_bounds:
        return grid_1d(-0.5, 359.5, 1, -90, 90, 1)
    else:
        return grid_1d(-0.5, 359.5, 1, -90.5, 90.5, 1).drop(["lon_b", "lat_b"])


def regrid_2d(
    ds_input: Union[xr.Dataset, xr.DataArray],
    method: str = "bilinear",
    periodic: bool = True,
) -> Union[xr.Dataset, xr.DataArray]:
    """
    Regrid to a 1x1 grid, using by default bilinear interpolation and periodic boundaries.

    Args:
        ds_input (Union[xr.Dataset, xr.DataArray]): The xarray object to preprocess.
        method (str, optional): The method to choose. By default "bilinear".
        periodic (bool, optional): Whether or not to treat the boundaries as periodic.
            Defaults to True.

    Returns:
        Union[xr.Dataset, xr.DataArray]: The preprocessed xarray object.
    """
    regridder = xe.Regridder(
        ds_input,
        xe.util.grid_global(1, 1),
        method,
        periodic=periodic,
        ignore_degenerate=True,
        extrap_method="nearest_s2d",
    )
    ds_output = regridder(ds_input, keep_attrs=True)
    return ds_output


def regrid_2d_to_standard(
    da: Union[xr.Dataset, xr.DataArray]
) -> Union[xr.Dataset, xr.DataArray]:
    """Fix weird da structure returned by xESMf."""
    return can_coords(
        da.rename({"x": "X", "y": "Y"})
        .assign_coords(
            {
                "X": ("X", da.isel(y=0).lon.values % 360),
                "Y": ("Y", da.isel(x=0).lat.values),
            }
        )
        .drop_vars(["lon", "lat"])
        .sortby("X")
    )


def regrid_1d(
    ds_input: Union[xr.Dataset, xr.DataArray],
    method: str = "bilinear",
    periodic: bool = True,
) -> Union[xr.Dataset, xr.DataArray]:
    """
    Regrid to a 1x1 grid, using by default bilinear interpolation and periodic boundaries.

    Args:
        ds_input (Union[xr.Dataset, xr.DataArray]): The xarray object to preprocess.
        method (str, optional): The method to choose. By default "bilinear".
        periodic (bool, optional): Whether or not to treat the boundaries as periodic.
            Defaults to True.

    Returns:
        Union[xr.Dataset, xr.DataArray]: The preprocessed xarray object.
    """
    regridder = xe.Regridder(
        ds_input,
        _regridding_ds_1d(with_bounds=False),
        method,
        periodic=periodic,
        ignore_degenerate=True,
        extrap_method="nearest_s2d",
    )
    ds_output = regridder(ds_input, keep_attrs=True)
    return ds_output


def regrid_1d_to_standard(
    da: Union[xr.Dataset, xr.DataArray]
) -> Union[xr.Dataset, xr.DataArray]:
    """Change names to standard names."""
    return can_coords(
        da.rename({"x": "X", "y": "Y"})
        .assign_coords(
            {
                "X": ("X", da.isel(y=0).lon.values),
                "Y": ("Y", da.isel(x=0).lat.values),
            }
        )
        .drop_vars(["lon", "lat"])
    )


# can_coords(da.)


if __name__ == "__main__":
    # python src/data_loading/regrid.py
    print(regrid_1d_to_standard(_regridding_ds_1d(with_bounds=False)))
    print(regrid_2d_to_standard(xe.util.grid_global(1, 1).drop(["lon_b", "lat_b"])))
