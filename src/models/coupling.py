"""Coupling between ocean and atmos models is through the sea surface stress.

Example:
    Import statement usage::
        from src.models.coupling import f_stress

"""
from typing import Tuple
import numpy as np
import pandas as pd
import xarray as xr
from typeguard import typechecked


rho_a: float = 1.225  # kg m-3
cd: float = 2.25e-3  # Pa m-1 s


@typechecked
def f_stress(
    rho_air: float,
    c_d: float,
    wind_speed_mean: float,
    u_wind: xr.DataArray,
    v_wind: xr.DataArray,
) -> Tuple[xr.DataArray, xr.DataArray]:
    """Wind stress flux.

    Args:
        rho_air (float): density of sst
        c_d (float): wind stress.
        wind_speed_mean: the climatological annual
            mean wind speed, which is taken from ECMWF
            reanalysis for our standard model and from the CMIP5
            multimodel mean when examining causes of bias in the
            CMIP5 model.
        u_wind (xr.DataArray): The zonal wind speed.
        v_wind (xr.DataArray): The meridional wind speed.

    Returns:
        Tuple[xr.DataArray, xr.DataArray]: [zonal wind stress, meridional wind stress]

    """
    stress_coeff = rho_air * c_d * wind_speed_mean
    return stress_coeff * u_wind, stress_coeff * v_wind


def stress_test() -> None:
    """Stress test."""
    time = pd.date_range("2014-09-06", periods=3)
    reference_time = pd.Timestamp("2014-09-05")
    u_vel = 15 + 8 * np.random.randn(2, 2, 3)
    v_vel = 10 * np.random.rand(2, 2, 3)
    lon = [[-99.83, -99.32], [-99.79, -99.23]]
    lat = [[42.25, 42.21], [42.63, 42.59]]

    u_da = xr.DataArray(
        data=u_vel,
        dims=["x", "y", "time"],
        coords=dict(
            lon=(["x", "y"], lon),
            lat=(["x", "y"], lat),
            time=time,
            reference_time=reference_time,
        ),
        attrs=dict(),
    )

    v_da = xr.DataArray(
        data=v_vel,
        dims=["x", "y", "time"],
        coords=dict(
            lon=(["x", "y"], lon),
            lat=(["x", "y"], lat),
            time=time,
            reference_time=reference_time,
        ),
        attrs=dict(),
    )

    tau_u, tau_v = f_stress(rho_a, cd, 0.5, u_da, v_da)

    assert isinstance(tau_u, xr.DataArray)
    assert isinstance(tau_v, xr.DataArray)


if __name__ == "__main__":
    # python3 src/models/coupling.py
    stress_test()
