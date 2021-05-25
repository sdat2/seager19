"""Coupling between ocean and atmos models is through the sea surface stress.

Example:
    Import statement usage::
        from src.models.coupling import f_stress

"""
from typing import Tuple
import xarray as xr
from typeguard import typechecked


class Couple:
    """
    Coupled model part.
    """

    def __init__(self):
        self.rho_air: float = 1.225  # kg m-3, density of sst
        self.c_d: float = 2.25e-3  # Pa m-1 s,  wind stress.

    @typechecked
    def f_stress(
        self,
        wind_speed_mean: float,
        u_wind: xr.DataArray,
        v_wind: xr.DataArray,
    ) -> Tuple[xr.DataArray, xr.DataArray]:
        """Wind stress flux.

        Args:
            wind_speed_mean (float): the climatological annual
                mean wind speed, which is taken from ECMWF
                reanalysis for our standard model and from the CMIP5
                multimodel mean when examining causes of bias in the
                CMIP5 model.
            u_wind (xr.DataArray): The zonal wind speed.
            v_wind (xr.DataArray): The meridional wind speed.

        Returns:
            Tuple[xr.DataArray, xr.DataArray]: [zonal wind stress,
                meridional wind stress]

        """
        stress_coeff = self.rho_air * self.c_d * wind_speed_mean
        return stress_coeff * u_wind, stress_coeff * v_wind
