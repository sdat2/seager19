"""Different metrics to calculate."""
from typing import Tuple
import numpy as np
import xarray as xr
from src.constants import NOAA_DATA_PATH
from src.xr_utils import sel, can_coords
from src.plot_utils import add_units


def nino_calculate(
    sst: xr.DataArray, reg: str = "nino3.4", roll_period: int = 3
) -> Tuple[xr.DataArray, xr.DataArray]:
    """
    Calculate the nino metric for a given region.

    https://rabernat.github.io/research_computing_2018/assignment-8-xarray-for-enso.html
    https://ncar.github.io/PySpark4Climate/tutorials/Oceanic-Ni%C3%B1o-Index/

    Can work on nino1+2, nino3, nino3.4, nino4.

    Args:
        sst (xr.DataArray): [description]
        reg (str, optional): [description]. Defaults to "nino3.4".
        roll_period (int, optional): The rolling period defined with respect to
            the time axes. Defaults to 3.

    Returns:
        Tuple[xr.DataArray, xr.DataArray]: metric, climateology
    """
    can_coords(sst)

    sst_reg = sel(sst, reg=reg)

    def mean_sst(sst_da: xr.DataArray) -> xr.DataArray:
        # Find mean temperature for each latitude
        mean_sst_lat = sst_da.mean(dim="X")

        # Find Weighted mean of those values
        # https://numpy.org/doc/stable/reference/generated/numpy.cos.html
        # https://numpy.org/doc/stable/reference/generated/numpy.radians.html
        num = (np.cos(np.radians(sst_da.Y)) * mean_sst_lat).sum(dim="Y")
        denom = np.sum(np.cos(np.radians(sst_da.Y)))

        # Find mean global temperature
        mean_temp = num / denom

        return mean_temp

    mean_timeseries = mean_sst(sst_reg)
    mean_timeseries.attrs["long_name"] = (
        "Sea surface temperature averaged over " + reg + " region"
    )
    mean_timeseries.attrs["units"] = r"$^{\circ}$C"
    mean_timeseries.coords["T"].attrs["long_name"] = "Month"
    time_coord = mean_timeseries.coords["T"].values
    time_str = time_coord[0].__str__()[0:4] + " to " + time_coord[-1].__str__()[0:4]
    climatology = mean_timeseries.groupby("T.month").mean("T")
    climatology.attrs["units"] = r"$^{\circ}$C"
    climatology.attrs["long_name"] = "Climateology for " + reg + " region " + time_str
    mean_state = mean_timeseries.mean(dim=["T"])
    mean_state.attrs["long_name"] = (
        "Average sea surface temperature over " + reg + " region"
    )
    mean_state.attrs["units"] = r"$^{\circ}$C"
    anomaly_timeseries = mean_timeseries.groupby("T.month") - climatology
    anomaly_timeseries.attrs["long_name"] = (
        "Sea surface temperature averaged over " + reg + " region"
    )
    anomaly_timeseries.attrs["units"] = r"$^{\circ}$C"
    metric = anomaly_timeseries.rolling(
        min_periods=1, T=roll_period, center=True
    ).mean()
    metric.attrs["long_name"] = (
        str(roll_period)
        + " month rolling average of sea surface "
        + "temperature averaged over "
        + reg
        + " region"
    )
    metric.attrs["units"] = r"$^{\circ}$C"
    metric.attrs["reg"] = reg
    metric.attrs["rolling_average"] = str(roll_period) + " months"
    metric.attrs["mean_state"] = mean_state.values
    metric.attrs["climatology"] = climatology.values
    metric.attrs["clim_months"] = climatology.month.values
    metric.attrs["clim_time_range"] = time_str
    return metric, climatology


def calculate_nino3_4_from_noaa() -> Tuple[xr.DataArray, xr.DataArray]:
    """
    Calculate the default nino3.4 region from noaa data.

    Returns:
        xr.DataArray: [description]
    """

    def load_noaa_data() -> xr.DataArray:
        noaa_da = add_units(can_coords(xr.open_dataarray(NOAA_DATA_PATH)))
        noaa_da.attrs["units"] = r"$^{\circ}$C"
        return noaa_da

    return nino_calculate(load_noaa_data())
