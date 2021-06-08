"""Different metrics to calculate."""
from typing import Tuple
import xarray as xr
from src.constants import NOAA_DATA_PATH
from src.xr_utils import sel, can_coords
from src.plot_utils import add_units


def nino_calculate(
    sst: xr.DataArray, reg: str = "nino3.4", roll_period: int = 3
) -> Tuple[xr.DataArray, xr.DataArray]:
    """
    Calculate the nino metric for a given region.

    Can work on nino.

    Args:
        sst (xr.DataArray): [description]
        reg (str, optional): [description]. Defaults to "nino3.4".
        roll_period (int, optional): The rolling period defined with respect to
            the time axes. Defaults to 3.

    Returns:
        Tuple[xr.DataArray, xr.DataArray]: [description]
    """
    can_coords(sst)

    sst_reg = sel(sst, reg=reg)

    mean_timeseries = sst_reg.mean(dim=["X", "Y"])
    mean_timeseries.attrs["long_name"] = (
        "Sea surface temperature averaged over " + reg + " region"
    )
    mean_timeseries.attrs["units"] = r"$^{\circ}$C"
    mean_timeseries.coords["T"].attrs["long_name"] = "Month"

    mean_state = mean_timeseries.mean(dim=["T"])
    mean_state.attrs["long_name"] = (
        "Average sea surface temperature over " + reg + " region"
    )
    mean_state.attrs["units"] = r"$^{\circ}$C"
    anomaly_timeseries = mean_timeseries - mean_state
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
    metric.attrs["rolling_average"] = str(roll_period) + " months"
    return metric, mean_state


def calculate_nino3_4_from_noaa() -> xr.DataArray:
    """
    Calculate the default nino3_4.

    Returns:
        xr.DataArray: [description]
    """

    def load_noaa_data() -> xr.DataArray:
        noaa_da = add_units(can_coords(xr.open_dataarray(NOAA_DATA_PATH)))
        noaa_da.attrs["units"] = r"$^{\circ}$C"
        return noaa_da

    return nino_calculate(load_noaa_data())
