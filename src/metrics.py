"""Different metrics to calculate

from src.visualisation.nino import plot_nino_trend
"""
import os
from typing import Tuple
import xarray as xr
from src.constants import (
    NOAA_DATA_PATH,
    NINO3_4_TEST_PATH,
    DATA_PATH,
    SEL_DICT,
    EPS_FRAC_LOGS,
)
from src.xr_utils import (
    sel,
    can_coords,
    spatial_mean,
    get_clim,
    min_clim,
    get_trend,
    open_dataset,
)
from src.utils import timeit
from src.plot_utils import add_units
from src.configs.load_config import load_config
from src.models.model_setup import ModelSetup


def nino_calculate(
    sst: xr.DataArray, reg: str = "nino3.4", roll_period: int = 3
) -> Tuple[xr.DataArray, xr.DataArray]:
    """
    Calculate the nino metric for a given region.

    https://rabernat.github.io/research_computing_2018/assignment-8-xarray-for-enso.html

    https://ncar.github.io/PySpark4Climate/tutorials/Oceanic-Ni%C3%B1o-Index/

    Can work on regions nino1+2, nino3, nino3.4, nino4 (or "pac").

    "pac" is a region defined by me mainly for plotting that
    includes most of the tropical pacific.

    Args:
        sst (xr.DataArray): Sea surface temperature datarray in standard format.
        reg (str, optional): The region to select for src.xr_utils.sel.
            Defaults to "nino3.4".
        roll_period (int, optional): The rolling period defined with respect to
            the time axes. Defaults to 3.

    Returns:
        Tuple[xr.DataArray, xr.DataArray]: metric timeseries, climatology
    """
    can_coords(sst)

    sst_reg = sel(sst, reg=reg)

    mean_timeseries = spatial_mean(sst_reg)
    mean_timeseries.attrs["long_name"] = (
        "Sea surface temperature averaged over " + reg + " region"
    )
    mean_timeseries.attrs["units"] = r"$^{\circ}$C"
    mean_timeseries.coords["T"].attrs["long_name"] = "Month"
    climatology = get_clim(mean_timeseries)
    time_coord = mean_timeseries.coords["T"].values
    time_str = time_coord[0].__str__()[0:4] + " to " + time_coord[-1].__str__()[0:4]
    mean_state = mean_timeseries.mean(dim=["T"])
    mean_state.attrs["long_name"] = (
        "Average sea surface temperature over " + reg + " region"
    )
    mean_state.attrs["units"] = r"$^{\circ}$C"
    anomaly_timeseries = min_clim(mean_timeseries, clim=climatology)
    anomaly_timeseries = mean_timeseries.groupby("T.month") - climatology
    anomaly_timeseries.attrs["long_name"] = (
        "Sea surface temperature anomaly averaged over " + reg + " region"
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
    metric.attrs["mean_state"] = float(mean_state.values)
    # metric.attrs["clim"] = climatology.values  # .tolist()
    # metric.attrs["clim_months"] = climatology.month.values  # .tolist()
    metric.attrs["clim_time_range"] = time_str
    return metric, climatology


def load_noaa_data() -> xr.DataArray:
    """
    Load the data from the noaa ERSSTv4.5 file.

    Returns:
        xr.DataArray: NOAA dataarray.
    """
    noaa_da = add_units(can_coords(xr.open_dataarray(NOAA_DATA_PATH)))
    noaa_da.attrs["units"] = r"$^{\circ}$C"
    return noaa_da


def calculate_nino3_4_from_noaa() -> Tuple[xr.DataArray, xr.DataArray]:
    """
    Calculate the default nino3.4 region from noaa data.

    Returns:
        Tuple[xr.DataArray, xr.DataArray]: metric timeseries, climatology
    """

    return nino_calculate(load_noaa_data())


def replace_nino3_4_from_noaa() -> None:
    """
    Calculate the default nino3.4 region from noaa data.
    """
    metric, clim = calculate_nino3_4_from_noaa()
    metric.to_netcdf(str(NINO3_4_TEST_PATH))
    print(NINO3_4_TEST_PATH)
    clim.to_netcdf(str(DATA_PATH / "nino3_4_noaa_clim.nc"))


@timeit
def get_other_trends(
    setup: ModelSetup,
) -> dict:
    """
    Get trends in nino regions for other variables other than sst.

    Args:
        setup (ModelSetup): the filespace object to find things using.

    Returns:
        dict: nino dict.
    """
    nino_dict = dict()

    for field in ["TDEEP_HMODEL", "SST_QNET", "SST_W1"]:

        output = can_coords(open_dataset(setup.om_run2f_nc())[field])
        output = output.where(output != 0.0)

        metric_l = list()
        clim_l = list()
        reg_l = list()

        for reg in reversed(sorted(SEL_DICT)):
            print(reg, field)

            metric, clim = nino_calculate(output, reg=reg)
            metric.attrs["long_name"] = "3 month rolling average SST anomaly"
            clim.attrs["long_name"] = clim.attrs["long_name"][0:29]
            # metric.plot(label=metric.attrs["reg"])

            rise = get_trend(metric, uncertainty=True)

            nino_dict["trend_" + field + "_" + reg] = rise.n
            nino_dict["trend_" + field + "_" + reg + "_unc"] = rise.s
            nino_dict["mean_" + field + "_" + reg] = metric.attrs["mean_state"]

            label = str(
                metric.attrs["reg"]
                + " "
                + field
                # + "\n"
                + r" $\Delta = $"
                # + "\n"
                "${:.2L}$".format(rise)
                # + tex_uf(rise)
            )
            print(label)

            metric_l.append(metric)
            clim_l.append(clim)
            reg_l.append(reg)

    for field in ["PRtrend", "utrend", "vtrend"]:

        tcam_output = can_coords(open_dataset(setup.tcam_output())[field])
        for reg in reversed(sorted(SEL_DICT)):
            nino_dict["trend_" + field + "_" + reg] = float(
                spatial_mean(sel(tcam_output, reg=reg)).values
            )

    return nino_dict


if __name__ == "__main__":
    # python src/metrics.py
    print(os.listdir(EPS_FRAC_LOGS))
    config = load_config()
    stp = ModelSetup(
        str(EPS_FRAC_LOGS / "k_days_10_eps_days_0.75_efrac_0.25_c_d_0.00225"),
        config,
        make_move=False,
    )
    print(get_other_trends(stp))
