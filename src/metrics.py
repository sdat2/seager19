"""Different metrics to calculate."""
from typing import Tuple
import matplotlib.pyplot as plt
import matplotlib
import xarray as xr
from src.constants import (
    NOAA_DATA_PATH,
    NINO3_4_TEST_PATH,
    DATA_PATH,
    SEL_DICT,
    FIGURE_PATH,
    CD_LOGS,
)
from src.xr_utils import (
    open_dataarray,
    sel,
    can_coords,
    spatial_mean,
    get_clim,
    open_dataset,
    get_trend,
    min_clim,
)
from src.plot_utils import (
    add_units,
    ps_defaults,
    get_dim,
    label_subplots,
    cmap,
)
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

    Args:
        sst (xr.DataArray): Sea surface temperature datarray in standard format.
        reg (str, optional): The region to select form src.xr_utils.sel.
            Defaults to "nino3.4".
        roll_period (int, optional): The rolling period defined with respect to
            the time axes. Defaults to 3.

    Returns:
        Tuple[xr.DataArray, xr.DataArray]: metric timeseries, climateology
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
    noaa_da = add_units(can_coords(xr.open_dataarray(NOAA_DATA_PATH)))
    noaa_da.attrs["units"] = r"$^{\circ}$C"
    return noaa_da


def calculate_nino3_4_from_noaa() -> Tuple[xr.DataArray, xr.DataArray]:
    """
    Calculate the default nino3.4 region from noaa data.

    Returns:
        Tuple[xr.DataArray, xr.DataArray]: metric timeseries, climateology
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


def plot_nino(ax: matplotlib.axes.Axes) -> None:
    """
    Plot nino.
    """

    def get_points(reg_dict: dict) -> Tuple[list]:
        """Get the rectangle."""
        x, y = [], []
        x.append(reg_dict["X"][0])
        y.append(reg_dict["Y"][0])
        x.append(reg_dict["X"][0])
        y.append(reg_dict["Y"][1])
        x.append(reg_dict["X"][1])
        y.append(reg_dict["Y"][1])
        x.append(reg_dict["X"][1])
        y.append(reg_dict["Y"][0])
        x.append(reg_dict["X"][0])
        y.append(reg_dict["Y"][0])
        return x, y

    for reg in reversed(sorted(SEL_DICT)):

        x, y = get_points(SEL_DICT[reg])
        # if False: # reg == "nino3.4":
        # metric
        # plt.fill(x, y, label=reg, alpha=0.5,
        # linewidth=1, color=SEL_DICT[reg]["color"])
        # else:
        ax.plot(x, y, label=reg, alpha=0.5, linewidth=2, color=SEL_DICT[reg]["color"])

    ax.set_title("")

    ax.set_xlim(95, 295)
    ax.set_ylim(-32, 32)

    ax.legend(
        bbox_to_anchor=(-0.02, 1.02, 1.15, 0.102),
        loc="lower left",
        mode="expand",
        ncol=5,
    )


def get_nino_trend(
    path_of_run2f: str,
    graph_path: str,
    nc_path: str,
) -> dict:
    """
    Get nino trend, mean, plot the graph.

    Args:
        path_of_run2f (str): path to the main output netcdf.
        graph_path (str): path to output plot.
        nc_path (str): path to output netcdf

    Returns:
        dict: nino dict.
    """
    plt.clf()
    ps_defaults(dpi=150)
    if "NOAA" not in path_of_run2f:
        sst_output = can_coords(open_dataset(path_of_run2f).SST_SST)
        sst_output = sst_output.where(sst_output != 0.0)
    else:
        sst_output = add_units(can_coords(open_dataarray(path_of_run2f)))

    _, axs = plt.subplots(3, 1, figsize=get_dim(ratio=1.2))

    metric_l = list()
    clim_l = list()
    reg_l = list()
    nino_dict = dict()

    add_units(get_trend(sst_output, min_clim_f=True, output="rise")).plot(
        ax=axs[0],
        cmap=cmap("delta"),
        cbar_kwargs={"label": r"$\Delta T_s$ [$\Delta$K]"},
    )
    plot_nino(axs[0])

    plt.xlim(95, 295)
    plt.ylim(-32, 32)

    for reg in reversed(sorted(SEL_DICT)):
        print(reg)

        metric, clim = nino_calculate(sst_output, reg=reg)
        metric.attrs["long_name"] = "3 month rolling average SST anomaly"
        clim.attrs["long_name"] = clim.attrs["long_name"][0:29]
        # metric.plot(label=metric.attrs["reg"])

        metric.plot(
            ax=axs[1],
            label=metric.attrs["reg"],
            color=SEL_DICT[reg]["color"],
            linewidth=0.5,
        )
        axs[1].set_title("")

        clim.plot(
            ax=axs[2],
            label=metric.attrs["reg"],
            color=SEL_DICT[reg]["color"],
            linewidth=0.5,
        )

        axs[2].set_title("")

        nino_dict["trend_" + reg] = get_trend(metric)
        nino_dict["mean_" + reg] = metric.attrs["mean_state"]
        metric_l.append(metric)
        clim_l.append(clim)
        reg_l.append(reg)

    axs[1].set_xlim(metric.coords["T"].values[0], metric.coords["T"].values[-1])

    plt.title("")
    axs[1].set_xlabel("")
    axs[2].set_xlabel("Month")
    label_subplots(axs, x_pos=-0.05, y_pos=1.24)
    axs[2].set_xlim(1, 12)
    axs[2].set_ylim(20, 30)
    plt.tight_layout()
    plt.savefig(graph_path)
    plt.clf()

    print(clim)

    if "NOAA" not in path_of_run2f:

        anom = (
            xr.concat(metric_l, "reg")
            .assign_coords(reg=reg_l)
            .isel(Z=0)
            .drop("Z")
            .to_dataset(name="anomaly")
        )
        clim = (
            xr.concat(clim_l, "reg")
            .assign_coords(reg=reg_l)
            .isel(Z=0)
            .drop("Z")
            .to_dataset(name="clim")
        )
        xr.merge([anom, clim]).to_netcdf(nc_path)

        return nino_dict


if __name__ == "__main__":
    # python src/metrics.py

    cfg = load_config()

    setup = ModelSetup(
        str("logs/it_1a"),
        cfg,
        make_move=False,
    )

    get_nino_trend(
        setup.om_run2f_nc(),
        str(FIGURE_PATH / "nino_it1_2.25_trend.png"),
        str(DATA_PATH / "cd.nc"),
    )
    get_nino_trend(
        setup.om_run2f_nc(),
        str(FIGURE_PATH / "nino_it1_2.25_trend.pdf"),
        str(DATA_PATH / "cd.nc"),
    )

    setup = ModelSetup(
        str(CD_LOGS / "cd_2.25" / "wandb" / "latest-run" / "files"),
        cfg,
        make_move=False,
    )

    get_nino_trend(
        setup.om_run2f_nc(),
        str(FIGURE_PATH / "nino_2.25_trend.png"),
        str(DATA_PATH / "cd.nc"),
    )
    get_nino_trend(
        setup.om_run2f_nc(),
        str(FIGURE_PATH / "nino_2.25_trend.pdf"),
        str(DATA_PATH / "cd.nc"),
    )

    # print("main")
    _, _ = calculate_nino3_4_from_noaa()
    get_nino_trend(
        str(NOAA_DATA_PATH),
        str(FIGURE_PATH / "nino_noaa_trend.png"),
        str(DATA_PATH / "noaa_trend.nc"),
    )
    get_nino_trend(
        str(NOAA_DATA_PATH),
        str(FIGURE_PATH / "nino_noaa_trend.pdf"),
        str(DATA_PATH / "noaa_trend.nc"),
    )
