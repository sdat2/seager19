"""Store plotting scripts for metrics.py"""
from typing import Tuple, Optional
import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
from src.constants import (
    SEL_DICT,
    NOAA_DATA_PATH,
    DATA_PATH,
    SEL_DICT,
    FIGURE_PATH,
    CD_LOGS,
)
from src.plot_utils import add_units, plot_defaults, cmap, label_subplots, get_dim
from src.xr_utils import (
    get_trend,
    open_dataarray,
    can_coords,
    open_dataset,
)
from src.metrics import nino_calculate, calculate_nino3_4_from_noaa
from src.configs.load_config import load_config
from src.models.model_setup import ModelSetup
from src.utils import timeit


def _get_points(reg_dict: dict) -> Tuple[list]:
    """Get the rectangle from the reg_dict."""
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


@timeit
def plot_nino_box(
    ax: matplotlib.axes.Axes, reg: str, color: Optional[str] = None
) -> None:
    """
    Plot a nino box.

    Args:
        ax (matplotlib.axes.Axes): axis to plot to.
        reg (str): region to plot
        color (Optional[str], optional): Override color option. Defaults to None.

    Example:
        Add Nino3.4 box to plot::

            import matplotlib.pyplot as plt
            from src.visualisation.nino import plot_nino_box

            ax = plt.gca()
            plot_nino_box(
                ax, reg="nino3.4", color=None
            )

    """
    x, y = _get_points(SEL_DICT[reg])
    if color is None:
        color = SEL_DICT[reg]["color"]
    ax.plot(x, y, label=reg, alpha=0.5, linewidth=2, color=color)


def plot_nino(ax: matplotlib.axes.Axes, legend: bool = False) -> None:
    """
    Plot all the nino boxes.

    Args:
        ax (matplotlib.axes.Axes): axis to plot to.
        legend (bool, optional): Whether to include legend. Defaults to False.
    """

    for reg in reversed(sorted(SEL_DICT)):
        plot_nino_box(ax, reg)
    ax.set_title("")

    ax.set_xlim(95, 295)
    ax.set_ylim(-32, 32)

    if legend:
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
    show_plots: bool = False,
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
    plot_defaults(dpi=150)
    # TODO: this is very hacky :'(  Fix this!
    if "NOAA" not in path_of_run2f and "ts" not in path_of_run2f:
        sst_output = can_coords(open_dataset(path_of_run2f).SST_SST)
        sst_output = sst_output.where(sst_output != 0.0)
    elif "ts" in path_of_run2f:
        sst_output = add_units(can_coords(open_dataarray(path_of_run2f))) - 273.15
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

        rise = get_trend(metric, uncertainty=True)

        nino_dict["trend_" + reg] = rise.n
        nino_dict["trend_" + reg + "_unc"] = rise.s
        nino_dict["mean_" + reg] = metric.attrs["mean_state"]

        label = str(
            metric.attrs["reg"]
            # + "\n"
            + r" $\Delta T_s = $"
            # + "\n"
            "${:.2L}$".format(rise)
            # + tex_uf(rise)
            + r"$^{\circ}$C"
        )

        metric.plot(
            ax=axs[1],
            label=label,
            color=SEL_DICT[reg]["color"],
            linewidth=0.5,
        )
        axs[1].set_title("")

        clim.plot(
            ax=axs[2],
            label=label,
            color=SEL_DICT[reg]["color"],
            linewidth=1.5,
        )

        axs[2].set_title("")

        metric_l.append(metric)
        clim_l.append(clim)
        reg_l.append(reg)

    axs[1].set_xlim(metric.coords["T"].values[0], metric.coords["T"].values[-1])

    axs[2].legend(
        bbox_to_anchor=(-0.02, 1.02, 1.0, 0.102),
        loc="lower left",
        mode="expand",
        ncol=2,
    )

    plt.title("")
    axs[1].set_xlabel("")
    axs[2].set_xlabel("Month")
    label_subplots(axs, x_pos=-0.05, y_pos=1.24)
    axs[2].set_xlim(1, 12)
    axs[2].set_ylim(20, 30)
    plt.tight_layout()
    plt.savefig(graph_path)

    if show_plots:
        plt.show()

    else:
        plt.clf()

    print(clim)

    def rem_z(da: xr.DataArray) -> xr.DataArray:
        if "Z" in da.dims:
            return da.isel(Z=0).drop("Z")
        else:
            return da

    if "NOAA" not in path_of_run2f:

        anom = rem_z(xr.concat(metric_l, "reg").assign_coords(reg=reg_l)).to_dataset(
            name="anomaly"
        )
        clim = rem_z(xr.concat(clim_l, "reg").assign_coords(reg=reg_l)).to_dataset(
            name="clim"
        )
        xr.merge([anom, clim]).to_netcdf(nc_path)

        return nino_dict


def make_plots() -> None:
    """
    Make sample plots for nino3.4 to check its still working.
    """
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


def multi_panel_nino(sst_output, graph_path, show_plots=False):
    """
    Multi panel nino plot.

    Args:
        sst_output (xr.DataArray): sst_output netcdf.
        graph_path (str): Graph path.
        show_plots (bool, optional): Show plots. Defaults to False.
    """
    _, axs = plt.subplots(3, 1, figsize=get_dim(ratio=1.2))

    metric_l = []
    clim_l = []
    reg_l = []
    nino_dict = []

    add_units(get_trend(sst_output, min_clim_f=True, output="rise")).plot(
        ax=axs[0],
        cmap=cmap("delta"),
        cbar_kwargs={"label": r"$\Delta T_s$ [$\Delta$K]"},
    )
    plot_nino(axs[0])

    plt.xlim(95, 295)
    plt.ylim(-32, 32)
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

        rise = get_trend(metric, uncertainty=True)

        nino_dict["trend_" + reg] = rise.n
        nino_dict["trend_" + reg + "_unc"] = rise.s
        nino_dict["mean_" + reg] = metric.attrs["mean_state"]

        label = str(
            metric.attrs["reg"]
            # + "\n"
            + r" $\Delta T_s = $"
            # + "\n"
            "${:.2L}$".format(rise)
            # + tex_uf(rise)
            + r"$^{\circ}$C"
        )

        metric.plot(
            ax=axs[1],
            label=label,
            color=SEL_DICT[reg]["color"],
            linewidth=0.5,
        )
        axs[1].set_title("")

        clim.plot(
            ax=axs[2],
            label=label,
            color=SEL_DICT[reg]["color"],
            linewidth=1.5,
        )

        axs[2].set_title("")

        metric_l.append(metric)
        clim_l.append(clim)
        reg_l.append(reg)

    axs[1].set_xlim(metric.coords["T"].values[0], metric.coords["T"].values[-1])

    axs[2].legend(
        bbox_to_anchor=(-0.02, 1.02, 1.0, 0.102),
        loc="lower left",
        mode="expand",
        ncol=2,
    )

    plt.title("")
    axs[1].set_xlabel("")
    axs[2].set_xlabel("Month")
    label_subplots(axs, x_pos=-0.05, y_pos=1.24)
    axs[2].set_xlim(1, 12)
    axs[2].set_ylim(20, 30)
    plt.tight_layout()
    plt.savefig(graph_path)

    if show_plots:
        plt.show()

    else:
        plt.clf()
