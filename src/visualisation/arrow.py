"""Arrow plots for mechanism."""
import os
from collections import OrderedDict
from typing import Optional
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from src.plot_utils import ps_defaults
from src.constants import FIGURE_PATH


def plot_arrow_plot(save_path: Optional[str] = None, show_plots: bool = False) -> None:
    """
    Plot the arrow plot to show that I have reproduced the paper.

    Args:
        save_path (Optional[str], optional): Where to save the plot to.
            Defaults to None. If None will not save.
        show_plots (bool, optional): Whether to show plots. Defaults to False.
    """
    ps_defaults(use_tex=False)

    color_d = {
        "EEEE": "blue",
        "EECE": "green",
        "EEEC": "orange",
        "EECC": "red",
    }

    def plot_error(x: float, y: float, yerr: float, mem: str) -> None:
        plt.fill_between(
            [x - 0.2, x + 0.2],
            [y + yerr, y + yerr],
            [y - yerr, y - yerr],
            color=color_d[mem],
            alpha=0.5,
        )
        plt.plot([x - 0.2, x + 0.2], [y, y], "black", linewidth=1)

    xlim = [0.5, 3.5]
    head_length = 0.02
    decrease_arrow = 0.01
    ax = plt.axes()
    ecmwf = 0.411
    # ax.arrow(0, 0, 0, 1, head_width=0.02, head_length=0.02, fc='k', ec='k')
    ax.arrow(
        1,
        ecmwf,
        0,
        0.054 - head_length - decrease_arrow,
        head_width=0.02,
        head_length=head_length,
        fc="k",
        ec="k",
    )
    plot_error(1, ecmwf + 0.054, 0.005, "EECE")
    ax.arrow(
        2,
        ecmwf,
        0,
        0.31 - head_length - decrease_arrow,
        head_width=0.02,
        head_length=head_length,
        fc="k",
        ec="k",
    )
    plot_error(2, ecmwf + 0.31, 0.03, "EEEC")
    ax.arrow(
        3,
        ecmwf,
        0,
        0.47 - head_length - decrease_arrow,
        head_width=0.02,
        head_length=head_length,
        fc="k",
        ec="k",
    )
    plot_error(3, ecmwf + 0.47, 0.04, "EECC")
    plt.plot(xlim, [ecmwf, ecmwf], color="blue", label="ECMWF/ORAS4 $= 0.411$ K ")
    plt.plot(
        xlim, [ecmwf + 0.478, ecmwf + 0.478], color="red", label="CMIP5 MMM $= 0.889$ K"
    )

    # plt.xticks([0, 1, 2, 3], ["ECMWF", "W", "RH", "RH+W"])
    plt.xticks(
        [1, 2, 3],
        [
            "W\n" + r"$+ 0.054 \pm 0.005$ K ",
            "RH\n " + r"$+ 0.31 \pm 0.03$ K",
            "RH+W\n " + r"$+ 0.47 \pm 0.04$ K",
        ],
    )

    plt.xlim(xlim)
    plt.ylabel("1958-2017, Trend in nino3.4 [K]")

    plt.legend(
        bbox_to_anchor=(0.0, 1.02, 1, 0.102),
        loc="lower left",
        mode="expand",
        ncol=2,
    )
    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path)
    if show_plots:
        plt.show()
    else:
        plt.clf()


def plot_arrow_plot_6(
    save_path: Optional[str] = None, show_plots: bool = False
) -> None:
    """
    Plot the arrow plot to show how it performs in cmip6.

    Args:
        save_path (Optional[str], optional): Where to save the plot to.
            Defaults to None. If None will not save.
        show_plots (bool, optional): Whether to show plots. Defaults to False.
    """
    ps_defaults(use_tex=False)

    color_d = {
        "EEEE": "blue",
        "EECE": "green",
        "EEEC": "orange",
        "EECC": "red",
    }

    def plot_error(x: float, y: float, yerr: float, mem: str) -> None:
        plt.fill_between(
            [x - 0.2, x + 0.2],
            [y + yerr, y + yerr],
            [y - yerr, y - yerr],
            color=color_d[mem],
            alpha=0.5,
        )
        plt.plot([x - 0.2, x + 0.2], [y, y], "black", linewidth=1)

    xlim = [0.5, 3.5]
    head_length = 0.02
    decrease_arrow = 0.01
    ax = plt.axes()
    ecmwf = 0.411
    # ax.arrow(0, 0, 0, 1, head_width=0.02, head_length=0.02, fc='k', ec='k')
    wind = 0.07
    wind_error = 0.01
    rh = 0.15
    rh_error = 0.02
    cmip6 = 0.772
    rh_and_wind = 0.29
    rh_and_wind_error = 0.04

    ax.arrow(
        1,
        ecmwf,
        0,
        wind - head_length - decrease_arrow,
        head_width=0.02,
        head_length=head_length,
        fc="k",
        ec="k",
    )
    plot_error(1, ecmwf + wind, wind_error, "EECE")
    ax.arrow(
        2,
        ecmwf,
        0,
        rh - head_length - decrease_arrow,
        head_width=0.02,
        head_length=head_length,
        fc="k",
        ec="k",
    )
    plot_error(2, ecmwf + rh, rh_error, "EEEC")
    ax.arrow(
        3,
        ecmwf,
        0,
        rh_and_wind - head_length - decrease_arrow,
        head_width=0.02,
        head_length=head_length,
        fc="k",
        ec="k",
    )
    plot_error(3, ecmwf + rh_and_wind, rh_and_wind_error, "EECC")
    plt.plot(xlim, [ecmwf, ecmwf], color="blue", label="ECMWF/ORAS4 $= 0.411$ K ")
    plt.plot(
        xlim,
        [cmip6, cmip6],
        color="red",
        label="CMIP6 MMM $= 0.772$ K",
    )

    # plt.xticks([0, 1, 2, 3], ["ECMWF", "W", "RH", "RH+W"])
    plt.xticks(
        [1, 2, 3],
        [
            "W\n"
            + r"$+ $"
            + str(wind)
            + r" $\pm$ "
            + r"$"
            + str(wind_error)
            + r"$"
            + " K ",
            "RH\n " + r"$+ $ $0.15$ $\pm$ $0.02$ K",
            "RH+W\n " + r"$+ $ $0.29$ $\pm$ $0.04$ K",
        ],
    )

    plt.xlim(xlim)
    plt.ylabel("1958-2017, Trend in nino3.4 [K]")

    plt.legend(
        bbox_to_anchor=(0.0, 1.02, 1, 0.102),
        loc="lower left",
        mode="expand",
        ncol=2,
    )
    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path)
    if show_plots:
        plt.show()
    else:
        plt.clf()


# Key Format: ${ts}${clt}${sfcwind}${rh}
RESULTS = OrderedDict(
    [
        ("EEEE", [0.436, 0.377, 0.461, 0.401]),  # All ECMWF inputs
        ("EECE", [0.492, 0.428, 0.522, 0.452]),  # cmip5 windspeed
        ("EEEC", [0.783, 0.646, 0.828, 0.686]),  # cmip5 relative humidity
        ("EECC", [0.948, 0.78, 1.002, 0.827]),  # cmip5 relative humidity and windspeed
        ("EE6E", [0.509, 0.436, 0.538, 0.462]),  # cmip6 windspeed
        ("EEE6", [0.601, 0.508, 0.638, 0.543]),  # cmip6 relative humidity
        ("EE66", [0.756, 0.627, 0.797, 0.664]),  # cmip6 relative humidity and windspeed
    ]
)
CMIP5_MMM = 0.889
CMIP6_MMM = 0.772
ECMWF = 0.411


def make_results_xr() -> xr.DataArray:
    """
    Make the results xarray from the dictionary above.

    Returns:
        xr.DataArray: The results xarray object.
    """
    key_list = []
    data_list = []
    for key in RESULTS:
        key_list.append(key)
        data_list.append(RESULTS[key])
    var_list = ["N", "A", "B", "C"]

    return xr.DataArray(
        data=np.array(data_list),
        dims=["mem", "var"],
        coords=dict(mem=(["mem"], key_list), var=(["var"], var_list)),
        attrs=dict(
            description="NINO3.4 trend 1958-2017",
            units=r"$\Delta$K",
            long_name="NINO3.4 change 1958-2017",
        ),
    ).rename("NINO3.4_TREND")


RESULTS_XR = make_results_xr()


def plot_results_xr() -> None:
    """
    Plot the `RESULTS_XR` object, to make a nice plot for AGU.
    """
    delta = 0.1
    points = [
        0.5 - 3 * delta / 2,
        0.5 - delta / 2,
        0.5 + delta / 2,
        0.5 + 3 * delta / 2,
    ]
    pairs = [("EEEE", "EEEE"), ("EECE", "EE6E"), ("EEEC", "EEE6"), ("EECC", "EE66")]
    colors = ["blue", "green", "orange", "red"]
    for pair_number, pair in enumerate(pairs):
        for i in [0, 1]:
            for count, letter in enumerate(RESULTS_XR.coords["var"].values):
                print(str(letter), type(str(letter)))
                plt.plot(
                    RESULTS_XR.sel(var=str(letter), mem=pair[i]).values,
                    1 - i + points[count],
                    color=colors[pair_number],
                    marker=r"$" + letter + r"$",
                )

    plt.plot([ECMWF, ECMWF], [0, 2], colors[0], label="ECMWF/ORAS4")
    plt.plot([CMIP5_MMM, CMIP5_MMM], [1, 2], "--", color=colors[3], label="CMIP5 MMM")
    plt.plot([CMIP6_MMM, CMIP6_MMM], [0, 1], ":", color=colors[3], label="CMIP6 MMM")
    plt.plot([ECMWF, CMIP5_MMM], [2, 2], "--", color="grey")
    plt.plot([ECMWF, max(CMIP5_MMM, CMIP6_MMM)], [1, 1], "--", color="grey")
    plt.plot([ECMWF, CMIP5_MMM], [0, 0], "--", color="grey")
    plt.text(0.35, 1.85, "(a)", color="black")
    plt.text(0.35, 0.85, "(b)", color="black")
    plt.xlim([0.3, 1.1])
    plt.ylim([-0.7, 2.5])
    plt.yticks([])
    plt.xlabel(r"NINO3.4 trend 1958-2017 [$\Delta$K]")
    plt.legend(loc="lower center")
    plt.savefig("result.png")


if __name__ == "__main__":
    # python src/visualisation/arrow.py
    print(RESULTS_XR)
    plot_results_xr()
    # plot_arrow_plot_6(save_path=os.path.join(FIGURE_PATH, "mech_arrow_cmip6.pdf"))
    # plot_arrow_plot_6(save_path=os.path.join(FIGURE_PATH, "mech_arrow_cmip6.png"))
