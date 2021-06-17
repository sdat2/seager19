"""Sets a consistent plotting settings across the project.

Example:
    Usage with simple plots::

        from src.plot_utils import (
            ps_defaults,
            label_subplots,
            get_dim,
            set_dim,
            PALETTE,
            STD_CLR_LIST,
            CAM_BLUE,
            BRICK_RED,
            OX_BLUE,
        )

        ps_defaults(use_tex=True)

        # ---- example set of graphs ---

        import numpy as np
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(2, 2)

        x = np.linspace(0, np.pi, num=100)
        axs[0, 0].plot(x, np.sin(x), color=STD_CLR_LIST[0])
        axs[0, 1].plot(x, np.cos(x), color=STD_CLR_LIST[1])
        axs[1, 0].plot(x, np.sinc(x), color=STD_CLR_LIST[2])
        axs[1, 1].plot(x, np.abs(x), color=STD_CLR_LIST[3])

        # set size
        set_dim(fig, fraction_of_line_width=1, ratio=(5 ** 0.5 - 1) / 2)

        # label subplots
        label_subplots(axs, start_from=0, fontsize=10)

"""
from typing import Sequence, Tuple, Union
import numpy as np
from sys import platform
import itertools
from distutils.spawn import find_executable
import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import seaborn as sns
import cartopy.crs as ccrs
import cftime
import cmocean
from src.constants import REPORT_WIDTH, DATE_TITLE_FORMAT
from src.utils import timeit


def label_subplots(
    axs: Sequence[matplotlib.axes.Axes],
    labels: Sequence[str] = [chr(ord("`") + z) for z in range(1, 27)],
    start_from: int = 0,
    fontsize: int = 10,
    x_pos: float = 0.02,
    y_pos: float = 0.95,
) -> None:
    """Adds e.g. (a), (b), (c) at the top left of each subplot panel.

    Labelling order achieved through ravelling the input `list` or `np.array`.

    Args:
        axs (Sequence[matplotlib.axes.Axes]): `list` or `np.array` of
            `matplotlib.axes.Axes`.
        labels (Sequence[str]): A sequence of labels for the subplots.
        start_from (int, optional): skips first `start_from` labels. Defaults to 0.
        fontsize (int, optional): Font size for labels. Defaults to 10.
        x_pos (float, optional): Relative x position of labels. Defaults to 0.02.
        y_pos (float, optional): Relative y position of labels. Defaults to 0.95.

    Returns:
        void; alters the `matplotlib.axes.Axes` objects

    Example:
        Here is an example of using this function::

            >>> from src.plot_utis import label_subplots
            >>> label_subplots(axs, start_from=0, fontsize=10)

    """
    if isinstance(axs, list):
        axs = np.asarray(axs)
    assert len(axs.ravel()) + start_from <= len(labels)
    subset_labels = []
    for i in range(len(axs.ravel())):
        subset_labels.append(labels[i + start_from])
    for i, label in enumerate(subset_labels):
        axs.ravel()[i].text(
            x_pos,
            y_pos,
            str("(" + label + ")"),
            color="black",
            transform=axs.ravel()[i].transAxes,
            fontsize=fontsize,
            fontweight="bold",
            va="top",
        )


def get_dim(
    width: float = REPORT_WIDTH,
    fraction_of_line_width: float = 1,
    ratio: float = (5 ** 0.5 - 1) / 2,
) -> Tuple[float, float]:
    """Return figure height, width in inches to avoid scaling in latex.

    Default width is `src.constants.REPORT_WIDTH`.
    Default ratio is golden ratio, with figure occupying full page width.

    Args:
        width (float, optional): Textwidth of the report to make fontsizes match.
            Defaults to `src.constants.REPORT_WIDTH`.
        fraction_of_line_width (float, optional): Fraction of the document width
            which you wish the figure to occupy.  Defaults to 1.
        ratio (float, optional): Fraction of figure width that the figure height
            should be. Defaults to (5 ** 0.5 - 1)/2.

    Returns:
        fig_dim (tuple):
            Dimensions of figure in inches

    Example:
        Here is an example of using this function::

            >>> from src.plot_utils import get_dim
            >>> dim_tuple = get_dim(fraction_of_line_width=1, ratio=(5 ** 0.5 - 1) / 2)

    """

    # Width of figure
    fig_width_pt = width * fraction_of_line_width

    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * ratio

    return (fig_width_in, fig_height_in)


def set_dim(
    fig: matplotlib.figure.Figure,
    width: float = REPORT_WIDTH,
    fraction_of_line_width: float = 1,
    ratio: float = (5 ** 0.5 - 1) / 2,
) -> None:
    """Set aesthetic figure dimensions to avoid scaling in latex.

    Default width is `src.constants.REPORT_WIDTH`.
    Default ratio is golden ratio, with figure occupying full page width.

    Args:
        fig (matplotlib.figure.Figure): Figure object to resize.
        width (float): Textwidth of the report to make fontsizes match.
            Defaults to `src.constants.REPORT_WIDTH`.
        fraction_of_line_width (float, optional): Fraction of the document width
            which you wish the figure to occupy.  Defaults to 1.
        ratio (float, optional): Fraction of figure width that the figure height
            should be. Defaults to (5 ** 0.5 - 1)/2.

    Returns:
        void; alters current figure to have the desired dimensions

    Example:
        Here is an example of using this function::

            >>> from src.plot_utils import set_dim
            >>> set_dim(fig, fraction_of_line_width=1, ratio=(5 ** 0.5 - 1) / 2)

    """
    fig.set_size_inches(
        get_dim(width=width, fraction_of_line_width=fraction_of_line_width, ratio=ratio)
    )


def ps_defaults(use_tex: bool = True, dpi: int = 600) -> None:
    """Apply plotting style to produce nice looking figures.

    Call this at the start of a script which uses `matplotlib`.
    Can enable `matplotlib` LaTeX backend if it is available.

    TODO: There seems to currently be a bug where markers
         don't appear with 'plt.plot([], [], "x").
         Need to figure out how I have caused this.

    Args:
        use_tex (bool, optional): Whether or not to use latex matplotlib backend.
            Defaults to True.
        dpi (int, optional): Which dpi to set for the figures.
            Defaults to 600 dpi (high quality). 150 dpi probably
            fine for notebooks. Largest dpi needed for presentations.

    Examples:
        Basic setting the plotting defaults::

            >>> from src.plot_utils import ps_defaults
            >>> ps_defaults()

        Setting defaults for a jupyter notebook::

            >>> from src.plot_utils import ps_defaults
            >>> ps_defaults(use_tex=False, dpi=150)

    """
    if platform == "darwin":
        matplotlib.use("TkAgg")

    p_general = {
        "font.family": "STIXGeneral",  # Nice alternative font.
        # "font.family": "serif",
        # "font.serif": [],
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 10,
        "font.size": 10,
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        # currently I don't use these autofomatters, but hopefully they
        # are sensible choices to have made.
        "date.autoformatter.year": "%Y",
        "date.autoformatter.month": "%Y-%m",
        "date.autoformatter.day": "%Y-%m-%d",
        "date.autoformatter.hour": "%m-%d %H",
        "date.autoformatter.minute": "%Y-%m-%d %H:%M:%S",
        "date.autoformatter.second": "%H:%M:%S",
        "date.autoformatter.microsecond": "%M:%S.%f",
        # Set the font for maths
        "mathtext.fontset": "cm",
        # "font.sans-serif": ["DejaVu Sans"],  # gets rid of error messages
        # "font.monospace": [],
        "figure.figsize": get_dim(),
        "lines.linewidth": 1.0,
        "scatter.marker": "x",
        "image.cmap": "viridis",
    }
    matplotlib.rcParams.update(p_general)
    matplotlib.style.use("seaborn-colorblind")

    if use_tex and find_executable("latex"):
        p_setting = {
            "pgf.texsystem": "pdflatex",
            "text.usetex": True,
            "pgf.preamble": (
                r"\usepackage[utf8x]{inputenc} \usepackage[T1]{fontenc}"
                + r"\usepackage[separate -uncertainty=true]{siunitx}"
            ),
        }
    else:
        p_setting = {
            "text.usetex": False,
        }
    matplotlib.rcParams.update(p_setting)


# Standard color list
STD_CLR_LIST = [
    "#4d2923ff",
    "#494f1fff",
    "#38734bff",
    "#498489ff",
    "#8481baff",
    "#c286b2ff",
    "#d7a4a3ff",
]
_paper_colors = sns.color_palette(STD_CLR_LIST)
# Note: To inspect colors, call `sns.palplot(_paper_colors)`
PALETTE = itertools.cycle(_paper_colors)
CAM_BLUE = "#a3c1ad"
OX_BLUE = "#002147"
BRICK_RED = "#CB4154"


def default_projection() -> ccrs.CRS:
    """
    Returns default projection.

    Returns:
        ccrs.CRS: Map projection.
    """
    return ccrs.Mollweide(central_longitude=180)


@timeit
def map_setup(ax: matplotlib.axes.Axes = None) -> matplotlib.axes.Axes:
    """Apply default map (Robinson centred on the Pacific).

    Args:
        ax (matplotlib.axes.Axes): the axes to format. Defaults to None.

    Returns:
        matplotlib.axes.Axes: axes object in Robinson config.

    Example:
        When using multiple subplots::

            from src.plot_utils import map_setup
            fig, axes = plt.subplots(
               2, 2, subplot_kw={"projection": default_projection()}
            )
            for ax in axes.ravel():
                ax = map_setup(ax=ax)

    """

    if ax is None:
        ax = plt.axes(projection=default_projection())

    ax.set_global()
    ax.coastlines()
    return ax


def time_title(
    ax: matplotlib.axes.Axes,
    time: Union[np.datetime64, float, cftime.Datetime360Day],
    date_time_formatter: str = DATE_TITLE_FORMAT,
) -> None:
    """Add time title to axes.

    Args:
        ax (matplotlib.axes.Axes): axis to add title to.
        time (Union[np.datetime64, float, cftime.Datetime360Day]): time string.
        date_time_formatter (str, optional): Default is
            `src.constants.DATE_TITLE_FORMAT`.

    Example:
        Usage with an xarray.Datarray object::

            >>> from src.plot_utils import time_title
            >>> time_title(ax, xr_da.time.values[index])

    """
    if isinstance(time, np.datetime64):
        # use pandas to format time
        ax.set_title(pd.to_datetime(str(time)).strftime(date_time_formatter))
    elif isinstance(time, cftime.Datetime360Day):
        ax.set_title(time.strftime()[0:10])
    elif isinstance(time, (float, np.floating)):
        # It would be better to have this as an option
        ax.set_title("%2.1f months after 1960" % time)
    else:
        print(
            "!Warning!: input of type "
            + str(type(time))
            + " does not lead to title plotting."
        )


def cmap(variable_name: str) -> matplotlib.colors.LinearSegmentedColormap:
    """Get cmap from a variable name string.

    Args:
        variable_name (str): name of variable to give colormap.

    Returns:
        matplotlib.colors.LinearSegmentedColormap: sensible colormap

    Example:
        Usage example for sea surface temperature::

            from src.plot_utils import cmap
            cmap_t = cmap("sst")

    """

    # make the function case insensitive
    variable_name = variable_name.lower()

    # collate the variables into a smaller number
    map_d = {
        "u": "v",
        "v": "v",
        "sst": "sst",
        "salt": "sss",
        "sss": "haline",
        "haline": "haline",
        "delta": "delta",
    }

    # map to cmocean colormaps
    cmap_map_d = {
        "sst": cmocean.cm.thermal,
        "haline": cmocean.cm.haline,
        "v": cmocean.cm.speed,
        "delta": cmocean.cm.balance,
    }

    # get cmap_t
    cmap_t = cmap_map_d[map_d[variable_name]]

    # make the map green-ish for nan values
    cmap_t.set_bad(color="#15b01a")

    return cmap_t


def add_units(
    xr_obj: Union[xr.DataArray, xr.Dataset], x_val="X", y_val="Y"
) -> Union[xr.DataArray, xr.Dataset]:
    """
    Adding good units to make axes plottable.

    Currently only for lat, lon axes, but could be improved to
    add degrees celsius and so on.

    Fails softly.

    Args:
        xr_da (Union[xr.DataArray, xr.Dataset]: Initial datarray/datset
            (potentially with units for axes).

    Returns:
        Union[xr.DataArray, xr.Dataset]: Datarray/Dataset with correct
            units/names for plotting. Assuming that you've given the
            correct x_val and y_val for the object.
    """
    if x_val in xr_obj.coords:
        xr_obj.coords[x_val].attrs["units"] = r"$^{\circ}$E"
        xr_obj.coords[x_val].attrs["long_name"] = "Longitude"
    if y_val in xr_obj.coords:
        xr_obj.coords[y_val].attrs["units"] = r"$^{\circ}$N"
        xr_obj.coords[y_val].attrs["long_name"] = "Latitude"
    return xr_obj
