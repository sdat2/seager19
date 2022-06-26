"""Plot and process biases in trends and means."""
import os
import pathlib
from typing import Optional, Union
import pandas as pd
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
from src.constants import (
    atmos_input_file_path,
    MODEL_NAMES,
    DATA_PATH,
    CMIP6_TS_PATH,
)
from src.plot_utils import plot_defaults, label_subplots, cmap, plot_defaults, set_dim
from src.xr_utils import (
    can_coords,
    sel,
    get_trend,
    spatial_mean,
    clip,
)
from src.visualisation.nino import plot_nino_box

# from src.data_loading.pangeo import load_mfds
# import numpy as np


def atmos_bias(var: str = "ts", model: str = "S", ending: str = "clim") -> xr.DataArray:
    """
    Atmospheric bias.

    Args:
        var (str, optional): Variable. Defaults to "ts".
        model (str, optional): Model. Defaults to "S".
        ending (str, optional): Ending. Defaults to "clim".

    Returns:
        xr.DataArray:
    """
    ecmwf = can_coords(
        xr.open_dataarray(
            atmos_input_file_path(var=var, model="E", ending=ending), decode_times=False
        )
    )
    data = can_coords(
        xr.open_dataarray(
            atmos_input_file_path(var=var, model=model, ending=ending),
            decode_times=False,
        )
    )
    bias = data - ecmwf
    bias.attrs["units"] = data.attrs["units"]
    bias.attrs["long_name"] = data.attrs["long_name"]
    bias.attrs["description"] = str(
        var + " bias " + MODEL_NAMES[model] + " above " + MODEL_NAMES["E"]
    )
    return bias


def get_trends(file_path: Union[str, pathlib.Path] = CMIP6_TS_PATH) -> xr.DataArray:
    """
    Get the 60 year trends over the nino3.4 region
    for differnt members of the nino3.4 ensemble.

    Args:
        file_path (Union[str, pathlib.Path], optional): Defaults to CMIP6_TS_PATH.

    Returns:
        xr.DataArray: xr trends.
    """
    file_list = os.listdir(file_path)
    tr_list = []
    fin_list = []
    for start, fin in [(x - 59, x) for x in range(2007, 2017)]:
        tr_list.append([])
        for i in range(len(file_list)):
            ex = can_coords(
                xr.open_dataarray(os.path.join(file_path, file_list[i])).sel(
                    T=slice(str(start), str(fin))
                )
            )
            tr_list[-1].append(
                get_trend(spatial_mean(sel(ex, reg="nino3.4")), min_clim_f=True)
            )
        fin_list.append(fin)
    return xr.DataArray(
        tr_list,
        dims=["end_year", "member"],
        coords=dict(end_year=fin_list, member=os.listdir(CMIP6_TS_PATH)),
        attrs=dict(units="K", long_name="60 year trend in nino3.4"),
    )


def trends_from_da(
    da: xr.DataArray, beginning: int = 2007, finish: int = 2017
) -> xr.DataArray:
    """
    Get 60 year trends over a few years.

    Args:
        da (xr.DataArray): Temperature dataarray.
        beginning (int): Defaults to 2007.
        finish (int): Defaults to 2017.

    Returns:
        xr.DataArray: Rise dataarray.
    """
    da_tr_list = []
    for start, end in [(str(x - 59), str(x)) for x in range(beginning, finish + 1)]:
        da_tr_list.append(
            get_trend(da.sel(T=slice(start, end)), keep_ds=True).expand_dims(
                {"end_year": [int(end)]}
            )
        )
    da_tr = xr.merge(da_tr_list)
    da_tr.end_year.attrs["long_name"] = "End year"
    da_tr.rise.attrs["units"] = "K"
    da_tr.rise.attrs["long_name"] = "Trend over 60 years up to end year"
    return da_tr.rise


def trends_from_csv(
    csv_path: Union[pathlib.Path, str] = DATA_PATH / "mmm-trends.csv"
) -> xr.DataArray:
    """
    Get trends csv.

    Args:
        csv_path (pathlib.Path, optional): Defaults to DATA_PATH / "mmm-trends.csv".

    Returns:
        xr.DataArray:
    """
    df = pd.read_csv(csv_path, index_col=0)
    print(df)
    da = xr.DataArray(
        data=df.values,
        dims=["end_year", "source"],
        coords=dict(
            end_year=(["end_year"], df.index.values),
            source=(["source"], df.columns.values),
        ),
        attrs=dict(long_name="nino3.4 60 year trend up to end year", units="K",),
    )
    da["end_year"].attrs["long_name"] = "End year"
    da["source"].attrs["long_name"] = "Multi-model-mean"
    return da


def plot_bias(
    ax: matplotlib.axes.Axes,
    var: str = "ts",
    ending: str = "trend",
    model: str = "S",
    cmap_loc: any = cmap("delta"),
    vmin: float = -1.5,
    vmax: float = 1.5,
    long_name: Optional[str] = None,
    units: Optional[str] = None,
):
    """
    Plot an individual bias.

    Args:
        ax (matplotlib.axes.Axes): axes.
        var (str, optional): Variable. Defaults to "ts".
        ending (str, optional): Ending. Defaults to "trend".
        model (str, optional): Model. Defaults to "S".
        cmap_loc (any, optional): cmap. Defaults to cmap("delta").
        vmin (float, optional): vmin. Defaults to -1.5.
        vmax (float, optional): vmax. Defaults to 1.5.
        long_name (str, optional): Long name - label for colormap. Defaults to None.
        units (str, optional): Units for colormap. Defaults to None.
    """
    bias = atmos_bias(var=var, ending=ending, model=model)
    clip_bias = clip(bias)
    if long_name is not None:
        clip_bias.attrs["long_name"] = long_name
    if units is not None:
        clip_bias.attrs["units"] = units
    clip_bias.plot(ax=ax, cmap=cmap_loc, vmin=vmin, vmax=vmax)
    plot_nino_box(ax, reg="nino3.4", color=None)


def calculate_cmip6_ts_trends():
    print("ok")


def line_plot(da: xr.DataArray, ax: matplotlib.axes.Axes, **kwargs) -> any:
    print("ok")
    da.plot.line(ax=ax, **kwargs)


def temperature_line_plot(ax: matplotlib.axes.Axes) -> None:
    """
    Temperature line plot.

    Args:
        ax (matplotlib.axes.Axes): which axis to add to.
    """
    da = trends_from_csv()
    reanal = ["NCEP NCAR", "ERSSTv5", "HadlSST"]
    mmm = ["CMIP5 MMM", "LENS MMM"]
    reanal_da = da.sel(source=reanal)
    reanal_da["source"].attrs["long_name"] = "Reanalysis product"
    reanal_da.plot.line(ax=ax, hue="source", linestyle="dashed", add_legend=False)
    da.sel(source=mmm).plot.line(ax=ax, hue="source", add_legend=False)
    # ax.legend(labels, bbox_to_anchor=(-1.15, 2.4, -1., .102),
    # loc='lower left', ncol = 6)
    ax.set_xlim([2008, 2017])


def multi_bias_plot(model: str = "S", vertical=True) -> None:
    """
    Plot the biases plot.

    Args:
        model (str, optional): Model char. Defaults to "S" for CMIP6.

    Example:
        plot cmip6::

            from src.visualisation.bias import multi_bias_plot

            multi_bias_plot("S")
    """
    plt.clf()
    if vertical:
        fig, axs = plt.subplots(4, 1)
    else:
        fig, axs = plt.subplots(2, 2)
    set_dim(fig, ratio=1.5)
    plot_defaults()
    axs = axs.ravel()
    temperature_line_plot(axs[0])
    i = 0
    cmap_axs = axs[1:4]
    for var, ending, lim, long_name, units in [
        ("ts", "trend", 1.5, "Temperature rise bias", r"$\Delta$K"),
        ("rh", "clim60", 10, "Relative humidity bias", "%"),
        ("sfcWind", "clim", 1.5, "Surface wind bias", r"m s$^{-1}$"),
    ]:
        try:
            plot_bias(
                cmap_axs[i],
                var=var,
                ending=ending,
                model=model,
                cmap_loc=cmap("delta"),
                vmin=-lim,
                vmax=lim,
                long_name=long_name,
                units=units,
            )
        # pylint: disable=broad-except
        except Exception as e:
            print(e)
        i += 1
    if vertical:
        axs[1].set_xlabel("")
        axs[1].set_xticks([])
        axs[2].set_xticks([])
        axs[2].set_xlabel("")
    else:
        axs[1].set_xticks([])
        axs[1].set_ylabel("")
        axs[1].set_xlabel("")
        axs[1].set_yticks([])
        axs[3].set_yticks([])
        axs[3].set_ylabel("")
    label_subplots(axs)
    plt.tight_layout()


if __name__ == "__main__":
    # python src/visualisation/bias.py
    print("ok")
