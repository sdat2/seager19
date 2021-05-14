"""Test the plot settings in `src.plot_settings`."""
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from src.plot_settings import (
    ps_defaults,
    label_subplots,
    set_dim,
    STD_CLR_LIST,
    BRICK_RED,
    OX_BLUE,
    map_setup,
    time_title,
    default_projection,
)
from src.constants import PROJECT_PATH


def test_plot() -> None:
    """Function to make example plot."""

    for use_tex in [False, True]:

        ps_defaults(use_tex=use_tex)

        fig, axs = plt.subplots(2, 2)
        x = np.linspace(0, np.pi, num=100)

        # plot the different axes.
        axs[0, 0].plot(x, np.sin(x), color=STD_CLR_LIST[0])
        axs[0, 1].plot(x, np.cos(x), color=STD_CLR_LIST[1])
        axs[1, 0].plot(x, np.sinc(x), color=BRICK_RED)
        axs[1, 1].plot(x, np.abs(x), color=OX_BLUE)

        # set size
        set_dim(fig, fraction_of_line_width=1, ratio=(5 ** 0.5 - 1) / 2)

        # label subplots
        label_subplots(axs, start_from=0, fontsize=10)

        plt.savefig(str(os.path.join(PROJECT_PATH, "gifs", "example.png")), dpi=800)
        plt.clf()


def test_map_plot() -> None:
    """Tests `src.plot_settings.map_setup`.

    Test takes a long time.

    452.55s call     src/test/test_plot.py::test_map_plot

    Example:
        More elongated example for facets:
            plt.tight_layout()
            plt.savefig(str(os.path.join(PROJECT_PATH, "gifs", "map_example.png")))
            plt.clf()

            _, axes = plt.subplots(2, 1, subplot_kw={"projection": default_projection()})

            for ax in axes.ravel():
                ax = map_setup(ax=ax)
                da.plot.imshow(
                    ax=ax, transform=default_projection(), cbar_kwargs={"shrink": 0.5}
                )
                time_title(ax, da.time.values)

            plt.tight_layout()
            plt.savefig(str(os.path.join(PROJECT_PATH, "gifs", "multi_map_example.png")))
            plt.clf()

            p = (
                xr.tutorial.open_dataset("rasm")
                .load()
                .Tair.isel(time=[0, 2])
                .plot(
                    transform=ccrs.PlateCarree(),
                    col="time",
                    subplot_kws={"projection": default_projection()},
                )
            )

            for ax in p.axes.flat:
                ax = map_setup(ax=ax)

            plt.tight_layout()
            plt.savefig(
                str(os.path.join(PROJECT_PATH, "gifs", "facet_map_example.png")),
                bbox_inches="tight",
            )
            plt.clf()

    """

    ps_defaults(use_tex=False)
    ax = map_setup()

    da = xr.tutorial.open_dataset("rasm").load().Tair.isel(time=0)
    da.plot.imshow(
        ax=ax,
        transform=default_projection(),
        cbar_kwargs={"shrink": 0.5},
    )
    time_title(ax, da.time.values)
