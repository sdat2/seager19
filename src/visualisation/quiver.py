"""Quiver plots."""
from scipy.interpolate import interp2d
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
from src.utils import timeit, get_default_setup
from src.models.model_setup import ModelSetup
from src.xr_utils import can_coords, clip
from src.plot_utils import cmap, add_units, get_dim, ps_defaults
from src.constants import FIGURE_PATH


@timeit
def prcp_quiver_plot(
    setup: ModelSetup = get_default_setup(),
    show_plots: bool = False,
    save_path=str(FIGURE_PATH / "uv_prcp.png"),
) -> None:
    """
    prcp wind trends quiver plot.

    Args:
        setup (ModelSetup, optional): setup object. Defaults to get_default_setup().
        show_plots (bool, optional): Whether to show plots. Defaults to False.
        save_path ([type], optional): Path to save fig to. Defaults to
            str(FIGURE_PATH / "uv_prcp.png").
    """
    ps_defaults(use_tex=False, dpi=200)
    _, ax = plt.subplots(1, figsize=get_dim(ratio=0.5))
    ads = xr.open_dataset(setup.tcam_output())
    pqp_part(ax, ads)
    # plt.title("$\Delta$ precipitation and wind velocities [m s$^{-1}$]")
    plt.tight_layout()
    plt.savefig(save_path)
    if show_plots:
        plt.show()
    else:
        plt.clf()


def pqp_part(ax: matplotlib.axes.Axes, ads: xr.Dataset) -> None:
    """
    Plot a panel of pqp figure.

    Args:
        ax (matplotlib.axes.Axes): axes to plot on.
        ads (xr.Dataset): Standard atmos dataset.
    """

    new_x = list(range(-100, 291, 5))
    new_y = list(range(-30, 31, 5))
    fvtrend = interp2d(ads.X, ads.Yv, ads.vtrend, kind="linear")
    futrend = interp2d(ads.X, ads.Yu, ads.utrend, kind="linear")
    new_ds = xr.Dataset(
        {
            "X": ("X", new_x),
            "Y": ("Y", new_y),
        }
    )
    new_ds.X.attrs = [("units", "degree_east")]
    new_ds.Y.attrs = [("units", "degree_north")]
    new_ds["utrend"] = (["Y", "X"], futrend(new_x, new_y))
    new_ds["vtrend"] = (["Y", "X"], fvtrend(new_x, new_y))
    clip(can_coords(ads.PRtrend)).plot(
        ax=ax,
        cmap=cmap("ranom"),
        cbar_kwargs={"label": "Precipitation [m s$^{-1}$]"},
    )
    quiver = add_units(new_ds).plot.quiver(
        ax=ax, x="X", y="Y", u="utrend", v="vtrend"
    )  # , normalize=matplotlib.colors.Normalize(vmin=0.01, vmax=1))#, scale=30)
    _ = plt.quiverkey(
        quiver,
        0.85,
        0.18,
        1,
        r"  $1$ m s$^{-1}$" + "\n" + r" $\Delta \vec{u}$",
        labelpos="S",
        coordinates="figure",
    )


if __name__ == "__main__":
    # python src/visualisation/quiver.py
    prcp_quiver_plot()
