"""Program to do; plot trends etc."""
import xarray as xr
import matplotlib.pyplot as plt
from src.utils import timeit
from src.models.model_setup import ModelSetup
from src.xr_utils import can_coords, clip, get_trend
from src.plot_utils import cmap, get_dim, plot_defaults, label_subplots
from src.utils import get_default_setup
from src.constants import FIGURE_PATH


@timeit
def up_therm_qnet(
    setup: ModelSetup = get_default_setup(),
    show_plots: bool = False,
    save_path: str =str(FIGURE_PATH / "tuq_trends.png"),
) -> None:
    """
    up trend qnet trend.

    Args:
        setup (ModelSetup, optional): setup object. Defaults to get_default_setup().
        show_plots (bool, optional): Whether to show the plots. Defaults to False.
        save_path (str, optional): path to save fig to.
            Defaults to str(FIGURE_PATH / "tuq_trends.png").
    """
    plot_defaults(use_tex=False, dpi=200)
    ds = xr.open_dataset(setup.om_run2f_nc(), decode_times=False)
    _, axs = plt.subplots(3, 1, figsize=get_dim(ratio=0.28 * 3.8))
    clip(get_trend(can_coords(ds.TDEEP_HMODEL))).plot(ax=axs[0], cmap=cmap("delta"))
    axs[0].set_title(r"$\Delta$ Thermocline depth over 58 years [m]")
    axs[0].set_xlabel("")
    clip(get_trend(can_coords(ds.SST_QNET))).plot(ax=axs[2], cmap=cmap("delta"))
    axs[2].set_title(r"$\Delta$ Net surface heat flux over 58 years [W m$^{-2}$]")
    clip(get_trend(can_coords(ds.SST_W1))).plot(ax=axs[1], cmap=cmap("delta"))
    axs[1].set_title(r"$\Delta$ Upwelling over 58 years [m s$^{-1}$]")
    axs[1].set_xlabel("")
    label_subplots(axs, y_pos=1.25, x_pos=-0.15)
    plt.tight_layout()
    plt.savefig(save_path)
    if show_plots:
        plt.show()
    else:
        plt.clf()


if __name__ == "__main__":
    # python src/visualisation/trends.py
    up_therm_qnet()
