"""Program to automate testing output fields against the paper.

Perhaps this module should be renamed 'comp_v_seager19'.
"""
from typing import Union
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from typing import List
from src.models.model_setup import ModelSetup
from src.xr_utils import open_dataset, get_trend, clip, can_coords, sel
from src.wandb_utils import setup_from_name
from src.utils import get_default_setup
from src.configs.load_config import load_config
from src.plot_utils import add_units, cmap, get_dim, label_subplots, ps_defaults
from src.constants import UC_LOGS, FIGURE_DATA_PATH, FIGURE_PATH
from src.visualisation.quiver import pqp_part


def comp_plot(
    ours: xr.DataArray,
    papers: xr.DataArray,
    default_cmap="delta",
    diff_cmap="delta",
    **kwargs
) -> None:
    """
    A comparison plot for two scalar fields.

    Args:
        ours (xr.DataArray): The output of the model.
        papers (xr.DataArray): The paper's values.
        default_cmap (str, optional): The colormap for the field.
            Defaults to "delta".
        diff_cmap (str, optional): The colormap for the difference
            between the two fields. Defaults to "delta".
    """
    ours, papers = add_units(ours), add_units(papers)
    _, axs = plt.subplots(4, figsize=get_dim(ratio=0.3 * 4), sharex=True)
    ours.plot(ax=axs[0], cmap=cmap(default_cmap), **kwargs)
    axs[0].set_xlabel("")
    papers.plot(ax=axs[1], cmap=cmap(default_cmap), **kwargs)
    axs[1].set_xlabel("")
    (ours - papers).plot(ax=axs[2], cmap=cmap(diff_cmap), **kwargs)
    axs[2].set_xlabel("")
    np.abs((ours - papers) / papers).plot(ax=axs[3], vmin=0, vmax=1, cmap=cmap("sst"))
    label_subplots(axs, y_pos=1.05, x_pos=-0.1)
    plt.tight_layout()


def comp_prcp_quiver_plot(
    ours: xr.Dataset, theirs: xr.Dataset, vmin=-5e-5, vmax=5e-5, x_pos=0.73, y_pos=-0.15
) -> None:
    """
    Compare the precipitation and windspeeds.

    Args:
        ours (xr.Dataset): Our dataset.
        theirs (xr.Dataset): Their dataset.
    """
    ps_defaults(use_tex=False, dpi=200)
    _, axs = plt.subplots(3, 1, figsize=get_dim(ratio=0.3 * 3), sharex=True)
    pqp_part(axs[0], ours, x_pos=x_pos, y_pos=y_pos, vmin=vmin, vmax=vmax)
    axs[0].set_xlabel("")
    pqp_part(axs[1], theirs, x_pos=x_pos, y_pos=y_pos, vmin=vmin, vmax=vmax)
    axs[1].set_xlabel("")
    diff = ours.copy()
    diff["utrend"] = ours["utrend"] - theirs["utrend"]
    diff["vtrend"] = ours["vtrend"] - theirs["vtrend"]
    diff["PRtrend"] = ours["PRtrend"] - theirs["PRtrend"]
    pqp_part(axs[2], diff, x_pos=x_pos, y_pos=-0.35, vmin=vmin, vmax=vmax)
    label_subplots(axs, y_pos=1.05, x_pos=-0.18)
    plt.tight_layout()


def return_var_list(num: Union[int, str]) -> List[str]:
    """
    Get a list of the variables from each figure.

    Args:
        num Union[int, str]: The figure number.
            Example input: int(4) or "2a".

    Returns:
        List[str]: A list of the variable names.
    """
    var_list = []
    for var in xr.open_dataset(FIGURE_DATA_PATH):
        if "Fig_" + str(num) in var:
            var_list.append(var)
    return var_list


def return_figure_ds(num: str) -> xr.Dataset:
    """
    Get the figure dataset.

    Args:
        num (str): The figure number e.g. "2c".

    Returns:
        xr.Dataset: the dataset with the standard names.
    """
    ps_defaults(use_tex=False, dpi=200)
    fig_data = xr.open_dataset(FIGURE_DATA_PATH)
    r_dict = {}
    for i in fig_data[return_var_list(num)]:
        r_dict[i] = i.split(".")[-1]

    return fig_data[return_var_list(num)].rename(r_dict)


def comp_uc_oc(setup: ModelSetup, panel="d", show_plots: bool = False) -> None:
    """
    Test to see if panel 1d is replicated.

    Args:
        setup (ModelSetup): The setup object.
        panel (str, optional): Which panel to test aginst. Defaults to "d".
    """
    ps_defaults(use_tex=False, dpi=200)
    fig_data = xr.open_dataset(FIGURE_DATA_PATH)
    uc_oc = xr.open_dataset(setup.om_run2f_nc(), decode_times=False)
    uc_oc_dt = add_units(get_trend(clip(can_coords(uc_oc.SST_SST))).isel(Z=0).drop("Z"))
    uc_oc_dt.attrs["units"] = r"$\Delta$ K"
    uc_oc_dt.attrs["long_name"] = r"$\Delta$ SST"
    ddata = add_units(
        sel(
            can_coords(fig_data["ForcedOceanModel.sst-trend-Fig_1" + panel + ".nc.SST"])
        )
    )
    ddata = ddata.where(ddata != 0.0).rename(r"$\Delta$ SST")
    ddata.attrs["units"] = r"$\Delta$ K"
    ddata.attrs["long_name"] = r"$\Delta$ SST"
    comp_plot(add_units(uc_oc_dt.interp_like(ddata)), ddata)

    if show_plots:
        plt.show()
    else:
        plt.clf()


def comp_uc_atm(setup: ModelSetup, panel="d", show_plots: bool = False) -> None:
    """
    Test to see if panel 2d is right.

    Args:
        setup (ModelSetup): The path object.
        panel (str, optional): Which panel to test against. Defaults to "d".
    """
    ps_defaults(use_tex=False, dpi=200)
    uc_atm = open_dataset(setup.tcam_output())
    ads = return_figure_ds("2" + panel)
    comp_prcp_quiver_plot(uc_atm, ads)
    plt.savefig("example.png")

    if show_plots:
        plt.show()
    else:
        plt.clf()


def comp_atm_prwnd(setup: ModelSetup, num: str, show_plots: bool = False) -> str:
    """
    Test to see if atm is right.

    Args:
        setup (ModelSetup): The path object.
        panel (str): Which panel to test against. E.g. 2d.
    """
    ps_defaults(use_tex=False, dpi=200)
    uc_atm = open_dataset(setup.tcam_output())
    ads = return_figure_ds(num)
    comp_prcp_quiver_plot(uc_atm, ads)
    plt.savefig(setup.rep_plot(num, "_prwnd"))

    if show_plots:
        plt.show()
    else:
        plt.clf()

    return setup.rep_plot(num, "_prwnd")


def comp_oc_sst(
    setup: ModelSetup, num: str, show_plots: bool = False, var="tstrend"
) -> str:
    """
    Compare the sea surface temperature trend of the final model iteration.

    Args:
        setup (ModelSetup): The setup.
        num (str): The number e.g. "2b".
    """
    if "1" in num:
        var = "SST"
    ps_defaults(use_tex=False, dpi=200)
    uc_oc = xr.open_dataset(setup.om_run2f_nc(), decode_times=False)
    oc_dt = add_units(get_trend(clip(can_coords(uc_oc.SST_SST))).isel(Z=0).drop("Z"))
    oc_dt.attrs["units"] = r"$\Delta$ K"
    oc_dt.attrs["long_name"] = r"$\Delta$ SST"
    ds = return_figure_ds(num)
    ddata = add_units(sel(can_coords(ds[var])))
    ddata = ddata.where(ddata != 0.0).rename(r"$\Delta$ SST")
    ddata.attrs["units"] = r"$\Delta$ K"
    ddata.attrs["long_name"] = r"$\Delta$ SST"
    comp_plot(add_units(oc_dt.interp_like(ddata)), ddata, vmin=-2, vmax=2)
    plt.savefig(setup.rep_plot(num, "_sst"))

    if show_plots:
        plt.show()
    else:
        plt.clf()

    return setup.rep_plot(num, "_sst")


def comp_oc_htherm(setup: ModelSetup, num: str, show_plots: bool = False) -> str:
    """
    Compare the sea surface temperature trend of the final model iteration.

    Args:
        setup (ModelSetup): The setup.
        num (str): The number e.g. "4b".
    """
    ps_defaults(use_tex=False, dpi=200)
    uc_oc = xr.open_dataset(setup.om_run2f_nc(), decode_times=False)
    oc_dt = add_units(
        get_trend(clip(can_coords(uc_oc.TDEEP_HMODEL))).isel(Z=0).drop("Z")
    )
    oc_dt.attrs["units"] = r"$\Delta$ m"
    oc_dt.attrs["long_name"] = r"$\Delta$ $H_T$"
    ds = return_figure_ds(num)
    ddata = add_units(sel(can_coords(ds["HTHERM"])))
    ddata = ddata.where(ddata != 0.0)  # .rename(r"$\Delta$ SST")
    ddata.attrs["units"] = r"$\Delta$ m"
    ddata.attrs["long_name"] = r"$\Delta$ $H_T$"
    comp_plot(add_units(oc_dt.interp_like(ddata)), ddata, vmin=-30, vmax=30)
    plt.savefig(setup.rep_plot(num, "_htherm"))

    if show_plots:
        plt.show()
    else:
        plt.clf()

    return setup.rep_plot(num, "_htherm")


if __name__ == "__main__":
    # python src/visualisation/comp.py
    # python src/visualisation/comp_v_seager19.py
    import os

    plot_dir = "/gws/nopw/j04/ai4er/users/sdat2/sensitivity/k_days_logs/k_days_10/plots"

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    uncoupled_run_dir = str(UC_LOGS / "it_1")
    cfg = load_config(test=False)
    uncoup_setup = ModelSetup(uncoupled_run_dir, cfg, make_move=False)
    coup_setup = get_default_setup()
    # comp_uc_atm(uncoup_setup)
    comp_oc_sst(coup_setup, "3")
    comp_atm_prwnd(coup_setup, "3")
    comp_oc_htherm(coup_setup, "4b")
    # comp_uc_oc(uncoup_setup)


def field_plot() -> None:
    """Plot the different input fields that are varied."""
    _, axs = plt.subplots(4, 1, figsize=get_dim(ratio=0.3 * 4), sharex=True)
    clip(
        can_coords(xr.open_dataarray(setup_from_name("ECMWF_coup").clim60_name(3)))
    ).plot(
        ax=axs[0],
        cmap=cmap("sst"),
        vmin=30,
        vmax=90,
        cbar_kwargs={"label": r"ECMWF   $\bar{r}$  [%]"},
    )
    axs[0].set_xlabel("")
    clip(
        can_coords(xr.open_dataarray(setup_from_name("C_RH_coup").clim60_name(3)))
    ).plot(
        ax=axs[1],
        cmap=cmap("sst"),
        vmin=30,
        vmax=90,
        cbar_kwargs={"label": r"CMIP5   $\bar{r}$  [%]"},
    )
    axs[1].set_xlabel("")
    clip(can_coords(return_figure_ds("5a")["wnspClim"])).plot(
        ax=axs[2],
        cmap=cmap("sst"),
        vmin=4,
        vmax=8,
        cbar_kwargs={"label": r"ECMWF   $\bar{W}$  [m s$^{-1}$]"},
    )
    axs[2].set_xlabel("")
    clip(can_coords(return_figure_ds("5d")["wnspClim"])).plot(
        ax=axs[3],
        cmap=cmap("sst"),
        vmin=4,
        vmax=8,
        cbar_kwargs={"label": r"CMIP5   $\bar{W}$  [m s$^{-1}$]"},
    )
    label_subplots(axs, y_pos=1.05, x_pos=-0.1)
    plt.tight_layout()
    plt.savefig(str(FIGURE_PATH / "fields.pdf"))
    plt.savefig(str(FIGURE_PATH / "fields.png"))
