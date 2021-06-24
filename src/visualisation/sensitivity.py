"""Find the sensitivity of the coupled omodel."""
import os
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from src.xr_utils import get_trend, sel
from src.plot_utils import add_units, cmap, get_dim, ps_defaults, axis_fomatter
from src.models.poly import plot
from src.constants import (
    SEL_DICT,
    FIGURE_PATH,
    ORIG_WANDB_DATA,
    NEW_WANDB_DATA,
    CD_LOGS,
)
from src.configs.load_config import load_config
from src.models.model_setup import ModelSetup
from src.utils import timeit
from src.wandb_utils import get_wandb_data

ps_defaults(use_tex=False, dpi=200)


@timeit
def cd_plots(show_plots: bool = False) -> None:
    """
    Generate the cd sensitivity plots.
    """

    df = get_wandb_data(str(NEW_WANDB_DATA))

    f_df = df[3:13]
    f_df = f_df.drop(labels=[11], axis=0)

    cd_list = list()

    for coup_dict_str in f_df["config_coup"]:
        # pylint: disable=eval-used
        coup_dict = eval(coup_dict_str)
        cd_list.append(coup_dict["c_d"])

    for reg in SEL_DICT:

        reg_list = f_df["end_trend_" + reg].tolist()

        if reg == "nino1+2":
            reg_type = "parab"
        else:
            reg_type = "lin"

        plot(
            cd_list,
            reg_list,
            reg_type=reg_type,
            x_label="$C_d$ [dimensionless]",
            y_label=r"$\Delta T_s$ over " + reg + r" region [$\Delta$ K]",
            fig_path=str(FIGURE_PATH / str("cd_" + reg + ".pdf")),
            ax_format="x",
        )
        if show_plots:
            plt.show()
        else:
            plt.clf()


@timeit
def cd_heatmaps(show_plots: bool = False) -> None:
    """
    Return the cd heatmaps

    Args:
        show_plots (bool, optional): Whether or not to show plots or just clear them.
            Defaults to False.
    """
    cfg = load_config(test=False)
    name_direc_l = []
    names = os.listdir(CD_LOGS)
    names.sort()
    for name in names:
        direc = str(CD_LOGS / name / "wandb" / "latest-run" / "files")
        name_direc_l.append(
            (
                name,
                float(name[3:]) * 1e-3,
                direc,
                ModelSetup(direc, cfg, make_move=False),
            )
        )

    print(name_direc_l)

    cd_and_file = list()

    for i in range(len(name_direc_l)):
        print(name_direc_l[i][1])
        print(name_direc_l[i][3].ts_trend(5))
        print(os.path.exists(name_direc_l[i][3].ts_trend(5)))
        cd_and_file.append([name_direc_l[i][1], name_direc_l[i][3].ts_trend(5)])

    da_list = [xr.open_dataarray(cd_and_file[x][1]) for x in range(len(cd_and_file))]

    cd_ts_da = xr.concat(da_list, r"$C_d$").assign_coords(
        {r"$C_d$": [cd_and_file[x][0] for x in range(len(cd_and_file))]}
    )

    setup = ModelSetup(direc, cfg, make_move=False)
    mask = xr.open_dataset(setup.om_mask()).mask

    def clip(da: xr.DataArray):
        return add_units(sel(da).where(sel(mask) != 0.0))

    clip(cd_ts_da).plot(
        col="$C_d$",
        aspect=2,
        col_wrap=2,
        cmap=cmap("delta"),
        cbar_kwargs={
            "aspect": 50,
            "label": (
                r"$\Delta T_s$, Change in sea surface "
                + r" temperature over 58 years [$\Delta K$] "
            ),
        },
        figsize=get_dim(ratio=1),
    )

    plt.savefig(FIGURE_PATH / "cd_facetplot.png")
    plt.savefig(FIGURE_PATH / "cd_facetplot.pdf")

    if show_plots:
        plt.show()
    else:
        plt.clf()

    slope, hatch_mask = get_trend(
        clip(cd_ts_da),
        min_clim_f=False,
        output="slope",
        t_var="$C_d$",
        make_hatch_mask=True,
    )

    sens_heatmap_kwargs = dict(
        {
            "label": r"$\Delta T_s / \Delta C_d$ [$\Delta K$]",
            # "aspect": 35,
            "format": axis_fomatter(),
        }
    )

    add_units(slope).plot(
        cmap=cmap("delta"),
        cbar_kwargs=sens_heatmap_kwargs.copy(),
        aspect=2,
        figsize=get_dim(ratio=0.5),
    )
    add_units(hatch_mask).where(hatch_mask != 0).plot(
        add_colorbar=False,
        cmap="Greys",
        alpha=0.3,
    )
    plt.tight_layout()
    plt.title("")
    plt.savefig(FIGURE_PATH / "ts_cd_sensitivity_hatched.png")
    plt.savefig(FIGURE_PATH / "ts_cd_sensitivity_hatched.pdf")
    if show_plots:
        plt.show()
    else:
        plt.clf()

    add_units(slope).plot(
        cmap=cmap("delta"),
        cbar_kwargs=sens_heatmap_kwargs.copy(),
        aspect=2,
        figsize=get_dim(ratio=0.5),
    )
    plt.title("")
    plt.tight_layout()
    plt.savefig(FIGURE_PATH / "ts_cd_sensitivity.png")
    plt.savefig(FIGURE_PATH / "ts_cd_sensitivity.pdf")

    if show_plots:
        plt.show()
    else:
        plt.clf()


@timeit
def nummode_plots(show_plots: bool = False) -> None:
    """
    Make a plot of the number of modes vs. time taken to run the ocean section.
    """
    df = pd.read_csv(ORIG_WANDB_DATA)

    uncoup_df = df[13:]
    uncoup_df = uncoup_df[uncoup_df["state"] == "finished"]
    nummode_list = list()

    for coup_dict_str in uncoup_df[uncoup_df["state"] == "finished"]["config_oc"]:
        # pylint: disable=eval-used
        coup_dict = eval(coup_dict_str)
        nummode_list.append(coup_dict["nummode"])

    ocean_run_list = uncoup_df[uncoup_df["state"] == "finished"][
        "end_ocean_run"
    ].tolist()

    plot(
        nummode_list,
        ocean_run_list,
        x_label="Number of modes in the ocean",
        y_label="Ocean model run time [s]",
        fig_path=str(FIGURE_PATH / "nummode_times.pdf"),
    )

    if show_plots:
        plt.show()
    else:
        plt.clf()


if __name__ == "__main__":
    # python src/visualisation/sensitivity.py
    cd_plots()
    nummode_plots()
    cd_heatmaps()
