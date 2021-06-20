"""Find the sensitivity of the coupled omodel."""
import pandas as pd
import matplotlib.pyplot as plt
import wandb_summarizer.download
from src.models.poly import plot
from src.constants import SEL_DICT, FIGURE_PATH, WANDB_DATA


def cd_plots(show_plots: bool = False) -> None:
    """
    Generate the cd sensitivity plots.
    """

    run_info = wandb_summarizer.download.get_results("sdat2/seager19")

    f_df = pd.DataFrame(run_info)[3:13]
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


def nummode_plots(show_plots: bool = False) -> None:
    """
    Make a plot of the number of modes vs. time taken to run the ocean section.
    """
    df = pd.read_csv(WANDB_DATA)

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
