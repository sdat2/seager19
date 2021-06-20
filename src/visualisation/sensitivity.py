"""Find the sensitivity of the coupled omodel."""
import pandas as pd
import matplotlib.pyplot as plt
import wandb_summarizer.download
from src.models.poly import plot
from src.constants import SEL_DICT, FIGURE_PATH


def cd_plots() -> None:
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
            fig_path=str(FIGURE_PATH / str("cd_" + reg + ".png")),
            ax_format="x",
        )
        plt.clf()
