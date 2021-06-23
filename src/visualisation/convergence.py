"""Look at the convergence of the coupling scheme."""
import numpy as np
import matplotlib.pyplot as plt
from src.plot_utils import ps_defaults
from src.wandb_utils import metric_conv_data
from src.constants import FIGURE_PATH

ps_defaults(use_tex=False, dpi=200)


def metric_conv_plot(
    metric_name: str = "mean_pac",
    long_name: str = "Mean Tropical Pacific (pac)",
    show_plots: bool = False,
):
    """
    Make the convergence plot for a particular metric.

    Args:
        metric_name (str, optional): The keyword to extract. Defaults to "mean_pac".
        long_name (str, optional): The long name for the ylabel.
            Defaults to "Mean Tropical Pacific (pac)".
    """

    metric_dict = metric_conv_data(metric_name=metric_name)

    for cd in sorted(metric_dict):
        plt.plot(
            metric_dict[cd][:, 0], metric_dict[cd][:, 1], label="{:.2e}".format(cd)
        )

    plt.legend(
        bbox_to_anchor=(-0.15, 1.02, 1.15, 0.102),
        loc="lower left",
        mode="expand",
        ncol=5,
    )
    plt.xlabel("Step")
    plt.ylabel(long_name + r" [$^{\circ}$C]")
    plt.savefig(FIGURE_PATH / str(metric_name + "_convergence.png"))
    plt.savefig(FIGURE_PATH / str(metric_name + "_convergence.pdf"))

    if show_plots:
        plt.show()
    else:
        plt.clf()

    for cd in metric_dict:
        plt.plot(
            np.abs(
                (metric_dict[cd][:, 1] - metric_dict[cd][5, 1]) / metric_dict[cd][5, 1]
            ),
            label="{:.2e}".format(cd),
        )
    plt.ylabel(long_name + r" $ \frac{|C-F|}{F}$")
    plt.yscale("log")
    plt.xlabel("Step")
    plt.legend(
        bbox_to_anchor=(-0.15, 1.02, 1.15, 0.102),
        loc="lower left",
        mode="expand",
        ncol=5,
    )
    plt.savefig(FIGURE_PATH / str(metric_name + "_convergence_log.png"))
    plt.savefig(FIGURE_PATH / str(metric_name + "_convergence_log.pdf"))

    if show_plots:
        plt.show()
    else:
        plt.clf()


if __name__ == "__main__":
    # python src/visualisation/convergence.py
    metric_conv_plot()
    metric_conv_plot(metric_name="trend_nino3.4", long_name="Trend nino3.4")
