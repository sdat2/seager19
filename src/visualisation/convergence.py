"""Look at the convergence of the coupling scheme."""
import numpy as np
import matplotlib.pyplot as plt
from src.wandb_utils import metric_conv_data


def metric_conv_plot(
    metric_name: str = "mean_pac", long_name: str = "Mean Tropical Pacific (pac)"
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
    plt.show()

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
    plt.savefig(f"")
    plt.show()
