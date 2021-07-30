"""Comp cd variation for various experiments with varying the inputs to the model."""
import os
import numpy as np
from typing import Optional
import matplotlib.pyplot as plt
from uncertainties import unumpy as unp
from src.constants import FIGURE_PATH
from src.plot_utils import ps_defaults
from src.models.poly import fit
from src.wandb_utils import cd_variation_comp


def plot_comp_cd(
    mem_d: dict, save_path: Optional[str] = None, show_plots: bool = False
) -> None:
    """
    Plot drag coefficient for different inputs.

    Args:
        mem_d (dict): dict with different inputs as keys.
        save_path (Optional[str], optional): Where to save the plot to.
            Defaults to None. If None will not save.
        show_plots (bool, optional): Whether to show plots. Defaults to False.
    """
    ps_defaults(use_tex=False)

    color_d = {
        "EEEE": "blue",
        "EECE": "green",
        "EEEC": "orange",
        "EECC": "red",
    }

    name_d = {
        "EEEE": "ECMWF",
        "EECE": "W",
        "EEEC": "RH",
        "EECC": "RH+W",
    }
    min_x = np.inf
    max_x = -np.inf

    for mem in name_d:
        if min(mem_d[mem][0]) < min_x:
            min_x = min(mem_d[mem][0])
        if max(mem_d[mem][0]) > max_x:
            max_x = max(mem_d[mem][0])

    ext = 0.05
    min_x_pred = min_x - (max_x - min_x) * ext
    max_x_pred = max_x + (max_x - min_x) * ext
    x_pred = np.linspace(min_x_pred, max_x_pred, num=50)

    for mem in name_d:
        param, func = fit(mem_d[mem][0], mem_d[mem][1], reg_type="lin")
        print("param", mem, param)
        y_pred = func(x_pred)
        y_pred_n = unp.nominal_values(y_pred)
        y_pred_s = unp.std_devs(y_pred)
        plt.fill_between(
            x_pred,
            y_pred_n + y_pred_s,
            y_pred_n - y_pred_s,
            alpha=0.25,
            color=color_d[mem],
        )
        plt.plot(x_pred, y_pred_n, color=color_d[mem], alpha=0.5)
        plt.scatter(
            mem_d[mem][0],
            mem_d[mem][1],
            marker="x",
            label=name_d[mem],
            c=color_d[mem],
        )

    plt.xlabel("Drag coefficient, $C_d$, [dimesionless]")
    plt.ylabel("1958-2017 nino3.4 trend [K]")
    plt.legend(
        bbox_to_anchor=(0.0, 1.02, 1, 0.102),
        loc="lower left",
        mode="expand",
        ncol=4,
    )
    plt.gca().ticklabel_format(
        axis="x", style="sci", scilimits=(0, 0), useMathText=True
    )

    plt.xlim(min_x_pred, max_x_pred)
    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path)
    if show_plots:
        plt.show()
    else:
        plt.clf()


if __name__ == "__main__":
    # python src/wandb_utils.py
    plot_comp_cd(
        cd_variation_comp(), save_path=os.path.join(FIGURE_PATH, "mech_sens_0.5.pdf")
    )
    plot_comp_cd(
        cd_variation_comp(e_frac=2),
        save_path=os.path.join(FIGURE_PATH, "mech_sens_2.pdf"),
    )
