"""Test the plot settings."""
import os
import numpy as np
import matplotlib.pyplot as plt
from src.plot_settings import (
    ps_defaults,
    label_subplots,
    set_dim,
    STD_CLR_LIST,
    BRICK_RED,
    OX_BLUE,
)
from src.constants import PROJECT_PATH


def test_plot() -> None:
    """Function to make example plot."""
    for use_tex in [False, True]:
        ps_defaults(use_tex=use_tex)

        fig, axs = plt.subplots(2, 2)
        x = np.linspace(0, np.pi, num=100)
        axs[0, 0].plot(x, np.sin(x), color=STD_CLR_LIST[0])
        axs[0, 1].plot(x, np.cos(x), color=STD_CLR_LIST[1])
        axs[1, 0].plot(x, np.sinc(x), color=BRICK_RED)
        axs[1, 1].plot(x, np.abs(x), color=OX_BLUE)

        # set size
        set_dim(fig, fraction_of_line_width=1, ratio=(5 ** 0.5 - 1) / 2)
        # label subplots
        label_subplots(axs, start_from=0, fontsize=10)

        plt.savefig(str(os.path.join(PROJECT_PATH, "gifs", "example.png")), dpi=800)
