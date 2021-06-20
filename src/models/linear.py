"""To apply polynomial fits, and propogate error."""
from typing import Callable, Tuple, Sequence, Union, Optional, Literal
import numpy as np
from uncertainties import unumpy as unp, ufloat
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from src.plot_utils import tex_uf, CAM_BLUE, BRICK_RED, OX_BLUE

Flt = Union[float, ufloat]


def _cubic(x: float, a: Flt, b: Flt, c: Flt, d: Flt) -> Flt:
    """Fit cubic curve to data"""
    return a * (x ** 3) + b * (x ** 2) + c * x + d


def _parab(x: float, a: Flt, b: Flt, c: Flt) -> Flt:
    """Fit parabola to data."""
    return a * (x ** 2) + b * x + c


def _lin(x: float, a: Flt, b: Flt) -> Flt:
    """ fit line to data using curve_fit"""
    return (a * x) + b


def _lin_0(x: float, a: Flt) -> Flt:
    """ fit line through zero data using curve_fit"""
    return a * x


def _return_func(param: unp.uarray, reg_type: str = "lin") -> Callable:
    """
    Return function so that the linear function only has to be referenced once.

    Args:
        param (np.ndarray): the param.
        reg_type (str, optional): Which fit occured. Defaults to "lin".

    Returns:
        Callable: Function
    """

    def lin(x: Sequence[Flt]) -> np.array:
        return _lin(np.array(x), param[0], param[1])

    def lin_0(x: Sequence[Flt]) -> np.array:
        return _lin_0(np.array(x), param[0])

    def parab(x: Sequence[Flt]) -> np.array:
        return _parab(np.array(x), param[0], param[1], param[2])

    def cubic(x: Sequence[Flt]) -> np.array:
        return _cubic(np.array(x), param[0], param[1], param[2], param[3])

    if reg_type == "lin":
        return lin
    elif reg_type == "lin0":
        return lin_0
    elif reg_type == "parab":
        return parab
    elif reg_type == "cubic":
        return cubic
    else:
        assert False


def fit(
    x_npa: Sequence[Union[float, int]],
    y_npa: Sequence[Union[float, int]],
    reg_type: Literal["lin_0", "lin", "parab", "cubic"] = "lin",
) -> Tuple[unp.uarray, Callable]:
    """
    Fit a line, with an estimate of the uncertainty.

    Args:
        x_npa (Sequence[Union[float, int]]): The x values to fit.
        y_npa (Sequence[Union[float, int]]): The y values to fit
        reg_type (str, optional): Which regression to do. Defaults to "lin".

    Returns:
        Tuple[unp.uarray, Callable]: Paramaters with uncertainty,
            function to put data into.
    """
    func_dict = {"lin": _lin, "lin0": _lin_0, "parab": _parab, "cubic": _cubic}

    assert reg_type in func_dict

    popt, pcov = curve_fit(func_dict[reg_type], x_npa, y_npa)
    perr = np.sqrt(np.diag(pcov))
    param = unp.uarray(popt, perr)

    return param, _return_func(param, reg_type=reg_type)


def plot(
    x_values: Sequence[Union[float, int]],
    y_values: Sequence[Union[float, int]],
    reg_type: Literal["lin_0", "lin", "parab", "cubic"] = "lin",
    x_label: str = "x label",
    y_label: str = "y label",
    fig_path: Optional[str] = None,
    ax_format: Optional[Literal["both", "x", "y"]] = "both",
) -> Tuple[unp.uarray, Callable]:
    """
    Plot the polynomial.

    Args:
        x_values (Sequence[Union[float, int]]): The x values to fit.
        y_values (Sequence[Union[float, int]]): The y values to fit.
        reg_type (str, optional): Which regression to do. Defaults to "lin".
        x_label (str, optional): X label for plot. e.g.
            r"$\\Delta \\bar{T}_s$ over tropical pacific (pac) region [$\\Delta$ K]"
        y_label (str): Y labelsfor plot. e.g.
            r"$\\Delta \\bar{T}_s$ over nino3.4 region [$\\Delta$ K]"
        fig_path (Optional[str], optional): Path to stor the figure in.
            Defaults to None.
        ax_format (Literal["both", "x", "y"], optional): which axes to format
            in scientific notation. Defaults to "both".

    Returns:
        Tuple[unp.uarray, Callable]: Paramaters with uncertainty,
            function to put data into.

    Example:
        Using plot::

            import wandb_summarizer.download
            from src.models.linear import plot

            run_info = wandb_summarizer.download.get_results("sdat2/seager19")
            f_df = pd.DataFrame(run_info)[3:13]
            f_df = f_df.drop(labels=[11], axis=0)
            nino_3_list = f_df["end_trend_nino3"].tolist()
            nino_4_list = f_df["end_trend_nino4"].tolist()
            plot(
                nino_3_list,
                nino_4_list,
                x_label=r"$\\Delta \\bar{T}_s$ over nino3 region [$\\Delta$ K]",
                y_label=r"$\\Delta \\bar{T}_s$ over nino4 region [$\\Delta$ K]",
                ax_format=None,
                reg_type="lin",
            )
    """

    ext = 0.05
    param, func = fit(x_values, y_values, reg_type=reg_type)
    min_x_data = min(x_values)
    max_x_data = max(x_values)
    min_x_pred = min_x_data - (max_x_data - min_x_data) * ext
    max_x_pred = max_x_data + (max_x_data - min_x_data) * ext

    x_pred = np.linspace(min_x_pred, max_x_pred, num=50)
    y_pred = func(x_pred)
    y_pred_n = unp.nominal_values(y_pred)
    y_pred_s = unp.std_devs(y_pred)

    if len(param) == 1:
        label = "y = " + tex_uf(param[0], bracket=True) + "$x$"
    elif len(param) == 2:
        label = "y = " + tex_uf(param[0], bracket=True) + "$x$  + " + tex_uf(param[1])
    elif len(param) == 3:
        label = (
            "y = "
            + tex_uf(param[0], bracket=True)
            + "$x^2$  + "
            + tex_uf(param[1], bracket=True)
            + "$x$  + "
            + tex_uf(param[2])
        )
    elif len(param) == 4:
        label = (
            "y = "
            + tex_uf(param[0], bracket=True)
            + "$x^3$ + "
            + tex_uf(param[1], bracket=True)
            + "$x^2$ + "
            + tex_uf(param[2], bracket=True)
            + "$x$ + "
            + tex_uf(param[3])
        )
    else:
        print("Too many parameters.")
        assert False

    plt.fill_between(
        x_pred, y_pred_n + y_pred_s, y_pred_n - y_pred_s, alpha=0.5, color=CAM_BLUE
    )
    plt.plot(x_pred, y_pred_n, label=label, color=BRICK_RED, alpha=0.7)
    plt.scatter(x_values, y_values, color=OX_BLUE, alpha=0.7)

    if ax_format is not None:
        plt.gca().ticklabel_format(
            axis=ax_format, style="sci", scilimits=(0, 0), useMathText=True
        )

    if len(param) >= 3:
        plt.legend(
            bbox_to_anchor=(-0.15, 1.02, 1, 0.102),
            loc="lower left",
            mode="expand",
        )
    else:
        plt.legend()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.xlim(min_x_pred, max_x_pred)
    plt.tight_layout()
    if fig_path is not None:
        plt.savefig(fig_path)

    return param, func
