"""To apply linear regression functions."""
from typing import Callable, Tuple, Sequence, Union
import numpy as np
from uncertainties import unumpy as unp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def _parab(x: float, a: float, b: float, c: float):
    """Fit parabola to data."""
    return a * (x ** 2) + b * x + c


def _lin(x: float, a: float, b: float) -> float:
    """ fit line to data using curve_fit"""
    return (a * x) + b


def _lin_0(x: float, a: float) -> float:
    """ fit line through zero data using curve_fit"""
    return a * x


def _return_func(popt: np.ndarray, reg_type: str = "lin") -> Callable:
    """
    Return function so that the linear function only has to be referenced once.

    Args:
        popt (np.ndarray): the param.
        reg_type (str, optional): Which fit occured. Defaults to "lin".

    Returns:
        Callable: Function
    """

    def lin(x):
        return _lin(np.array(x), popt[0], popt[1])

    def lin_0(x):
        return _lin_0(np.array(x), popt[0])

    def parab(x):
        return _parab(np.array(x), popt[0], popt[1], popt[2])

    if reg_type == "lin":
        return lin
    elif reg_type == "lin0":
        return lin_0
    elif reg_type == "parab":
        return parab
    else:
        assert False


def fit(
    x_npa: np.ndarray, y_npa: np.ndarray, reg_type: str = "lin"
) -> Tuple[unp.uarray, Callable]:
    """
    Fit a line, with an estimate of the uncertainty.

    Args:
        x_npa (np.ndarray): The x values to fit.
        y_npa (np.ndarray): The y values to fit
        reg_type (str, optional): [description]. Defaults to "lin".

    Returns:
        Tuple[unp.uarray, Callable]: Paramaters with uncertainty,
            function to put data into.
    """
    func_dict = {"lin": _lin, "lin0": _lin_0, "parab": _parab}

    assert reg_type in func_dict

    popt, pcov = curve_fit(func_dict[reg_type], x_npa, y_npa)
    perr = np.sqrt(np.diag(pcov))

    return unp.uarray(popt, perr), _return_func(popt, reg_type=reg_type)


def plot(
    x_values: Sequence[Union[float, int]],
    y_values: Sequence[Union[float, int]],
    x_label: str,
    y_label: str,
    fig_path: str,
) -> Tuple[unp.uarray, Callable]:
    """
    [summary]

    [extended_summary]

    Args:
        x_values (Sequence[Union[float, int]]): The x values to fit.
        y_values (Sequence[Union[float, int]]): The y values to fit.
        x_label (str): X label for plot. e.g.
            r"$\\Delta \\bar{T}_s$ over tropical pacific (pac) region [$\\Delta$ K]"
        y_label (str): Y labelsfor plot. e.g.
            r"$\\Delta \\bar{T}_s$ over nino3.4 region [$\\Delta$ K]"
        fig_path (str): Path to stor the figure in.

    Returns:
        Tuple[unp.uarray, Callable]: Paramaters with uncertainty,
            function to put data into.
    """
    plt.scatter(x_values, y_values)
    ext = 0.05

    param, func = fit(x_values, y_values)
    min_x_data = min(x_values)
    max_x_data = max(x_values)
    min_x_pred = min_x_data - (max_x_data - min_x_data) * ext
    max_x_pred = max_x_data + (max_x_data - min_x_data) * ext

    x_pred = np.linspace(min_x_pred, max_x_pred, num=50)
    y_pred = func(x_pred)
    plt.plot(
        x_pred,
        y_pred,
        label="y = ({:2.1f}".format(param[0]) + ")$x$  +  {:2.1f}".format(param[1]),
        color="red",
    )
    plt.legend()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.tight_layout()
    plt.savefig(fig_path)

    return param, func
