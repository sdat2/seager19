"""To apply linear regression functions."""
from typing import Callable, Tuple
import numpy as np
from uncertainties import unumpy as unp
from scipy.optimize import curve_fit


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
        Callable: [description]
    """

    def lin(x):
        return _lin(np.array(x), popt[0], popt[1])

    def lin_0(x):
        return _lin_0(np.array(x), popt[0])

    if reg_type == "lin":
        return lin
    elif reg_type == "lin0":
        return lin_0
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
    func_dict = {"lin": _lin, "lin0": _lin_0}

    assert reg_type in func_dict

    popt, pcov = curve_fit(_lin, x_npa, y_npa)
    perr = np.sqrt(np.diag(pcov))

    return unp.uarray(popt, perr), _return_func(popt, reg_type=reg_type)
