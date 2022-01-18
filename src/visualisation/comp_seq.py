"""Comp seq function."""
import numpy as np
import numpy.ma as ma
import scipy
import matplotlib.pyplot as plt


def comp_seq(
    name: str,
    varx: np.array,
    vary: np.array,
    varx_name: str = "psl",
    vary_name: str = "our code",
    show_plot: bool = True,
) -> None:
    """
    Function to compare two different sequences to see how closely they match

    Args:
        name (str): The name of the sequence you're comparing.
        varx (np.array): The name of the target variable.
        vary (np.array): The name of the output variable.
        varx_name (str, optional): Name of target variable. Defaults to "psl".
        vary_name (str, optional): Name of output variable. Defaults to "our code".
        show_plot (bool, optional): Whether to show the plot. Defaults to True.
    """
    print("\n\n=====", name, "=====")
    mask = ~np.isnan(varx) & ~np.isnan(vary)
    print(
        "correlation matrix",
        ma.corrcoef(ma.masked_invalid(varx), ma.masked_invalid(vary)),
    )
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(
        varx[mask], vary[mask]
    )
    print(f"slope: {slope}")
    print(f"intercept: {intercept}")
    print(f"r_value: {r_value}")
    print(f"p_value: {p_value}")
    print(f"std_err: {std_err}")
    print("==================\n\n\n")
    plt.scatter(varx, vary)
    plt.xlabel(f"{name} target ({varx_name})")
    plt.ylabel(f"{name} result ({vary_name})")
    if show_plot:
        plt.show()
    plt.clf()
