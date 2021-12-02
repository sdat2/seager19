"""A module to contain my pairplot function."""
import numpy.ma as ma
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns


def pairplot(df_tmp: pd.DataFrame) -> None:
    """
    Improved seaborn pairplot from:

    https://stackoverflow.com/questions/50832204/show-correlation-values-in-pairplot-using-seaborn-in-python/50835066

    Args:
        df_tmp (pd.DataFrame): A data frame.

    Example:
        Usage with psl indices::

            import matplotlib.pyplot as plt
            from src.visualisation.pairplot import pairplot
            from src.data_loading.psl import get_psl_indices

            ds = get_psl_indices()
            df = ds.to_dataframe()
            df1 = df[["tni", "nino34", "meiv2", "soi"]]
            pairplot(df1)
            plt.show()
    """

    def corrfunc(x, y, ax=None, **kws):
        """Plot the correlation coefficient in the top left hand
            corner of a plot.

        A function to use with seaborn's `map_lower` api."""
        corr = ma.corrcoef(ma.masked_invalid(x), ma.masked_invalid(y))
        corr_coeff = corr[0, 1]
        ax = ax or plt.gca()
        ax.annotate(f"ρ = {corr_coeff:.2f}", xy=(0.1, 0.9), xycoords=ax.transAxes)

    g = sns.pairplot(df_tmp, corner=True)
    g.map_lower(corrfunc)
    plt.show()
