"""transforms.py - Script to store transformations to datasets files."""


def rdict(index: int) -> dict:
    """Returns renaming dict for xarray.DataArray output of ocean model.

    Made to reformat the output datarrays of the Fortran
    ocean model used.

    Args:
        index (int): index on coords.

    Returns:
        dict: renaming dict.

    """
    return {
        "T_0" + str(index): "time",
        "Y_0" + str(index): "y",
        "X_0" + str(index): "x",
        "L_0" + str(index): "Z",
    }


def qflx_rdict() -> dict:
    """
    To plot the heat flux it is useful to transform to the standard dict.

    Returns:
        dict: The rename dictionary to apply to the xr.DataArray.
    """
    return {
        "T": "time",
        "Y": "y",
        "X": "x",
        "Z": "Z",
    }
