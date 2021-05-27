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
