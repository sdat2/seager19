"""PSL data reader.

You can find all of these indices indexed here:

<https://psl.noaa.gov/data/climateindices/list/>"""
import os
import datetime
import urllib
import numpy as np
import xarray as xr
from src.constants import PSL_INDICES_PATH


# All the indices which I think are related to ENSO.
url_d = {
    "nino1+2": "https://psl.noaa.gov/data/correlation/nina1.anom.data",
    "nino3": "https://psl.noaa.gov/data/correlation/nina3.anom.data",
    "nino34": "https://psl.noaa.gov/data/correlation/nina34.anom.data",
    "nino4": "https://psl.noaa.gov/data/correlation/nina4.anom.data",
    "oni": "https://psl.noaa.gov/data/correlation/oni.data",
    "tni": "https://psl.noaa.gov/data/correlation/tni.data",
    "soi": "https://psl.noaa.gov/data/correlation/soi.data",
    "meiv2": "https://psl.noaa.gov/enso/mei/data/meiv2.data",
    # "mei": "https://psl.noaa.gov/enso/mei.old/table.html",
    "best": "https://psl.noaa.gov/data/correlation/censo.data",
    "hurr": "https://psl.noaa.gov/data/correlation/hurr.data",
    "prcp": "https://psl.noaa.gov/data/correlation/espi.data",
    "pacwarm": "https://psl.noaa.gov/data/correlation/pacwarm.data",
}


def index_da(var: str = "tni") -> xr.DataArray:
    """
    Read an index from the psl website,
    and reformat it into an xarray.Datarray object.

    Args:
        var (str, optional): Which index to read. Defaults to "tni".

    Returns:
        xr.DataArray: A standard xarray 1D datarray with the index values,
            and a datetime axis.
            Assumes that the value corresponds to the 15th of the month,
            and nans out missing values.
    """
    url = url_d[var]
    file = urllib.request.urlopen(url)

    def parse_line(decoded_line):
        return [x for x in decoded_line.strip("\n").split(" ") if x != ""]

    stage = 0
    index_list = []
    date_list = []
    description = ""
    for line in file:
        decoded_line = line.decode("utf-8")
        list_line = parse_line(decoded_line)
        if stage == 0:
            first_year = int(list_line[0])
            last_year = int(list_line[1])
            stage = 1
        elif stage == 1:
            year = int(list_line.pop(0))
            date_line = [datetime.datetime(year, i, 15, 0, 0, 0) for i in range(1, 13)]
            index_list.extend(list_line)
            date_list.extend(date_line)
            if year == last_year:
                stage = 2
        elif stage == 2:
            missing_value = float(list_line[0])
            stage = 3
        else:
            description += decoded_line

    index_np = np.asarray(index_list, dtype="float32")
    date_np = np.asarray(date_list)
    var_xr = xr.DataArray(
        data=index_np,
        dims=["time"],
        coords=dict(time=(["time"], date_np)),
        attrs=dict(
            name=var,
            long_name=var.upper(),
            units="dimensionless",
            description=description,
            missing_value=missing_value,
        ),
    ).rename(var)
    return var_xr.where(var_xr != missing_value)


def get_psl_indices(reload: bool = False) -> xr.Dataset:
    """
    A function to return the psl indices, either from memory,
    or by reading them from the web.

    Args:
        reload (bool, optional): Whether to prefer to reload them.
             Defaults to False.

    Returns:
        xr.Dataset: Dataset with indices from the PSL over
            all available datapoints. Missing values are marked as nan.
    """
    if os.path.exists(PSL_INDICES_PATH) and not reload:
        ds = xr.open_dataset(PSL_INDICES_PATH)
    else:
        ds = xr.merge([index_da(var) for var in url_d], join="outer")
        ds.to_netcdf(PSL_INDICES_PATH)

    return ds
