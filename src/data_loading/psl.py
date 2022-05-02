"""PSL data reader.

You can find all of these indices indexed here:

<https://psl.noaa.gov/data/climateindices/list/>"""
import os
import numpy.ma as ma
import datetime
import shutil
import urllib
import numpy as np
import xarray as xr
from src.utils import timeit, time_stamp
from src.constants import PSL_INDICES_PATH, ERSSTV5_PATH
from src.xr_utils import can_coords, fix_calendar
from src.metrics import nino_calculate

# All the indices which I think are related to ENSO.
# perhaps this should be used to `src.constants`.
# and also replaced with dict comprehension.
url_d = {
    "nino1+2": "https://psl.noaa.gov/data/correlation/nina1.anom.data",
    "nino3": "https://psl.noaa.gov/data/correlation/nina3.anom.data",
    "nino3.4": "https://psl.noaa.gov/data/correlation/nina34.anom.data",
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


def psl_index_da(var: str = "tni") -> xr.DataArray:
    """
    Read an index from the psl website,
    and reformat it into an xarray.Datarray object.

    Args:
        var (str, optional): Which index to read. Defaults to "tni".

    Returns:
        xr.DataArray: A standard xarray 1D datarray with the index values,
            and a datetime axis.
            Assumes that the value corresponds to the 1st of the month,
            and nans out missing values.
    """

    def parse_line(decoded_line):
        return [x for x in decoded_line.strip("\n").split(" ") if x != ""]

    url = url_d[var]
    with urllib.request.urlopen(url) as file:

        stage = 0
        index_list = []
        date_list = []
        description = ""
        for line in file:
            decoded_line = line.decode("utf-8")
            list_line = parse_line(decoded_line)
            if stage == 0:
                last_year = int(list_line[1])
                stage = 1
            elif stage == 1:
                year = int(list_line.pop(0))
                date_line = [
                    datetime.datetime(year, i, 1, 0, 0, 0) for i in range(1, 13)
                ]
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
            units="dimensionless",  # TODO: This is wrong for most metrics; change.
            description=description,
            missing_value=missing_value,
        ),
    ).rename(var)
    return var_xr.where(var_xr != missing_value)


@timeit
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
        ds = xr.merge([psl_index_da(var) for var in url_d], join="outer")
        ds.attrs = {
            "name": "PSL data climate indices",
            "units": "Variable",
            "Description": "PSL data scraped from the website on " + time_stamp(),
        }
        ds.to_netcdf(PSL_INDICES_PATH)

    return ds


@timeit
def get_ersstv5(reload: bool = False) -> xr.DataArray:
    """
    Get ERSSTV5 datarray.

    Args:
        reload (bool, optional): Whether to prefer redownloading. Defaults to False.

    Returns:
        xr.DataArray: Straight from website.
    """
    if os.path.exists(ERSSTV5_PATH) and not reload:
        da = xr.open_dataset(ERSSTV5_PATH).sst
    else:
        # SHOULD this url not be in constants
        url = "https://downloads.psl.noaa.gov/Datasets/noaa.ersst.v5/sst.mnmean.nc"
        name = "sst.mnmean.nc"
        urllib.request.urlretrieve(url, name)
        shutil.move(name, ERSSTV5_PATH)
        da = xr.open_dataset(ERSSTV5_PATH).sst

    return da


@timeit
def psl_metric_test() -> None:
    """Test to see if the psl and I agree on the metrics given the
    ERSSTv5 SST data (i.e. to debug metrics etc.)"""
    ds = get_psl_indices()
    ersstv5 = get_ersstv5()
    ersstv5_slim = ersstv5.sel(time=slice("1948", "2021"))
    ersstv5_can = fix_calendar(can_coords(ersstv5_slim)).isel(variable=0)
    for x, y in [("nino34", "nino3.4"), ("nino4", "nino4"), ("nino1+2", "nino1+2")]:
        m, _ = nino_calculate(ersstv5_can, reg=y)
        my_nino = m.values
        psl_nino = ds[x].sel(time=slice("1948", "2021-10")).values
        print(ma.corrcoef(ma.masked_invalid(my_nino), ma.masked_invalid(psl_nino)))


if __name__ == "__main__":
    # python src/data_loading/psl.py
    # get_ersstv5(reload=True)
    print(get_psl_indices(reload=True))
