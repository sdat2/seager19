"""Script to cut out spin up conditions.

Seeks to replicate::

    ```
    \\begin{ingrid}

    (DATA/qflx-0.nc)readCDF .X /XM exch def
    (DATA/qflx-0.nc)readCDF .Y /YM exch def

    (output/om_diag.nc)readCDF .SST_QFLX
    T last 11 sub last RANGE
    T /T (months since 1960-01-01) periodic 0.5 1 11.5 NewEvenGRID replaceGRID
    /X X periodic setgridtype def
    Y -91 1 91 evengridAverage
    0 replaceNaN
    L /Z renameGRID
    (qflx)rn
    (DATA/qflx.nc)writeCDF
    \\end{ingrid}
    ```

Example:

    Usage within the `ocean.RUN.run-model.sh` script::
        >>> conda activate ../../env
        >>> python3 ingrid.py

"""
import numpy as np
import xarray as xr
from src.constants import OCEAN_DATA_PATH, OCEAN_OUTPUT_PATH
from src.utils import timeit


@timeit
def linear_qflx_replacement(output_file_name: str = "qflx.nc") -> None:
    """Uses `xarray` linear interpolation to replace netcdf with qflx in."""
    sst_qflx = xr.open_dataset(
        OCEAN_OUTPUT_PATH / "om_diag.nc", decode_times=False
    ).SST_QFLX.rename({"L_01": "Z", "T_01": "T", "X_01": "X", "Y_01": "Y"})
    lent = len(sst_qflx.coords["T"])
    sst_qflx_subset = sst_qflx.isel(T=slice(lent - 12, lent + 1)).rename("qflx")
    sst_qflx_subset.coords["T"] = [x + 0.5 for x in range(12)]
    sst_qflx_subset.attrs = {
        "units": "unitless",
        "file_missing_value": -987654.0,
        "long_name": "qflx",
        "longname": "qflx",
    }
    sst_qflx_subset.coords["T"].attrs = {
        "modulus": 12.0,
        "modulo": 12.0,
        "pointwidth": 1.0,
        "calendar": "360",
        "gridtype": 1,
        "units": "months since 1960-01-01",
    }
    sst_qflx_subset.coords["Z"].attrs = {
        "long_name": "Level",
        "gridtype": 0,
        "units": "level",
    }
    sst_qflx_subset.coords["Y"].attrs = {
        "uniquename": "Y",
        "pointwidth": 1.0,
        "gridtype": 0,
        "units": "degree_north",
    }
    sst_qflx_subset.coords["X"].attrs = {
        "modulus": 360.0,
        "uniquename": "X",
        "pointwidth": 1.0,
        "gridtype": 1,
        "units": "degree_east",
    }
    sst_qflx_subset = sst_qflx_subset.interp(
        Y=np.array(list(range(-90, 91)), dtype="float32"),
        X=np.array(list(range(0, 360)), dtype="float32"),
        kwargs={"fill_value": 0.0},
        method="linear",
    ).fillna(0.0)
    sst_qflx_subset = sst_qflx_subset.astype("float32")
    sst_qflx_subset.to_netcdf(OCEAN_DATA_PATH / output_file_name)


def test() -> None:
    """Test the qflx replacement function."""
    linear_qflx_replacement(output_file_name="qflx-test.nc")
    qflx_test = xr.open_dataarray(OCEAN_DATA_PATH / "qflx-test.nc", decode_times=False)
    qflx_old = xr.open_dataarray(OCEAN_DATA_PATH / "qflx-0.nc", decode_times=False)

    print(qflx_test)

    print(qflx_old)

    # xr.testing.assert_equal(qflx_test, qflx_old)
    xr.testing.assert_allclose(qflx_test, qflx_old, atol=1e-2)


if __name__ == "__main__":
    # main
    linear_qflx_replacement()
