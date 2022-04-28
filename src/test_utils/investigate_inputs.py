""""""
from typing import Union
import pathlib
from src.xr_utils import open_dataarray, sel, spatial_mean


def da_diagnostics(path: Union[pathlib.Path, str], reg: str = "nino3.4") -> None:
    """
    Check a datarray by printing it's values over the region.

    Args:
        path (Union[pathlib.Path, str]): Path to xarray datarray netcdf.
        reg (str, optional): Region. Defaults to "nino3.4".
    """
    da = open_dataarray(path)
    da_avg = spatial_mean(sel(da, reg=reg))
    print(da_avg)


if __name__ == "__main__":
    sst_path = "/home/users/sithom/seager19/atmos/DATA/sst-ECMWF-clim.nc"
    da_diagnostics(sst_path, reg="nino3.4")
    sst_path = "/home/users/sithom/seager19/atmos/DATA/sst-CMIP6-clim.nc"
    da_diagnostics(sst_path, reg="nino3.4")
    taux_path = "/home/users/sithom/seager19/ocean/DATA/tau-ECMWF-clim.x"
    da_diagnostics(taux_path, reg="nino3.4")
    tauy_path = "/home/users/sithom/seager19/ocean/DATA/tau-CMIP6-clim.y"
    da_diagnostics(tauy_path, reg="nino3.4")
