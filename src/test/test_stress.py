"""
Test stress movement.
"""
import numpy as np
import xarray as xr
from src.constants import TEST_DIREC


def test_inversibility() -> None:
    """Test inversibility of the stress files using NETCDF3_CLASSIC."""
    tau_clim_obj = xr.open_dataarray(
        str(TEST_DIREC / "ocean" / "DATA" / "tau-ECMWF-clim.x"), decode_times=False
    )
    tau_clim_obj.to_netcdf(str(TEST_DIREC / "test.x"), format="NETCDF3_CLASSIC")
    byte_num = 100

    with open(str(TEST_DIREC / "ocean" / "DATA" / "tau-ECMWF-clim.x"), "r+b") as f:
        byt1 = f.read(byte_num)
        print(byt1)

    with open("test.x", "r+b") as f:
        byt2 = f.read(byte_num)
        print(byt2)

    np.all([byt1[x] == byt2[x] for x in range(len(byt1))])
