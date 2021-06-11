"""Test the stress script."""
import numpy as np
import pandas as pd
import xarray as xr
from src.models.coupling import Coupling, ModelSetup
from src.configs.load_config import load_config


def test_stress() -> None:
    """Stress test."""

    time = pd.date_range("2014-09-06", periods=3)
    reference_time = pd.Timestamp("2014-09-05")

    u_vel = 15 + 8 * np.random.randn(2, 2, 3)
    v_vel = 10 * np.random.rand(2, 2, 3)
    lon = [[-99.83, -99.32], [-99.79, -99.23]]
    lat = [[42.25, 42.21], [42.63, 42.59]]

    u_da = xr.DataArray(
        data=u_vel,
        dims=["x", "y", "time"],
        coords=dict(
            lon=(["x", "y"], lon),
            lat=(["x", "y"], lat),
            time=time,
            reference_time=reference_time,
        ),
        attrs=dict(),
    )

    v_da = xr.DataArray(
        data=v_vel,
        dims=["x", "y", "time"],
        coords=dict(
            lon=(["x", "y"], lon),
            lat=(["x", "y"], lat),
            time=time,
            reference_time=reference_time,
        ),
        attrs=dict(),
    )
    cfg = load_config()
    setup = ModelSetup("")

    couple = Coupling(
        cfg,
        setup,
    )

    tau_u, tau_v = couple.f_stress(0.5, u_da, v_da)

    assert isinstance(tau_u, xr.DataArray)
    assert isinstance(tau_v, xr.DataArray)

    couple.replace_stress(1)
    couple.replace_dq(1)
