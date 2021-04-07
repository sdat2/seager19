"""Test `src.visualisation.ani.py`."""
import pandas as pd
import xarray as xr
from src.visualisation.ani import animate_prediction, animate_xr_da


def test_animate_prediction() -> None:
    """Test `src.visualisation.ani.animate_prediction`."""
    da = (
        xr.tutorial.load_dataset("air_temperature")
        .air.isel(time=slice(0, 10))
        .rename({"time": "year", "lat": "y", "lon": "x"})
        - 150
    ) / 150
    x_da_band = xr.concat(
        [da for _ in range(4)], pd.Index(["JFM", "AMJ", "JAS", "OND"], name="mn")
    )
    x_da = xr.concat(
        [x_da_band for _ in range(6)],
        pd.Index(["R", "G", "B", "nir", "SWIR1", "SWIR2"], name="band"),
    )
    y_da = da
    pred_da = da
    animate_prediction(x_da, y_da, pred_da)


def test_animate_xr_da() -> None:
    """Test `src.visualisation.ani.animate_xr_da`."""
    da = (
        xr.tutorial.load_dataset("air_temperature")
        .air.isel(time=slice(0, 10))
        .rename({"time": "time", "lat": "y", "lon": "x"})
        - 150
    ) / 150
    animate_xr_da(da)
    animate_xr_da(da, video_path="output.gif")
