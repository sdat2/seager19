"""test `src.visualisation.ani.py`"""
import pandas as pd
import xarray as xr
from src.visualisation.ani import animate_prediction


def test_animate_prediction() -> None:
    """test `src.visualisation.ani.animate_prediction` """
    print("testing animate prediction")
    da = (
        xr.tutorial.load_dataset("air_temperature")
        .air.isel(time=slice(0, 7))
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
