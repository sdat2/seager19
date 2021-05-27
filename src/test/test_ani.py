"""Test `src.visualisation.ani.py`.

pytest src/test/test_ani.py
"""
import os
import pandas as pd
import xarray as xr
from src.visualisation.ani import animate_prediction, animate_xr_da, animate_ds
from src.constants import OCEAN_OUTPUT_PATH, GIF_PATH
from src.data_loading.download import get_data


def test_animate_prediction() -> None:
    """Test `src.visualisation.ani.animate_prediction`."""
    da = (
        xr.tutorial.load_dataset("air_temperature")
        .air.isel(time=slice(0, 5))
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
    animate_prediction(x_da, y_da, pred_da, video_path="gifs/tepr_output.mp4")
    animate_prediction(x_da, y_da, pred_da, video_path="gifs/tepr_output.gif")


def test_animate_xr_da() -> None:
    """Test `src.visualisation.ani.animate_xr_da`."""
    da = (
        xr.tutorial.load_dataset("air_temperature")
        .air.isel(time=slice(0, 5))
        .rename({"time": "time", "lat": "y", "lon": "x"})
        - 150
    ) / 150
    animate_xr_da(da, video_path="gifs/test_output.mp4")
    animate_xr_da(da, video_path="gifs/test_output.gif")


def test_animate_ds() -> None:
    """Animate the sst into gifs."""
    get_data()

    for x in [
        "om_diag",
    ]:
        animate_ds(
            xr.open_dataset(
                str(os.path.join(OCEAN_OUTPUT_PATH, x)) + ".nc", decode_times=False
            ),
            x,
            GIF_PATH,
            dpi=200,
            front_trim=1,
        )
