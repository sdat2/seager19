"""
Unit test for downloads.
"""
from src.data_loading.download import (
    get_member_data,
    get_original_models,
    get_figure_data,
    get_noaa_data
)


def test_download() -> None:
    """
    Test download of data.
    """
    get_member_data()
    get_original_models()
    get_figure_data(force_refresh=True)
    get_noaa_data()
