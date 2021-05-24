"""Test the ingrid replacement python functions."""
import os
import xarray as xr
import recursive_diff
from hydra.experimental import initialize, compose
from src.data_loading.ingrid import linear_qflx_replacement
from src.constants import OCEAN_DATA_PATH
from src.constants import PROJECT_PATH, CONFIG_NAME, CONFIG_PATH
from src.data_loading.download import get_data
from src.models.model_setup import ModelSetup
from src.configs.config import format_config


def test_ingrid() -> None:
    """Test the qflx replacement function."""

    # get_data if it does not exists
    get_data()

    with initialize(
        config_path="../../" + str(CONFIG_PATH).replace(str(PROJECT_PATH) + "/", "")
    ):
        # config is relative to a module
        cfg = compose(
            config_name=CONFIG_NAME, overrides=["user=test_user", "name=test_run"]
        )
        print(cfg)
        cfg = format_config(cfg)
        print(cfg)
    setup = ModelSetup(cfg)

    # make a qflx-test file.
    linear_qflx_replacement(setup, output_file_name="qflx-test.nc")

    print("files in paths:\t\t\n", str(os.listdir(str(OCEAN_DATA_PATH))))

    # load different qflx files.
    qflx_test = xr.open_dataarray(
        str(OCEAN_DATA_PATH / "qflx-test.nc"), decode_times=False
    )
    qflx_old = xr.open_dataarray(str(OCEAN_DATA_PATH / "qflx-0.nc"), decode_times=False)

    # look at stuff
    print(qflx_test)
    print(qflx_old)

    # xr.testing.assert_equal(qflx_test, qflx_old)
    xr.testing.assert_allclose(qflx_test, qflx_old, atol=1e-2)

    for x in recursive_diff.recursive_diff(qflx_test, qflx_old, abs_tol=1e-2):
        print(x)

        if "-999.0" not in str(x):

            # otherwise it would catch errors about what
            # the missing value attribute is.
            assert False
