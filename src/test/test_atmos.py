"""Test the atmospheric model script.

Example:
    Test using::
        pytest src/test/test_atmos.py

"""
from src.models.atmos import Atmos
from src.data_loading.download import get_data
from src.configs.load_config import load_config
from src.models.model_setup import ModelSetup
from src.constants import TEST_DIREC
from src.utils import delete_folder_contents


def test_atmos():
    """test all functions in document."""
    get_data()
    delete_folder_contents(str(TEST_DIREC))
    cfg = load_config()
    setup = ModelSetup(str(TEST_DIREC))
    atmos = Atmos(cfg, setup)
    atmos.run_all()
    delete_folder_contents(str(TEST_DIREC))
    # atmos.make_figure()
    # atmos.output_trends()
    # atmos.output_dq()
