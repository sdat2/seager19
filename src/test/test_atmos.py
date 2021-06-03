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


def test_atmos():
    """test all functions in document."""
    get_data()
    cfg = load_config()
    setup = ModelSetup(str(TEST_DIREC))
    atmos = Atmos(cfg, setup)
    atmos.run_all()
    # atmos.make_figure()
    # atmos.output_trends()
    # atmos.output_dq()
