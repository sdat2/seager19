"""Test the atmospheric model script.

Example:
    Test using::
        pytest src/test/test_atmos.py

"""
from src.models.atmos import Atmos
from src.data_loading.download import get_data
from src.test.test_hydra import test_with_initialize
from src.models.model_setup import ModelSetup


def test_atmos():
    """test all functions in document."""
    get_data()
    cfg = test_with_initialize()
    setup = ModelSetup(cfg)
    atmos = Atmos(cfg, setup)
    atmos.run_all()
    # atmos.make_figure()
    # atmos.output_trends()
    # atmos.output_dq()
