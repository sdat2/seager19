"""Test the atmospheric model script.

Example:
    Test using::
        pytest src/test/test_atmos.py

"""
from omegaconf import OmegaConf
from src.models.atmos import Atmos
from src.data_loading.download import get_data


def test_make_figure():
    """test all functions in document"""
    get_data()
    atmos = Atmos(OmegaConf.create({"atm": "v", "list": [1, {"a": "1", "b": "2"}]}))
    atmos.make_figure()
    atmos.output_trends()
    atmos.output_dq()
