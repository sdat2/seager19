"""Test the atmospheric model script.

Example:
    Test using::
        pytest src/test/test_atmos.py

"""
from src.models import atmos
from src.data_loading.download import get_data


def test_make_figure():
    """test all functions in document"""
    get_data()
    atmos.make_figure()
    atmos.output_trends()
    atmos.output_dq()
