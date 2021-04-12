"""Test the atmospheric model script.

Example:
    Test using::
        pytest src/test/test_atmos.py

"""
from src.models import atmos


def test_make_figure():
    """test all functions in document"""
    atmos.make_figure()
    atmos.output_trends()
    atmos.output_dq()
