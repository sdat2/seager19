"""Test that f2py is working.

pytest src/test/test_f2py.py
"""
import os
from src.constants import SRC_PATH


def test_f2py():
    """test2py."""
    os.system(
        "conda activate ./env\n"
        + "cd "
        + os.path.join(SRC_PATH, "test")
        + "\n"
        + "f2py -c primes.f95 -m primes"
    )
    # os.system("cd "+ os.path.join(SRC_PATH, "test") +"\n" +
    #           "python -m numpy.f2py -c primes.f95 -m primes")
    # pylint: disable=import-outside-toplevel
    from src.test import primes

    print(primes.__doc__)
    print(primes.logical_to_integer.__doc__)
    sieve_array = primes.sieve(100)
    prime_numbers = primes.logical_to_integer(sieve_array, sum(sieve_array))
    print(prime_numbers)
