"""Program to do."""
from typing import Tuple
import numpy as np
from src.utils import timeit


@timeit
def default() -> Tuple[np.array]:
    """
    Default function.

    Returns:
        Tuple[np.array]: Two random arrays.
    """
    return np.array([0, 0, 0]), np.array([0.5, 0.5, 0.5])
