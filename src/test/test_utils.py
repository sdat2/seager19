"""Test utils.py"""
import numpy as np
from src.utils import timeit, get_byte_size


def test_timeit() -> None:
    """test `src.utils.timeit` function"""

    @timeit
    def loop(**kwargs) -> None:
        """quick loop function"""
        total = 0
        for _ in range(int(10e2)):
            for _ in range(int(10e2)):
                total += 1

    tmp_log_d = {}
    loop(log_time=tmp_log_d)
    print(tmp_log_d["loop"])
    assert isinstance(tmp_log_d["loop"], float)
    loop()


def test_get_byte_size() -> None:
    """test `src.utils.get_byte_size` function"""
    print(get_byte_size(np.zeros(int(10e4))))
    print(get_byte_size(list(range(int(10e4)))))
    print(get_byte_size("test string"))
    assert isinstance(get_byte_size(np.zeros(int(10e4))), str)
