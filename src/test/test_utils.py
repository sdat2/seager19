"""Test `src.utils.py`"""
import numpy as np
from src.utils import timeit, get_byte_size, human_readable_size


def test_timeit() -> None:
    """test `src.utils.timeit` function"""

    @timeit
    def loop(**kwargs) -> None:
        """quick loop function"""
        total = 0
        for _ in range(int(10e2)):
            for _ in range(int(10e2)):
                total += 1
        for kwarg in kwargs:  # to remove pylint error
            print(kwarg)

    tmp_log_d = {}
    loop(log_time=tmp_log_d)
    print(tmp_log_d["loop"])
    assert isinstance(tmp_log_d["loop"], float)
    loop()

    def loop_2(**kwargs) -> None:
        """quick loop function"""
        total = 0
        for _ in range(int(10e2)):
            for _ in range(int(10e2)):
                total += 1
        for kwarg in kwargs:  # to remove pylint error
            print(kwarg)

    timed_loop = timeit(loop_2)
    timed_loop()


def test_get_byte_size() -> None:
    """test `src.utils.get_byte_size` function"""
    print(get_byte_size(np.zeros(int(10e4))))
    print(get_byte_size(list(range(int(10e4)))))
    print(get_byte_size("test string"))
    print(get_byte_size({"great": 1, "test": 0, "ok": -1}))
    # print(get_byte_size(xr.tutorial.load_dataset("air_temperature").air))
    assert isinstance(get_byte_size(np.zeros(int(10e4))), str)


def test_human_readable_size() -> None:
    """test `src.utils.get_byte_size` function"""
    assert human_readable_size(int(10e5)) == "977 KB"
    assert human_readable_size(int(10e13)) == "91 TB"
