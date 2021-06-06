"""Test `src.utils.py`.

pytest src/test/test_utils.py
"""
import numpy as np
from src.utils import timeit, get_byte_size, human_readable_size, hr_time


def test_timeit() -> None:
    """Test `src.utils.timeit` function."""

    @timeit
    def loop(**kwargs) -> None:
        """Quick loop function."""
        total = 0
        for _ in range(int(10e2)):
            for _ in range(int(10e2)):
                total += 1
        for kwarg in kwargs:
            # to remove pylint error
            print(kwarg)

    tmp_log_d = {}
    loop(log_time=tmp_log_d)
    print(tmp_log_d["loop"])
    assert isinstance(tmp_log_d["loop"], float)
    loop()

    def loop_2(**kwargs) -> None:
        """Quick loop function."""
        total = 0
        for _ in range(int(10e2)):
            for _ in range(int(10e2)):
                total += 1
        for kwarg in kwargs:  # to remove pylint error
            print(kwarg)

    timed_loop = timeit(loop_2)
    timed_loop()


def test_get_byte_size() -> None:
    """Test `src.utils.get_byte_size` function."""
    print(get_byte_size(np.zeros(int(10e4))))
    print(get_byte_size(list(range(int(10e4)))))
    print(get_byte_size("test string"))
    print(get_byte_size({"great": 1, "test": 0, "ok": -1}))
    assert isinstance(get_byte_size(np.zeros(int(10e4))), str)

    class ExampleClass:
        def __init__(self) -> None:
            self.prop_dict = {"great": 1, "test": 0, "ok": -1}
            self.list = list(range(int(10e4)))

    example = ExampleClass()

    print(get_byte_size(example))


def test_human_readable_size() -> None:
    """Test `src.utils.get_byte_size` function."""
    assert human_readable_size(int(10e5)) == "977 KB"
    assert human_readable_size(int(10e13)) == "91 TB"


def test_hr_time() -> None:
    """Test `src.utils.hr_time`."""
    assert hr_time(2) == "%2.5f s" % 2
    assert hr_time(10e6) == "%2.5f s" % 10e6
