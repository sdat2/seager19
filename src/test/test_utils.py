"""Test `src.utils.py`."""
import xarray as xr
import numpy as np
from src.utils import timeit, get_byte_size, human_readable_size, fix_calendar
from src.data_loading.download import get_data
from src.constants import OCEAN_DATA_PATH


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


def test_fix_calendar() -> None:
    """Test `src.utils.fix_calendar`."""
    get_data()
    ds_new = fix_calendar(
        xr.open_dataset(OCEAN_DATA_PATH / "qflx.nc", decode_times=False)
    )
    da_new = fix_calendar(
        xr.open_dataset(OCEAN_DATA_PATH / "qflx.nc", decode_times=False).qflx
    )
    print(ds_new)
    print(da_new)


def test_get_byte_size() -> None:
    """Test `src.utils.get_byte_size` function."""
    print(get_byte_size(np.zeros(int(10e4))))
    print(get_byte_size(list(range(int(10e4)))))
    print(get_byte_size("test string"))
    print(get_byte_size({"great": 1, "test": 0, "ok": -1}))
    # print(get_byte_size(xr.tutorial.load_dataset("air_temperature").air))
    assert isinstance(get_byte_size(np.zeros(int(10e4))), str)


def test_human_readable_size() -> None:
    """Test `src.utils.get_byte_size` function."""
    assert human_readable_size(int(10e5)) == "977 KB"
    assert human_readable_size(int(10e13)) == "91 TB"
