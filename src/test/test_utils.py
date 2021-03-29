"""Test utils.py"""
import numpy as np
from src.utils import timeit, get_byte_size


def test_timeit():
    @timeit
    def loop(**kwargs):
        total = 0
        for i in range(int(10e2)):
            for j in range(int(10e2)):
                total += 1
    tmp_log_d={}
    part = loop(log_time=tmp_log_d)
    print(tmp_log_d["loop"])
    loop()


def test_get_size():
    print(get_byte_size(np.zeros(int(10e4))))
