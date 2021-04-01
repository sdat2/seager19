"""General project util functions"""
from typing import Callable
import inspect
import time
from functools import wraps
from sys import getsizeof


def timeit(method: Callable) -> Callable:
    """timeit is a wrapper for performance analysis which should
    return the time taken for a function to run. Alters `log_time` `dict` if fed in.
    Add @timeit to the function you want to time. Function needs `**kwargs` if
    you want it to be able to feed in `log_time` `dict`.

    Args:
        method (Callable):  the function that it takes as an input

    Examples:
        Here is an example with the tracking functionality and without::

            >>> @timeit
            ... def loop(**kwargs):
            ...     total = 0
            ...     for i in range(int(10e2)):
            ...         for j in range(int(10e2)):
            ...             total += 1
            >>> tmp_log_d = {}
            >>> loop(log_time=tmp_log_d)
            >>> print(tmp_log_d["loop"])
            >>> loop()
    """

    @wraps(method)
    def timed(*args, **kw):
        ts = time.perf_counter()
        result = method(*args, **kw)
        te = time.perf_counter()
        if "log_time" in kw:
            name = kw.get("log_name", method.__name__.lower())
            kw["log_time"][name] = te - ts
        else:
            print("%r  %2.5f s\n" % (method.__name__, (te - ts)))
        return result

    return timed


def human_readable_size(num: int, suffix: str = "B") -> str:
    """
    Convert a number of bytes into human readable format.
    This function is meant as a helper function for `get_byte_size`.
    Args:
        num (int): The number of bytes to convert
        suffix (str, optional): The suffix to use for bytes. Defaults to 'B'.
    Returns:
        str: A human readable version of the number of bytes.
    """
    assert num >= 0, "Size cannot be negative."
    for unit in ["", "K", "M", "G", "T", "P", "E", "Z"]:
        if num < 1024:
            return f"{num:.0f} {unit}{suffix}"
        num /= 1024.0
    return f"{num:1f} Y{suffix}"


def calculate_byte_size_recursively(obj: object, seen: set = None) -> int:
    """
    Recursively calculate size of objects in memory in bytes.
    From: https://github.com/bosswissam/pysize. Meant as a helper function for
    `get_byte_size`.
    Args:
        obj (object): The python object to get the size of
        seen (set, optional): This variable is needed to for the recusrive
            function evaluations, to ensure each object only gets counted once.
            Leave it at "None" to get the full byte size of an object. Defaults to None.
    Returns:
        int: The size of the object in bytes
    """

    # Note: getsizeof alone is not enough, as it only returns the size of the top
    #  level object, not of its member variables/objects. Hence the recursive calls.
    size = getsizeof(obj)
    if seen is None:
        # At first iteration (top level object), initialize 'seen' as empty set
        seen = set()

    obj_id = id(obj)
    if obj_id in seen:
        # If object was already counted, return 0 size to avoid double counting.
        return 0

    # Important: Mark as seen *before* entering recursion to handle
    # self-referential objects
    seen.add(obj_id)

    if hasattr(obj, "__dict__"):
        # handles class objects
        for cls in obj.__class__.__mro__:
            if "__dict__" in cls.__dict__:
                d = cls.__dict__["__dict__"]
                if inspect.isgetsetdescriptor(d) or inspect.ismemberdescriptor(d):
                    # Recursively calculate size of member objects & variables
                    size += calculate_byte_size_recursively(obj.__dict__, seen)
                break
    if isinstance(obj, dict):
        # handles dictionaries
        size += sum((calculate_byte_size_recursively(v, seen) for v in obj.values()))
        size += sum((calculate_byte_size_recursively(k, seen) for k in obj.keys()))
    elif hasattr(obj, "__iter__") and not isinstance(obj, (str, bytes, bytearray)):
        # handles array like objects (need to exclude str, bytes bytearray since they
        #  also implement __iter__)
        size += sum((calculate_byte_size_recursively(i, seen) for i in obj))

    if hasattr(obj, "__slots__"):  # can have __slots__ with __dict__
        size += sum(
            calculate_byte_size_recursively(getattr(obj, s), seen)
            for s in obj.__slots__
            if hasattr(obj, s)
        )
    return size


def get_byte_size(obj: object) -> str:
    """
    Return human readable size of a python object in bytes.
    Args:
        obj (object): The python object to analyse
    Returns:
        str: Human readable string with the size of the object
    """

    return human_readable_size(calculate_byte_size_recursively(obj))
