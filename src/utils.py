"""General project utility functions."""
from typing import Callable
import inspect
import time
from functools import wraps
from sys import getsizeof
import signal
from contextlib import contextmanager
from src.models.model_setup import ModelSetup
from src.configs.load_config import load_config
from src.constants import K_LOGS


class TimeoutException(Exception):
    pass


@contextmanager
def time_limit(seconds: int) -> None:
    """Time limit manager.

    Function taken from:

    https://stackoverflow.com/questions/366682/
    how-to-limit-execution-time-of-a-function-call

    Args:
        seconds (int): how  many seconds to wait until timeout.

    Example:
        Call a function which will take longer than the time limit::

            import time
            from src.utils import time_limit, TimeoutException

            def long_function_call():
                for t in range(5):
                    print("t=", t, "seconds")
                    time.sleep(1)
            try:
                with time_limit(3):
                    long_function_call()
                    assert False
            except TimeoutException as e:
                print("Timed out!")
            except:
                print("A different exception")

    """

    def _signal_handler(signum, frame):
        raise TimeoutException("Timed out!")

    signal.signal(signal.SIGALRM, _signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


def time_stamp() -> str:
    """
    Return the current local time.

    Returns:
        str: Time string format "%Y-%m-%d %H:%M:%S".
    """
    current_time = time.localtime()
    return time.strftime("%Y-%m-%d %H:%M:%S", current_time)


def hr_time(time_in: float) -> str:
    """
    Return human readable time as string.

    I got fed up with converting the number in my head.
    Probably runs very quickly.

    Args:
        time (float): time in seconds

    Returns:
        str: string to print.

    Example:
        120 seconds to human readable string::

            >>> from src.utils import hr_time
            >>> hr_time(120)
                "2 min 0 s"
    """
    if time_in < 60:
        return "%2.5f s" % time_in
    elif 60 < time_in < 60 * 60:
        return time.strftime("%M min %S s", time.gmtime(time_in))
    elif 60 * 60 < time_in < 24 * 60 * 60:
        return time.strftime("%H hr %M min %S s", time.gmtime(time_in))
    else:
        return "%2.5f s" % time_in


def timeit(method: Callable) -> Callable:
    """`src.timeit` is a wrapper for performance analysis.

    It should return the time taken for a function to run. Alters `log_time` `dict`
    if fed in. Add @timeit to the function you want to time. Function needs
    `**kwargs` if you want it to be able to feed in `log_time` `dict`.

    Args:
        method (Callable):  the function that it takes as an input

    Examples:
        Here is an example with the tracking functionality and without::

            >>> from src.utils import timeit
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
        # time.gmtime()
        if "log_time" in kw:
            name = kw.get("log_name", method.__name__.lower())
            kw["log_time"][name] = te - ts
            print("%r " % method.__name__, hr_time(te - ts), "\n")
        else:
            print("%r " % method.__name__, hr_time(te - ts), "\n")
        return result

    return timed


def human_readable_size(num: int, suffix: str = "B") -> str:
    """Convert a number of bytes into human readable format.

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
    """Recursively calculate size of objects in memory in bytes.

    From: https://github.com/bosswissam/pysize.
    Meant as a helper function for `get_byte_size`.

    Args:
        obj (object): The python object to get the size of
        seen (set, optional): This variable is needed to for the recusrive
            function evaluations, to ensure each object only gets counted once.
            Leave it at "None" to get the full byte size of an object. Defaults to None.

    Returns:
        int: The size of the object in bytes.

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
    """Return human readable size of a python object in bytes.

    Args:
        obj (object): The python object to analyse

    Returns:
        str: Human readable string with the size of the object

    """

    return human_readable_size(calculate_byte_size_recursively(obj))


def get_default_setup() -> ModelSetup:
    """Return the default run setup to get the data."""
    run_dir = str(K_LOGS / "k_days_10")
    cfg = load_config(test=False)
    setup = ModelSetup(run_dir, cfg, make_move=False)
    return setup


def in_notebook() -> bool:
    """
    Check if in notebook.

    Taken from this answer:
    https://stackoverflow.com/a/22424821

    Returns:
        bool: whether in notebook.
    """
    try:
        # pylint: disable=import-outside-toplevel
        from IPython import get_ipython

        if "IPKernelApp" not in get_ipython().config:  # pragma: no cover
            return False
    except ImportError:
        return False
    except AttributeError:
        return False
    return True
