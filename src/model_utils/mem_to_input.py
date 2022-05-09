"""Convert mem sting into input function"""
from typing import List
import pandas as pd
from src.constants import (
    MODEL_NAMES,
    VAR_DICT,
)


def mems_to_df(mem_list: List[str]) -> pd.DataFrame:
    """
    Turn a list of mems into a dataframe of inputs.

    Args:
        mem_list (List[str]): List of mem to turn into corresponding
            dataframe of inputs.

    Returns:
        pd.DataFrame: A dataframe

    Example:
        Work out what inputs a list of runs got::

            mems_to_df(["EEEE", "CCCC", "66E6"])
    """
    results_lol = []
    for _ in VAR_DICT:
        results_lol.append([])
    for _, mem in enumerate(mem_list):
        for j in VAR_DICT:
            if len(mem) <= j:
                # ECMWF inputs are the default.
                results_lol[j].append(MODEL_NAMES["E"])
            else:
                results_lol[j].append(MODEL_NAMES[mem[j]])
    return pd.DataFrame(
        data={VAR_DICT[i]: results_lol[i] for i in range(len(results_lol))},
        index=mem_list,
    )


def mem_to_dict(mem: str) -> dict:
    """
    Single mem variable to dictionary of inputs.

    Uses logic in `mems_to_df`.

    Args:
        mem (str): The mem input e.g "EEEE"

    Returns:
        dict: dictionary of inputs, e.g.
    """
    return {var: input_d[mem] for (var, input_d) in mems_to_df([mem]).to_dict().items()}
