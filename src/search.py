"""search.py"""
from typing import List
import os
import numpy as np
from itertools import product
import hydra
from hydra.experimental import initialize, compose
from omegaconf import DictConfig
from pyparsing import Char
from src.constants import (
    CONFIG_PATH,
    CONFIG_NAME,
    SENS_RANGES,
    SENS_SETTINGS,
    PROJECT_PATH,
)
from src.utils import timeit


def rand(low: float, high: float) -> float:
    """
    Get a single random float from the range.

    A wrapper round `np.random.uniform`, so that I can quickly change this
    if need be from here. Could replace this with a latin hypercube search,
    or equivalently Sobel indices etc.
    These alternatives that they would more fairly sample the space.

    Args:
        low (float): low end.
        high (float): high end.

    Returns:
        float: random float.
    """
    return float(np.random.uniform(low, high, 1))


@timeit
@hydra.main(config_path=CONFIG_PATH, config_name=SENS_SETTINGS)
def main(settings: DictConfig) -> None:
    """The main function to run the model and animations.

    Takes the src/configs/config.yaml file as input alongside any command line
    arguments.

    Args:
        settings (DictConfig): The hyrda dict config from the wrapper.

    """
    print(CONFIG_PATH)
    print(PROJECT_PATH)
    print(os.getcwd())
    # rel_path = CONFIG_PATH.replace(os.getcwd() + "/src/", "")
    rel_path = CONFIG_PATH.replace(str(PROJECT_PATH), "../../..")
    print(rel_path)

    print(settings)

    with initialize(config_path=rel_path):
        sens = compose(
            config_name=SENS_RANGES,
            # overrides=override_list,
        )

    print(sens)
    override_list = list()

    for _ in range(settings.runs):

        for i in sens:
            override_list.append(i + "={:.3e}".format(rand(sens[i][0], sens[i][1])))

            # print(override_list)

        with initialize(config_path=rel_path):
            cfg = compose(
                config_name=CONFIG_NAME,
                overrides=override_list,
            )

        # print(cfg)
        print(cfg.name)


def between_two(choices: List[Char] = ["C", "E"], length: int = 4) -> List[str]:
    """
    All possible string sequences betweeen two characters for some length.

    Args:
        choices (List[Char], optional): Characters to choose between. Defaults to ["C", "E"].
        length (int, optional): _description_. Defaults to 4.

    Returns:
        (List[str]): list of possible sequences.
    """
    assert len(choices) == 2
    output = []
    for x in map(
        "".join,
        product(
            *zip(
                [choices[0] for _ in range(length)],
                [choices[1] for _ in range(length)],
            )
        ),
    ):
        output.append(x)
    return output


def variable_combinations(
    control: Char = "E", exps: List[Char] = ["C", "6"]
) -> List[str]:
    """
    Get the full set of options to try if there is one control
    set and multiple experiments deviations.

    Args:
        control (Char, optional): _description_. Defaults to "E".
        exps (List[Char], optional): _description_. Defaults to ["C", "6"].

    Returns:
        List[str]: _description_
    """

    def _union(lst1: list, lst2: list) -> list:
        final_list = list(set(lst1) | set(lst2))
        return final_list

    output = []
    for i in exps:
        output = _union(output, between_two(choices=[control, i]))
    return output


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    # python src/search.py
    # main()
    print(variable_combinations())
