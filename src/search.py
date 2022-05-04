"""search.py"""
from typing import List
import os
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
    ARCHIVE_DIR,
)
from src.utils import timeit


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
        choices (List[Char], optional): Characters to choose between.
            Defaults to ["C", "E"].
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
    control: Char = "E",
    exps: List[Char] = ["C", "6"],
    vary: List[bool] = [True, True, True, True],
) -> List[str]:
    """
    Get the full set of options to try if there is one control
    set and multiple experiments deviations.

    Args:
        control (Char, optional): _description_. Defaults to "E".
        exps (List[Char], optional): _description_. Defaults to ["C", "6"].

    Returns:
        List[str]: List of combinations to try.
    """
    length = len([x for x in vary if x])

    def _union(lst1: list, lst2: list) -> list:
        final_list = list(set(lst1) | set(lst2))
        return final_list

    output = []
    for i in exps:
        output = _union(output, between_two(choices=[control, i], length=length))

    for i in range(len(vary)):
        if not vary[i]:
            output = [x[:i] + control + x[i:] for x in output]
    return output


def remainder_combinations() -> List[str]:
    """
    Work out which combinations are still to do.

    Returns:
        List[str]: list to-do.
    """
    big_list = variable_combinations()
    small_list = variable_combinations(vary=[False, False, True, True])
    return [x for x in big_list if x not in small_list]


def var_clt_combinations() -> List[str]:
    """
    Work out which combinations are still to do.

    Returns:
        List[str]: list to-do.
    """
    big_list = variable_combinations(vary=[False, True, True, True], exps=["6"])
    small_list = variable_combinations(vary=[False, False, True, True], exps=["6"])
    return [x for x in big_list if x not in small_list]


def var_ts_combinations() -> List[str]:
    """
    Work out which combinations are still to do.

    Returns:
        List[str]: list to-do.
    """
    big_list = variable_combinations(vary=[True, True, True, True], exps=["6"])
    small_list = variable_combinations(vary=[False, True, True, True], exps=["6"])
    return [x for x in big_list if x not in small_list]


def list_to_hydra_input(comb_list: List[str]) -> str:
    """
    List to hydra.

    Args:
        comb_list (List[str]): List to go through.

    Returns:
        str: string to add to terminal input.
    """
    output = ""
    for i in comb_list:
        output += i
        output += ","
    return output[:-1]


def terminal_call(
    e_frac: str = "0.5,2",
    clouds: str = "true,false",
    mem: str = list_to_hydra_input(remainder_combinations()),
) -> str:
    """
    Return terminal call.

    Args:
        e_frac (str, optional): Defaults to "0.5,2".
        clouds (str, optional): Defaults to "true,false".
        mem (str, optional): Defaults to
            `list_to_hydra_input(remainder_combinations())`.

    Returns:
        str: Terminal call to run model some number of time.
    """
    comp = which_comp(mem)
    return str(
        f"python src/main.py -m atm.e_frac={e_frac} "
        + f"atm.vary_cloud_const={clouds} atm.mem={mem} "
        + f"archive_dir={ARCHIVE_DIR} "
        + f"comp.sst={comp} comp.prwnd={comp}"
    )


def which_comp(mem: str) -> str:
    """
    Which figure to compare with.

    Args:
        mem (str): variable string.

    Returns:
        str: Figure string.
    """
    if mem[3] != "E":
        if mem[2] != "E":
            return "5b"
        else:
            return "5c"
    else:
        return "5a"


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    # python src/search.py
    # main()
    """
    print(list_to_hydra_input(remainder_combinations()))
    print(remainder_combinations())
    print(len(remainder_combinations()))
    print(terminal_call())
    # for comb in remainder_combinations():
    #    print(terminal_call(mem=comb))
    for comb in variable_combinations(
        vary=[False, False, True, True],
    ):
        print(terminal_call(mem=comb))
    """

    for comb in var_ts_combinations():
        print(terminal_call(mem=comb))
