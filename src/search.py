"""search.py"""
import os
import numpy as np
import hydra
from hydra.experimental import initialize, compose
from omegaconf import DictConfig
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


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()
