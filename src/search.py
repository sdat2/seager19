"""search.py"""
import os
import numpy as np
from hydra.experimental import initialize, compose
from src.constants import CONFIG_PATH, CONFIG_NAME, SENS_NAME


def rand(low: float, high: float) -> float:
    return float(np.random.uniform(low, high, 1))


if __name__ == "__main__":
    # python src/search.py

    rel_path = CONFIG_PATH.replace(os.getcwd() + "/src/", "")

    with initialize(config_path=rel_path):
        sens = compose(
            config_name=SENS_NAME,
            # overrides=override_list,
        )

    print(sens)
    override_list = list()

    for _ in range(30):

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
