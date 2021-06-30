"""search.py"""
import os
from hydra.experimental import initialize, compose
from src.constants import CONFIG_PATH, CONFIG_NAME, SENS_NAME


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

    for i in sens:
        override_list.append(i + "={:.3e}".format((sens[i][0] + sens[i][1]) / 2))
        print(override_list)

    with initialize(config_path=rel_path):
        cfg = compose(
            config_name=CONFIG_NAME,
            overrides=override_list,
        )

    print(cfg)
    print(cfg.name)
