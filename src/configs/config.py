"""File for formatting config.py file.

from src.config.configs import format_config

"""
from omegaconf import DictConfig, OmegaConf


def format_config(cfg: DictConfig) -> DictConfig:
    """cfg reformatting.

    Args:
        cfg (DictConfig): cfg to reformat.

    """

    print("OmegaConf.to_yaml(cfg)", OmegaConf.to_yaml(cfg))
    print(cfg.__repr__())

    for i in cfg.atm:
        print("i", i)
        item = cfg.atm[i]
        if isinstance(item, str):
            if "/" in item:
                fl_list = item.split("/")
                print(fl_list)
                if False:  # "$" in fl_list[0]:
                    print(fl_list[0].strip(r"${}"))
                    total = cfg.atm[fl_list[0].strip(r"${}")]
                else:
                    total = float(fl_list[0])
                fl_list.pop(0)
                for j in range(len(fl_list)):
                    if False:  # "$" in fl_list[j]:
                        print(fl_list[j].strip(r"${}"))
                        total = total / cfg.atm[fl_list[j].strip(r"${}")]
                    else:
                        print(float(fl_list[j]))
                        total = total / float(fl_list[j])
                cfg.atm[i] = total

    return cfg
