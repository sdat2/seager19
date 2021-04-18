"""File for formatting config.py file.

from src.config.configs import format_config

"""
from omegaconf import DictConfig, OmegaConf


def format_config(cfg: DictConfig):
    """cfg reformatting.

    Args:
        cfg (DictConfig): cfg to reformat.

    """

    print("OmegaConf.to_yaml(cfg)", OmegaConf.to_yaml(cfg))
    print(cfg.__repr__())

    for i in cfg.atm:
        item = cfg.atm[i]
        if isinstance(item, str):
            if "/" in item:
                fl_list = item.split("/")
                total = float(fl_list[0])
                fl_list.pop(0)
                for j in range(len(fl_list)):
                    total = total / float(fl_list[j])
                cfg.atm[i] = total

    return cfg
