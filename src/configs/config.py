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
            if "/" in item or "*" in item:
                # pylint: disable=eval-used
                cfg.atm[i] = eval(item)
    return cfg
