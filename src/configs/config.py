"""File for formatting config.py file.


Example:
    Usage at top of program::

        from src.config.configs import format_config

    Instantiation somewhere in the progam::

        cfg = format_config(cfg)

"""
from omegaconf import DictConfig


def format_config(cfg: DictConfig) -> DictConfig:
    """cfg reformatting.

    Currently just evaluates the arithemtic put into the program.
    This isn't the best idea in the world, but it works, and
    it's broadly your own sorry fault if you put in
    dangerous code to the command line and it gets evaluated.

    Args:
        cfg (DictConfig): cfg to reformat.

    """

    # print("OmegaConf.to_yaml(cfg)", OmegaConf.to_yaml(cfg))
    # print(cfg.__repr__())

    for i in cfg.atm:
        item = cfg.atm[i]
        if isinstance(item, str):
            if "/" in item or "*" in item:
                # pylint: disable=eval-used
                cfg.atm[i] = eval(item)
    return cfg
