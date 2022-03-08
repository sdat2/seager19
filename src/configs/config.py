"""File for formatting config.py file.


Example:
    Usage at top of program::

        from src.config.configs import format_config

    Instantiation somewhere in the progam::

        cfg = format_config(cfg)

"""
from dataclasses import dataclass
from omegaconf import DictConfig
import hydra
from hydra.core.config_store import ConfigStore


@dataclass
class Atm:
    mem: str = "EEEE"


@dataclass
class Config:
    atm: Atm = Atm()


cs = ConfigStore.instance()
cs.store(name="config", node=Config)

# database_lib registers its configs
# in database_lib/db
# database_lib.register_configs()


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


@hydra.main(config_path=None, config_name="config")
def test(cfg: Config) -> None:
    print(cfg)


if __name__ == "__main__":
    # python src/configs/config.py
    # pylint: disable=no-value-for-parameter
    test()
