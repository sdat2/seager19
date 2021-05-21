"""A file to run model runs from with wandb.

Basically a wrapper for bash commands.

Example:
   Usage of script::
       python3 main.py name=test26

"""
import os
import wandb
import logging
import hydra
from omegaconf import DictConfig, OmegaConf
from subprocess import PIPE, run
from src.constants import PROJECT_PATH, CONFIG_PATH, CONFIG_NAME
from src.utils import timeit
from src.models.ocean import compile_all, copy_all, run_all, animate_all
from src.models.atmos import Atmos
from src.configs.config import format_config


log = logging.getLogger(__name__)


def start_wandb(cfg: DictConfig, unit_test=False) -> None:
    """
    Intialises wandb for run.

    Args:
        cfg (DictConfig): [description]
        unit_test (bool, optional): Whether or not this is a unit-test.
            Defaults to False.

    """
    if not unit_test:
        run_dir = os.path.join(PROJECT_PATH, "logs", cfg.name)
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
    else:
        run_dir = ""
    wandb.init(
        project=cfg.project,
        entity=cfg.user,
        dir=run_dir,
        save_code=True,
        name=cfg.name,
        notes=cfg.notes,
        # pylint: disable=protected-access
        config=cfg._content,
    )

    def get_v(inp: str) -> None:
        """Get version of compilers"""
        # pylint: disable=subprocess-run-check
        return run(
            inp, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True
        ).stderr.split("\n")[-2]

    wandb.config.update({"gfortran": get_v("gfortran -v"), "gcc": get_v("gcc -v")})


@timeit
@hydra.main(config_path=CONFIG_PATH, config_name=CONFIG_NAME)
def main(cfg: DictConfig) -> None:
    """The main function to run the model and animations.

    Args:
        cfg (DictConfig): The hyrda dict config from the wrapper.

    """

    cfg = format_config(cfg)

    print("OmegaConf.to_yaml(cfg)", OmegaConf.to_yaml(cfg))

    start_wandb(cfg)

    # ocean model
    compile_all()
    if cfg.run:
        run_all(cfg)
    copy_all(cfg)
    if cfg.animate:
        animate_all(cfg)

    # atmos model.
    if cfg.atmos:
        atmos = Atmos(cfg)
        atmos.run_all(direc=wandb.run.dir)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()
    # 2. Save model inputs and hyperparameters
    # wandb.config
    # config.learning_rate = 0.01
    # 3. Log metrics over time to visualize performance
    # wandb.log({"loss": loss})
    # wandb.log(cfd)
