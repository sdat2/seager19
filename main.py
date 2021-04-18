"""A file to run model runs from with wandb.

Basically a wrapper for bash commands.

Example:
   Usage of script::
       python3 main.py run_name=test22

"""
import os
import wandb
import logging
import hydra
from omegaconf import DictConfig, OmegaConf
from src.constants import (
    SRC_PATH,
    PROJECT_PATH,
)
from src.utils import timeit
from src.models.ocean import compile_all, copy_all, run_all, animate_all


log = logging.getLogger(__name__)


def start_wandb(cfg: DictConfig) -> None:
    """Intialises wandb for run."""
    run_dir = os.path.join(PROJECT_PATH, "logs", cfg.run_name)
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)
    wandb.init(
        project=cfg.project,
        entity=cfg.user,
        dir=run_dir,
        save_code=True,
        name=cfg.run_name,
        notes=cfg.notes,
        config=cfg,
    )


@timeit
@hydra.main(config_path=os.path.join(SRC_PATH, "configs"), config_name="config")
def main(cfg: DictConfig):
    """The main function to run the model and animations.

    Args:
        cfg (DictConfig): The hyrda dict config from the wrapper.

    """

    for i in cfg.atm:
        # print(i)
        # print(type(i))
        item = cfg.atm[i]
        if isinstance(item, str):
            if "/" in item:
                # print(item.split("/"))
                fl_list = item.split("/")
                total = float(fl_list[0])
                fl_list.pop(0)
                for j in range(len(fl_list)):
                    total = total / float(fl_list[j])
                cfg.atm[i] = total

    print("OmegaConf.to_yaml(cfg)", OmegaConf.to_yaml(cfg))
    print(cfg.__repr__())
    # print(type(cfg.__repr__()))

    if not cfg.test:
        start_wandb(cfg)
        compile_all()
        if cfg.run:
            run_all(cfg)
        copy_all()
        if cfg.animate:
            animate_all()


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()
    # 2. Save model inputs and hyperparameters
    # wandb.config
    # config.learning_rate = 0.01
    # 3. Log metrics over time to visualize performance
    # wandb.log({"loss": loss})
    # wandb.log(cfd)
