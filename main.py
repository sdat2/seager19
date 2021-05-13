"""A file to run model runs from with wandb.

Basically a wrapper for bash commands.

Example:
   Usage of script::
       python3 main.py run_name=test26

"""
import os
import wandb
import logging
import hydra
from omegaconf import DictConfig, OmegaConf
from src.constants import PROJECT_PATH, CONFIG_PATH, CONFIG_NAME
from src.utils import timeit
from src.models.ocean import compile_all, copy_all, run_all, animate_all
from src.models.atmos import output_trends, output_dq, make_figure
from src.configs.config import format_config


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
        # pylint: disable=protected-access
        config=cfg._content,
    )


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
    output_trends(wandb.run.dir)
    output_dq(wandb.run.dir)
    make_figure(wandb.run.dir)


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()
    # 2. Save model inputs and hyperparameters
    # wandb.config
    # config.learning_rate = 0.01
    # 3. Log metrics over time to visualize performance
    # wandb.log({"loss": loss})
    # wandb.log(cfd)
