"""A file to run model runs from with wandb.

Basically a wrapper for bash commands.

Example:
   Usage of script::
       python3 main.py name=test26

"""
import wandb
import logging
import hydra
from omegaconf import DictConfig
from src.constants import CONFIG_PATH, CONFIG_NAME
from src.utils import timeit
from src.models.ocean import Ocean
from src.models.atmos import Atmos
from src.models.model_setup import ModelSetup
from src.configs.config import format_config
from src.configs.wandb_util import start_wandb

log = logging.getLogger(__name__)


@timeit
@hydra.main(config_path=CONFIG_PATH, config_name=CONFIG_NAME)
def main(cfg: DictConfig) -> None:
    """The main function to run the model and animations.

    Args:
        cfg (DictConfig): The hyrda dict config from the wrapper.

    """

    cfg = format_config(cfg)

    # print("OmegaConf.to_yaml(cfg)", OmegaConf.to_yaml(cfg))

    start_wandb(cfg)

    setup = ModelSetup(str(wandb.run.dir))

    # ocean model
    ocean = Ocean(cfg, setup)
    ocean.compile_all()
    if cfg.run:
        ocean.run_all()
    if cfg.animate:
        ocean.animate_all()

    # atmos model.
    if cfg.atmos:
        # atmos takes in cfg
        atmos = Atmos(cfg, setup)
        atmos.run_all()

    wandb.finish()


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()
    # 2. Save model inputs and hyperparameters
    # wandb.config
    # config.learning_rate = 0.01
    # 3. Log metrics over time to visualize performance
    # wandb.log({"loss": loss})
    # wandb.log(cfd)
