"""A file to run model runs from with wandb.

Basically a wrapper for bash commands.

Example:
   Usage of script::

       python3 src/main.py name=test26

"""
import wandb
import logging
import hydra
from omegaconf import DictConfig
from src.constants import CONFIG_PATH, CONFIG_NAME, run_dir
from src.utils import timeit
from src.models.coupling import Coupling
from src.models.model_setup import ModelSetup
from src.configs.config import format_config
from src.wandb_utils import start_wandb
from src.data_loading.download import get_data

log = logging.getLogger(__name__)


@timeit
@hydra.main(config_path=CONFIG_PATH, config_name=CONFIG_NAME)
def main(cfg: DictConfig) -> None:
    """The main function to run the model and animations.

    Takes the src/configs/config.yaml file as input alongside any command line
    arguments.

    Args:
        cfg (DictConfig): The hyrda dict config from the wrapper.

    """
    get_data()  # try to get the data at this point if not already downloaded.
    cfg = format_config(cfg)
    sub_main(cfg)


def sub_main(cfg: DictConfig, unit_test: bool = False) -> None:
    """
    Subsection of main to run from a unit test.

    Args:
        cfg (DictConfig): The config from whichever method.
        unit_test (bool): Whether or not this is run
            from a unit test. Defaults to False.

    """

    # print("OmegaConf.to_yaml(cfg)", OmegaConf.to_yaml(cfg))

    if cfg.wandb:
        start_wandb(cfg, unit_test=unit_test)

    setup = ModelSetup(run_dir(), cfg)

    couple = Coupling(cfg, setup)
    couple.run()

    if cfg.wandb:
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
