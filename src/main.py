"""A file to run model runs from with hydra/wandb.

Basically a wrapper for bash commands that run the ocean model (fortran/C),
and calls the atmospheric model.

Example:
   Usage of script::

       python3 src/main.py name=test26

"""
import shutil
import wandb
import hydra
from omegaconf import DictConfig
from src.constants import (
    CONFIG_PATH,
    CONFIG_NAME,
    FILES_TO_DELETE,
    run_path,
)
from src.utils import timeit
from src.models.coupling import Coupling
from src.models.model_setup import ModelSetup
from src.configs.config import format_config
from src.wandb_utils import start_wandb
from src.data_loading.download import get_data


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
    Subsection of main to run from a unit test or a sensitivity search.

    Args:
        cfg (DictConfig): The config from whichever method.
        unit_test (bool): Whether or not this is run
            from a unit test. Defaults to False.

    """

    # print("OmegaConf.to_yaml(cfg)", OmegaConf.to_yaml(cfg))
    run_p = run_path(cfg, unit_test=unit_test)

    if cfg.wandb:
        start_wandb(cfg, unit_test=unit_test)

    setup = ModelSetup(run_p, cfg)
    couple = Coupling(cfg, setup)
    couple.run()

    if cfg.wandb:
        wandb.finish()

    if cfg.archive:

        @timeit
        def archive() -> None:
            try:
                shutil.move(
                    run_p,
                    str(cfg.archive_dir),
                )
            # pylint: disable=bare-except
            except:
                print("files not deleted sucessfully.", "run:  rm " + run_p)
                with open(FILES_TO_DELETE, "a") as f:
                    f.write("rm " + run_p)

        archive()


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()
