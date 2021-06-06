"""Sets up the weights and biases script."""
import os
import wandb
import logging
from omegaconf import DictConfig
from subprocess import PIPE, run
from src.constants import PROJECT_PATH, TEST_DIREC


log = logging.getLogger(__name__)


def get_v(inp: str) -> str:
    """Get version of compilers.

    Args:
        inp (str): input string. E.g. "gfortran -v".

    Returns:
        str: selects the fortran version
            part of the output.
    """
    # pylint: disable=subprocess-run-check
    return run(
        inp, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True
    ).stderr.split("\n")[-2]


def start_wandb(cfg: DictConfig, unit_test: bool = False) -> None:
    """
    Intialises wandb for run.

    Weights and biases provides the run tracking for the model runs at different
    parameter settings.

    TODO: Need to improve the ability to initialise wandb from a unit test.

    Args:
        cfg (DictConfig): The config settings to pass in to the wandb syncing.
        unit_test (bool, optional): Whether or not this is a unit-test.
            Defaults to False. If this is a unit test will currently not initialise
            wandb, but will call related functions.

    """
    if not unit_test:
        run_dir = os.path.join(PROJECT_PATH, "logs", cfg.name)
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
    else:
        run_dir = str(TEST_DIREC)

    if not unit_test:

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

        wandb.config.update({"gfortran": get_v("gfortran -v"), "gcc": get_v("gcc -v")})
    else:
        print({"gfortran": get_v("gfortran -v"), "gcc": get_v("gcc -v")})
