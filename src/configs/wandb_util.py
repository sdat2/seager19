"""Sets up the weights and biases script."""
import os
from typing import Optional
import wandb_summarizer.download
import pandas as pd
import wandb
import logging
from omegaconf import DictConfig
from subprocess import PIPE, run
from src.constants import LOG_PATH, TEST_DIREC, DATA_PATH


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
        run_dir = os.path.join(LOG_PATH, cfg.name)
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


def get_full_csv() -> pd.DataFrame:
    """
    Get the full csv.
    """
    api = wandb.Api()

    # Project is specified by <entity/project-name>
    runs = api.runs("sdat2/seager19")
    summary_list = []
    config_list = []
    name_list = []

    for run in runs:

        # run.summary are the output key/values like accuracy.
        # We call ._json_dict to omit large files
        summary_list.append(run.summary._json_dict)

        # run.config is the input metrics.
        # We remove special values that start with _.
        config = {k: v for k, v in run.config.items() if not k.startswith("_")}
        config_list.append(config)

        # run.name is the name of the run.
        name_list.append(run.name)

    summary_df = pd.DataFrame.from_records(summary_list)
    config_df = pd.DataFrame.from_records(config_list)
    name_df = pd.DataFrame({"name": name_list})
    all_df = pd.concat([name_df, config_df, summary_df], axis=1)

    all_df.to_csv(DATA_PATH / "project.csv")

    return all_df


def get_wandb_data(save_path: Optional[str] = None) -> pd.DataFrame:
    """
    Get wandb data (and save it?).

    Args:
        save_path (Optional[str], optional): Path to new csv file. Defaults to None.
            If it is None then doesn't try to save.

    Returns:
        pd.DataFrame: The pandas dataframe of final results.
    """
    run_info = wandb_summarizer.download.get_results("sdat2/seager19")
    df = pd.DataFrame(run_info)
    if save_path is not None:
        df.to_csv(save_path)

    return df
