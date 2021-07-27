"""Sets up the weights and biases script and
provides functionality to get data from wandb."""
import os
from typing import Optional, List, Union, Tuple
import math
import numpy as np
import pandas as pd
import wandb
import logging
from omegaconf import DictConfig, OmegaConf
import matplotlib.pyplot as plt
from subprocess import PIPE, run
from uncertainties import unumpy as unp
from src.constants import DATA_PATH, run_path, FIGURE_PATH
from src.plot_utils import ps_defaults
from src.models.model_setup import ModelSetup
from src.models.poly import fit


log = logging.getLogger(__name__)


def get_v(inp: str) -> str:
    """Get version of compilers.

    Args:
        inp (str): input string. E.g. "gfortran -v".

    Returns:
        str: selects the fortran version (or gcc version)
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

        wandb.init(
            project=cfg.project,
            entity=cfg.user,
            dir=run_path(cfg, unit_test=unit_test),
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

    for rn in runs:

        # run.summary are the output key/values like accuracy.
        # We call ._json_dict to omit large files
        # pylint: disable=protected-access
        summary_list.append(rn.summary._json_dict)

        # run.config is the input metrics.
        # We remove special values that start with _.
        config = {k: v for k, v in rn.config.items() if not k.startswith("_")}
        config_list.append(config)

        # run.name is the name of the run.
        name_list.append(rn.name)

    summary_df = pd.DataFrame.from_records(summary_list)
    config_df = pd.DataFrame.from_records(config_list)
    name_df = pd.DataFrame({"name": name_list})
    all_df = pd.concat([name_df, config_df, summary_df], axis=1)

    all_df.to_csv(DATA_PATH / "project.csv")

    return all_df


def get_wandb_data(save_path: Optional[str] = None) -> pd.DataFrame:
    """
    Get wandb data (and save it?) now doesn't redownload.

    Args:
        save_path (Optional[str], optional): Path to new csv file. Defaults to None.
            If it is None then doesn't try to save.

    Returns:
        pd.DataFrame: The pandas dataframe of final results.
    """
    # run_info = wandb_summarizer.download.get_results("sdat2/seager19")
    # df = pd.DataFrame(run_info)
    # if save_path is not None:
    #    df.to_csv(save_path)

    df = pd.read_csv(save_path)

    return df


def finished_names() -> List[str]:
    """
    Return all the finished run names.

    Returns:
        List[str]: list of run names.
    """
    api = wandb.Api()
    # Project is specified by <entity/project-name>
    runs = api.runs("sdat2/seager19")
    name_list = [rn.name for rn in runs if rn.state == "finished"]
    return name_list


# pylint: disable=dangerous-default-value
def metric_conv_data(
    metric_name: str = "mean_pac",
    prefix: str = "cd_",
    ex_list: List[str] = [
        "cd_norm",
        "nummode",
    ],
    control_variable_list=[(("atm", "k_days"), 10), (("atm", "e_frac"), 2)],
    index_by: tuple = ("coup", "c_d"),
) -> Tuple[dict, dict]:
    """
    Generate the data for the convergence of a particular item.

    Used in src.visualisation.convergence.metric_conv_plot

    Args:
        metric_name (str, optional): Which keyword to use. Defaults to "mean_pac".

    Returns:
        dict: A dictionary with the relevant metrics in.
    """
    api = wandb.Api()
    # Project is specified by <entity/project-name>
    runs = api.runs("sdat2/seager19")

    def check_controls(config: DictConfig) -> bool:
        truth_list = []
        for i in control_variable_list:
            if i[0][1] in config[i[0][0]]:
                # print(i)
                # print([i[0][0], i[0][1]])
                # pylint: disable=eval-used
                af = config[i[0][0]][i[0][1]]
                # print(af, i[1], af == i[1])
                truth_list.append(af == i[1])
            else:
                truth_list.append(False)
        return np.all(truth_list)

    metric_dict = {}
    setup_dict = {}

    # pylint: disable=unnecessary-comprehension
    for rn in runs:  # [x for x in runs][0:13]:
        if prefix in rn.name and np.all([x not in rn.name for x in ex_list]):
            config = {k: v for k, v in rn.config.items() if not k.startswith("_")}
            config["name"] = rn.name
            cfg = fix_config(config)
            if check_controls(cfg):
                print(cfg.name)
                # print(rn, config["name"])
                pair_list = []
                for _, row in rn.history().iterrows():
                    step = row["_step"]
                    metric = row[metric_name]
                    if not math.isnan(metric):
                        # print(step, metric, type(metric))
                        pair_list.append([step, metric])

                # enforce needing 6 iterations.
                if len(pair_list) == 6:
                    metric_dict[cfg[index_by[0]][index_by[1]]] = np.array(pair_list)
                    setup_dict[cfg[index_by[0]][index_by[1]]] = setup_from_config(cfg)

    return metric_dict, setup_dict


def fix_config(config: Union[dict, DictConfig]) -> DictConfig:
    """
    Turn the config dict back into a DictConfig object.

    Args:
        config (Union[dict, DictConfig]): config dictionary.

    Returns:
        DictConfig:
    """
    if isinstance(config, dict):
        for i in config:
            if "{" in config[i] or "[" in config[i]:
                # pylint: disable=eval-used
                config[i] = eval(config[i])
        return OmegaConf.create(config)
    else:
        return config


def archive_dir_from_config(cfg: Union[DictConfig, dict]) -> str:
    """
    Get the archived folder from the names stored online.

    Args:
        cfg (Union[DictConfig, dict]): The cfg

    Returns:
        str: [description]
    """
    if "archive_dir" not in cfg:
        cfg["archive_dir"] = "unknown"
        print("no archive dir recorded for", cfg.name)
    if isinstance(cfg, DictConfig):
        return os.path.join(cfg.archive_dir, cfg.name)
    elif isinstance(cfg, dict):
        return os.path.join(cfg["archive_dir"], cfg["name"])


def setup_from_config(cfg: Union[DictConfig]) -> ModelSetup:
    """
    Gets the setup object for the archived run from the config.

    Args:
        cfg (Union[DictConfig, dict]): Either the dictconfig or the dict.

    Returns:
        ModelSetup: The model setup object.
    """
    return ModelSetup(archive_dir_from_config(cfg), cfg, make_move=False)


def setup_from_name(name: str) -> ModelSetup:
    """Get the model setup from a name."""
    api = wandb.Api()
    # Project is specified by <entity/project-name>
    runs = api.runs("sdat2/seager19")
    for rn in runs:  # [x for x in runs][0:13]:
        if name == rn.name:
            config = {k: v for k, v in rn.config.items() if not k.startswith("_")}
            config["name"] = rn.name
            cfg = fix_config(config)
            setup = setup_from_config(cfg)
    return setup


def _other_tests():
    """Private function to store test."""
    metric_d, _ = metric_conv_data(
        metric_name="trend_nino3.4",
        prefix="k_days_",
        control_variable_list=[
            (("atm", "eps_days"), 0.75),
            (("atm", "e_frac"), 2),
            (("coup", "c_d"), 2.25e-3),
        ],
        index_by=("atm", "k_days"),
    )
    metric_d, _ = metric_conv_data(
        metric_name="trend_nino3.4",
        prefix="cd_",
        control_variable_list=[
            (("atm", "k_days"), 10),
            (("atm", "e_frac"), 2),
            (("atm", "eps_days"), 0.75),
        ],
        index_by=("atm", "e_frac"),
    )
    metric_d, _ = metric_conv_data(
        metric_name="trend_nino3.4",
        prefix="k_days_",
        control_variable_list=[
            (("atm", "k_days"), 10),
            (("coup", "c_d"), 2.25e-3),
            (("atm", "eps_days"), 0.75),
        ],
        index_by=("atm", "e_frac"),
    )
    print(metric_d)


def cd_variation_comp(e_frac=0.5) -> dict:
    """
    Vary drag coefficient and get the final metric.
    """
    mem_dict = {}
    for mem in ["EEEE", "EECE", "EEEC", "EECC"]:
        print(mem)
        metric_d, _ = metric_conv_data(
            metric_name="trend_nino3.4",
            prefix="k_days_10_eps_days",
            ex_list=["ingrid_True"],
            control_variable_list=[
                (("atm", "k_days"), 10),
                # (("coup", "c_d"), 2.25e-3),
                (("atm", "e_frac"), e_frac),
                (("atm", "eps_days"), 0.75),
                (("atm", "mem"), mem),
                (("atm", "vary_cloud_const"), True),
            ],
            index_by=("coup", "c_d"),
        )
        # pylint: disable=condider-using-dict-items
        mem_dict[mem] = [[], []]
        for val in metric_d:
            print(val, "\t", float(metric_d[val][5, 1]))
            mem_dict[mem][0].append(float(val))
            mem_dict[mem][1].append(float(metric_d[val][5, 1]))

    print(mem_dict)
    return mem_dict


def plot_comp_cd(
    mem_d: dict, save_path: Optional[str] = None, show_plots: bool = False
) -> None:
    """
    Plot drag coefficient for different inputs.

    Args:
        mem_d (dict): dict with different inputs as keys.
        save_path (Optional[str], optional): Where to save the plot to.
            Defaults to None. If None will not save.
        show_plots (bool, optional): Whether to show plots. Defaults to False.
    """
    ps_defaults(use_tex=False)

    color_d = {
        "EEEE": "blue",
        "EECE": "green",
        "EEEC": "orange",
        "EECC": "red",
    }

    name_d = {
        "EEEE": "ECMWF",
        "EECE": "W",
        "EEEC": "RH",
        "EECC": "RH+W",
    }
    min_x = np.inf
    max_x = -np.inf

    for mem in name_d:
        if min(mem_d[mem][0]) < min_x:
            min_x = min(mem_d[mem][0])
        if max(mem_d[mem][0]) > max_x:
            max_x = max(mem_d[mem][0])

    ext = 0.05
    min_x_pred = min_x - (max_x - min_x) * ext
    max_x_pred = max_x + (max_x - min_x) * ext
    x_pred = np.linspace(min_x_pred, max_x_pred, num=50)

    for mem in name_d:
        param, func = fit(mem_d[mem][0], mem_d[mem][1], reg_type="lin")
        print("param", mem, param)
        y_pred = func(x_pred)
        y_pred_n = unp.nominal_values(y_pred)
        y_pred_s = unp.std_devs(y_pred)
        plt.fill_between(
            x_pred,
            y_pred_n + y_pred_s,
            y_pred_n - y_pred_s,
            alpha=0.25,
            color=color_d[mem],
        )
        plt.plot(x_pred, y_pred_n, color=color_d[mem], alpha=0.5)
        plt.scatter(
            mem_d[mem][0],
            mem_d[mem][1],
            marker="x",
            label=name_d[mem],
            c=color_d[mem],
        )

    plt.xlabel("Drag coefficient, $C_d$, [dimesionless]")
    plt.ylabel("1958-2017 nino3.4 trend [K]")
    plt.legend(
        bbox_to_anchor=(0.0, 1.02, 1, 0.102),
        loc="lower left",
        mode="expand",
        ncol=4,
    )
    plt.gca().ticklabel_format(
        axis="x", style="sci", scilimits=(0, 0), useMathText=True
    )

    plt.xlim(min_x_pred, max_x_pred)
    plt.tight_layout()

    if save_path is not None:
        plt.savefig(save_path)
    if show_plots:
        plt.show()
    else:
        plt.clf()


if __name__ == "__main__":
    # python src/wandb_utils.py
    _other_tests()
    plot_comp_cd(
        cd_variation_comp(), save_path=os.path.join(FIGURE_PATH, "mech_sens_0.5.pdf")
    )
    plot_comp_cd(
        cd_variation_comp(e_frac=2),
        save_path=os.path.join(FIGURE_PATH, "mech_sens_2.pdf"),
    )
