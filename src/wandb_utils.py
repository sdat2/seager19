"""Sets up the weights and biases script and
provides functionality to get data from wandb."""
import os
from typing import Optional, List, Union, Tuple
import math
import numpy as np
import pandas as pd
import wandb
import logging
from uncertainties import ufloat
from omegaconf import DictConfig, OmegaConf
from subprocess import PIPE, run
from src.utils import timeit, in_notebook
from src.constants import DATA_PATH, run_path, DEFAULT_PROJECT
from src.models.model_setup import ModelSetup
from src.model_utils.mem_to_input import mems_to_df

log = logging.getLogger(__name__)


def _get_runs(project: str = DEFAULT_PROJECT) -> wandb.Api.runs:
    """
    Get the wandb runs using the api.

    Args:
        project (str, optional): Wandb run. E.g. "sdat2/seager19".
            Defaults to DEFAULT_PROJECT.

    Returns:
        wandb.Api.runs: _description_
    """
    api = wandb.Api(timeout=20)
    # Project is specified by <entity/project-name>
    return api.runs(project)


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


def get_full_csv(project: str = DEFAULT_PROJECT) -> pd.DataFrame:
    """
    Get the full csv.
    """
    runs = _get_runs(project=project)
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


def finished_names(project: str = DEFAULT_PROJECT) -> List[str]:
    """
    Return all the finished run names.

    Returns:
        List[str]: list of run names.
    """
    runs = _get_runs(project=project)
    name_list = [rn.name for rn in runs if rn.state == "finished"]
    return name_list


# pylint: disable=dangerous-default-value
def metric_conv_data(
    metric_name: str = "mean_pac",
    prefix: str = "cd_",
    ex_list: List[str] = ["cd_norm", "nummode",],
    control_variable_list=[(("atm", "k_days"), 10), (("atm", "e_frac"), 2)],
    index_by: tuple = ("coup", "c_d"),
    project: str = DEFAULT_PROJECT,
) -> Tuple[dict, dict]:
    """
    Generate the data for the convergence of a particular item.

    Used in src.visualisation.convergence.metric_conv_plot

    Args:
        metric_name (str, optional): Which keyword to use. Defaults to "mean_pac".

    Returns:
        Tuple[dict, dict]: metric_dict, setup_dict.
    """
    runs = _get_runs(project=project)

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
                    if didnt_blow_up(rn):
                        if len(index_by) == 2:
                            metric_dict[cfg[index_by[0]][index_by[1]]] = np.array(
                                pair_list
                            )
                            setup_dict[
                                cfg[index_by[0]][index_by[1]]
                            ] = setup_from_config(cfg)
                        else:
                            metric_dict[cfg[index_by]] = np.array(pair_list)
                            setup_dict[cfg[index_by]] = setup_from_config(cfg)

    return metric_dict, setup_dict


def didnt_blow_up(rn: wandb.apis.public.Run) -> bool:
    """
    Test if the run blew up. True if it didn't blow up.

    Args:
        rn (wandb.apis.public.Run): run.

    Returns:
        bool: whether there was any blow up during the run.
    """
    # limits (degrees celsius).
    limits = {"mean_pac": [15, 30], "mean_nino3.4": [20, 30]}
    results = []
    rn_hist = rn.scan_history(keys=list(limits.keys()))
    for region in limits:
        results.append(
            np.all(
                [
                    row[region] > limits[region][0] and row[region] < limits[region][1]
                    for row in rn_hist
                ]
            )
        )
    return np.all(results)


def fix_config(config: Union[dict, DictConfig]) -> DictConfig:
    """
    Turn the config dict back into a DictConfig object.

    Args:
        config (Union[dict, DictConfig]): config dictionary.

    Returns:
        DictConfig: original configuration.
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
        cfg (Union[DictConfig, dict]): The config from the run.

    Returns:
        str: The archive directory path string.
    """
    if "archive_dir" not in cfg:
        cfg["archive_dir"] = "unknown"
        print("no archive dir recorded for", cfg.name)
    if isinstance(cfg, DictConfig):
        return os.path.join(cfg.archive_dir, cfg.name)
    elif isinstance(cfg, dict):
        return os.path.join(cfg["archive_dir"], cfg["name"])


def setup_from_config(cfg: DictConfig) -> ModelSetup:
    """
    Gets the setup object for the archived run from the config.

    Args:
        cfg (DictConfig): Either the dictconfig or the dict.

    Returns:
        ModelSetup: The model setup object.
    """
    return ModelSetup(archive_dir_from_config(cfg), cfg, make_move=False)


def setup_from_name(name: str, project: str = DEFAULT_PROJECT) -> ModelSetup:
    """Get the model setup from a name.

    Args:
        name (str): model name.

    Returns:
        ModelSetup: The model setup object.
    """
    runs = _get_runs(project=project)
    for rn in runs:  # [x for x in runs][0:13]:
        if name == rn.name:
            config = {k: v for k, v in rn.config.items() if not k.startswith("_")}
            config["name"] = rn.name
            cfg = fix_config(config)
            setup = setup_from_config(cfg)
    return setup


def _other_tests() -> None:
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


def cd_variation_comp(e_frac: float = 0.5) -> dict:
    """
    Vary drag coefficient and get the final metric.

    Args:
        e_frac (float): Defaults to 0.5.

    Returns:
        dict: mem_dict.
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


def _get_cfg(rn) -> DictConfig:
    config = {k: v for k, v in rn.config.items() if not k.startswith("_")}
    config["name"] = rn.name
    return fix_config(config)


PARAM = ["c_d", "eps_days", "eps_frac", "vary_cloud_const"]
PARAM_HYDRA = ["coup.c_d", "atm.eps_days", "atm.e_frac", "atm.vary_cloud_const"]

RESULTS = [
    "trend_nino3.4 [degC]",
    "mean_nino3.4 [degC]",
    "mean_pac [degC]",
]


@timeit
def summary_table(project: str = DEFAULT_PROJECT) -> pd.DataFrame:
    """
    Key input parameters, key output parameters, in a simple dataframe.

    index=number

    paramters=mem, ${ts}${clt}${sfcwind}${rh}${pr}${ps}${tau}, c_d, eps_frac, eps,

    Key indexes: trend_nino3.4, mean_nino3.4,  mean_pac

    Args:
        project (str, optional): Which weights and biases project to scan.
            Defaults to DEFAULT_PROJECT.

    Returns:
        pd.DataFrame: A dataframe.
    """
    runs = _get_runs(project=project)
    mem_list = []
    metric_l = PARAM + RESULTS
    metric_d = {key: [] for key in metric_l}

    for rn in runs:
        cfg = _get_cfg(rn)
        summary = rn.summary
        try:
            res = [
                cfg.coup.c_d,
                cfg.atm.eps_days,
                cfg.atm.e_frac,
                cfg.atm.vary_cloud_const,
                summary["trend_nino3.4"],
                summary["mean_nino3.4"],
                summary["mean_pac"],
            ]
            for i, key in enumerate(metric_l):
                metric_d[key].append(res[i])
            mem_list.append(cfg.atm.mem)
        # pylint: disable=broad-except
        except Exception as e:
            if not in_notebook:
                print(e)

    metric_df = pd.DataFrame(metric_d)
    mem_df = mems_to_df(mem_list).reset_index(0)
    combined_df = pd.concat([mem_df, metric_df], axis=1, join="inner")

    return combined_df


@timeit
def _aggregate_matches(
    summary_df: pd.DataFrame, filter_df: pd.DataFrame
) -> List[pd.DataFrame]:
    df_list = []
    for _, row in filter_df.iterrows():
        # print(i, row)]
        df = summary_df
        for column_number, entry in enumerate(row):
            column = filter_df.columns[column_number]
            # print(i, filter_df.columns[column_number], entry)
            df = df[df[column] == entry]
        df_list.append(df)
    return df_list


def find_missing(df_list: List[pd.DataFrame], param: List[str] = PARAM) -> None:
    """
    Find which runs are missing from the project and print the commands to add them in.

    Args:
        df_list (List[pd.DataFrame]): list of dataframes by initial aggregation.
        param (List[str], optional): Parameters to compare. Defaults to PARAM.
    """
    missing_list = []
    prefix = "python src/main.py "
    new_df_list = []
    for df in df_list:
        new_df_list.append(df[param])
    big_df = pd.concat(new_df_list)
    unique = big_df.groupby(param).size().reset_index().rename(columns={0: "count"})
    # print(unique)
    for df in df_list:
        for _, row in unique.iterrows():
            new_df = df.copy()
            for column_number, entry in enumerate(row):
                column = unique.columns[column_number]
                # print(i, filter_df.columns[column_number], entry)
                if column != "count":
                    new_df = new_df[new_df[column] == entry]
            if len(new_df) == 0 and len(df) != 0:
                # print("\n MISSING \n", df["index"].to_numpy()[0], "\n", row)
                command = prefix + "atm.mem=" + df["index"].to_numpy()[0]
                for i, par in enumerate(param):
                    command += " " + PARAM_HYDRA[i] + "=" + str(row[par])
                missing_list.append(command)

    for item in missing_list:
        print(item)
    print("MISSING: ", len(missing_list))


def aggregate_matches(
    summary_df: pd.DataFrame,
    filter_df: pd.DataFrame,
    results: List[str] = RESULTS,
    include_std_dev: bool = True,
    print_missing: bool = False,
) -> pd.DataFrame:
    """
    Aggregate the matches between two dataframes to find the
    mean and std devation of a set of results.

    Args:
        summary_df (pd.DataFrame): The summary df create by summary_table.
        filter_df (pd.DataFrame): The dataframe to filter by.
        results (List[str], optional): _description_. Defaults to RESULTS.
        include_std_dev (bool, optional): Whether to calculate standard devation. Defaults to True.
        print_missing (bool, optional): Whether to highlight missing runs from ensemble. Defaults to False.

    Returns:
        pd.DataFrame: Includes uncertainty.ufloat values if include_std_dev=True.
    """
    df_list = _aggregate_matches(summary_df, filter_df)
    results_d = {result: [] for result in results}
    member_l = []
    for df in df_list:
        members = len(df)
        member_l.append(members)
        for result in results_d:
            if include_std_dev and members != 0:
                results_d[result].append(ufloat(df[result].mean(), df[result].std()))
            else:
                results_d[result].append(df[result].mean())
    for result in results_d:
        filter_df[result] = results_d[result]
    filter_df["N"] = member_l
    if print_missing:
        find_missing(df_list)
    return filter_df


DEFAULT_MEM_LIST: List[str] = ["EEEE", "EECE", "EEEC", "EECC"]


def aggregate_table(
    project: str = DEFAULT_PROJECT, mem_list: List[str] = DEFAULT_MEM_LIST,
) -> pd.DataFrame:
    """
    Make aggregate table.

    Args:
        project (str, optional): _description_. Defaults to DEFAULT_PROJECT.
        mem_list (List[str], optional): _description_. Defaults to DEFAULT_MEM_LIST.

    Returns:
        pd.DataFrame: _description_
    """
    return aggregate_matches(summary_table(project=project), mems_to_df(mem_list))


def _add_change_column(
    df: pd.DataFrame, initial_variable: str = "trend_nino3.4 [degC]"
) -> Tuple[pd.DataFrame, str]:
    new_variable = "change " + initial_variable
    df[new_variable] = df[initial_variable] - df[initial_variable].to_numpy()[0]
    return df, new_variable


def _remove_row(df: pd.DataFrame, row_index: str = "EEEE") -> pd.DataFrame:
    df_new = df.copy()
    return df_new.drop(row_index)


def change_table(
    project: str = DEFAULT_PROJECT, mem_list: List[str] = DEFAULT_MEM_LIST,
) -> Tuple[pd.DataFrame, str]:
    """
        Return a table with the differences between ECMWF run and the different inputs

        Args:
            project (str, optional): Which project to read.
            Defaults to DEFAULT_PROJECT.
            mem_list (List[str], optional): What list of inputs to compare.
            Defaults to DEFAULT_MEM_LIST.

        Returns:
            Tuple[pd.DataFrame, str]: The change table,
    and the name of the new variable column.
    """
    table = aggregate_table(project=project, mem_list=mem_list)
    table, new_variable = _add_change_column(table)
    table = _remove_row(table)
    # print(table[new_variable])
    return table[[x for x in table.columns if x not in RESULTS]], new_variable


def output_fig_2_data(project: str = DEFAULT_PROJECT) -> Tuple[List[pd.DataFrame], str]:
    """
    Output the figure 2 data for plotting.

    Args:
        project (str, optional): Wandb project to read. Defaults to DEFAULT_PROJECT.

    Returns:
        Tuple[List[pd.DataFrame], str]: The change table,
        and the name of the new variable column.
    """
    table_list = []
    for mem_list in [
        ["EEEE", "EECE", "EEEC", "EECC"],
        ["EEEE", "EESE", "EEES", "EESS"],
    ]:
        table, new_variable = change_table(project=project, mem_list=mem_list)
        table_list.append(table)
    return table_list, new_variable


if __name__ == "__main__":
    # python src/wandb_utils.py
    # _other_tests()
    output_fig_2_data("sdat2/ENSOTrend-gamma")
    # print(summary_table(project="sdat2/ENSOTrend-beta"))
    # print(summary_table(project="sdat2/ENSOTrend-gamma"))
