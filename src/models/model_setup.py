"""Set up the model, copy the files, get the names."""
from typing import List
import os
import pandas as pd
from omegaconf import DictConfig
from src.constants import (
    OCEAN_RUN_PATH,
    OCEAN_SRC_PATH,
    OCEAN_DATA_PATH,
    ATMOS_DATA_PATH,
    ATMOS_TMP_PATH,
    MODEL_NAMES,
    VAR_DICT,
)


def mems_to_df(mem_list: List[str]) -> pd.DataFrame:
    """
    Turn a list of mems into a dataframe of inputs.

    Args:
        mem_list (List[str]): List of mem to turn into corresponding
            dataframe of inputs.

    Returns:
        pd.DataFrame: A dataframe

    Example:
        Work out what inputs a list of runs got::

            mems_to_df(["EEEE", "CCCC", "66E6"])
    """
    results_lol = []
    for i in VAR_DICT:
        results_lol.append([])
    for i, mem in enumerate(mem_list):
        for j in VAR_DICT:
            if len(mem) <= j:
                results_lol[j].append(MODEL_NAMES["E"])
            else:
                results_lol[j].append(MODEL_NAMES[mem[j]])
    return pd.DataFrame(
        data={VAR_DICT[i]: results_lol[i] for i in range(len(results_lol))},
        index=mem_list,
    )


def mem_to_dict(mem: str) -> dict:
    """
    Single mem variable to dictionary of inputs.

    Uses logic in `mems_to_df`.

    Args:
        mem (str): the mem input e.g "EEEE"

    Returns:
        dict: dictionary of inputs, e.g.
    """
    return {var: input_d[mem] for (var, input_d) in mems_to_df([mem]).to_dict().items()}


class ModelSetup:
    """Initialise, store, and setup directories for models."""

    def __init__(self, direc: str, cfg: DictConfig, make_move: bool = True) -> None:
        """
        Model directories init.

        Copies directories, relevant files, and creates

        Args:
            direc (str): The main model directory.
            cfg (DictConfig): The config object.
            make_move (bool): whether to move the files and make the folders.
                Defaults to True.
        """

        # setup ocean paths
        self.direc = direc
        self.cfg = cfg

        # setup general output paths
        self.gif_path = os.path.join(direc, "gifs")
        self.nino_data_path = os.path.join(direc, "nino_data")
        self.nino_plot_path = os.path.join(direc, "nino_plot")
        self.plot_path = os.path.join(direc, "plots")

        # setup ocean paths
        self.ocean_path = os.path.join(direc, "ocean")
        self.ocean_run_path = os.path.join(self.ocean_path, "RUN")
        self.ocean_src_path = os.path.join(self.ocean_path, "SRC")
        self.ocean_data_path = os.path.join(self.ocean_path, "DATA")
        self.ocean_output_path = os.path.join(self.ocean_path, "output")
        self.ocean_old_io_path = os.path.join(self.ocean_path, "old_io")

        # setup atmospheric paths
        self.atmos_path = os.path.join(direc, "atmos")
        self.atmos_data_path = os.path.join(self.atmos_path, "DATA")
        self.atmos_tmp_path = os.path.join(self.atmos_path, "tmp")

        # the different model names in a dict? - used by key from self.mem.
        # most of this data isn't available.
        self.names: dict = MODEL_NAMES

        # dict of variables that are read in.
        self.var: dict = VAR_DICT
        # temperature of the surface, cloud area fraction, surface wind, rel humidity.

        self.input_dict = mem_to_dict(self.cfg.atm.mem)

        if make_move:

            for i in [
                # make general paths
                self.gif_path,
                self.nino_data_path,
                self.nino_plot_path,
                self.plot_path,
                # make ocean paths
                self.ocean_path,
                self.ocean_run_path,
                self.ocean_src_path,
                self.ocean_data_path,
                self.ocean_output_path,
                self.ocean_old_io_path,
                # make atmos paths
                self.atmos_path,
                self.atmos_data_path,
                self.atmos_tmp_path,
            ]:
                if not os.path.exists(i):
                    os.mkdir(i)

            # make symlinks in ocean model

            for i, j in [
                [self.ocean_data_path, os.path.join(self.ocean_run_path, "DATA")],
                [self.ocean_data_path, os.path.join(self.ocean_src_path, "DATA")],
                [self.ocean_output_path, os.path.join(self.ocean_run_path, "output")],
                [self.ocean_output_path, os.path.join(self.ocean_src_path, "output")],
            ]:
                if not os.path.exists(j):
                    os.symlink(i, j)

            self._init_ocean()
            self._init_atmos()

    def _init_ocean(self) -> None:
        """initialise the ocean model by copying files over."""

        for file_ending in ["*.F", "*.c", "*.h", "*.inc", "*.mod", ".tios"]:

            os.system(
                "cd "
                + str(OCEAN_SRC_PATH)
                + " \n cp "
                + file_ending
                + " "
                + str(self.ocean_src_path)
            )

        os.system(
            "cd "
            + str(OCEAN_SRC_PATH)
            + " \n cp "
            + "Makefile"
            + " "
            + str(self.ocean_src_path)
        )

        os.system(
            "cd "
            + str(OCEAN_RUN_PATH)
            + " \n cp "
            + "*"
            + " "
            + str(self.ocean_run_path)
        )

        os.system(
            "cd "
            + str(OCEAN_DATA_PATH)
            + " \n cp "
            + "*"
            + " "
            + str(self.ocean_data_path)
        )

        os.system("cd " + str(self.ocean_data_path) + " \n make all")

    def _init_atmos(self) -> None:
        """Creating atmos by copying files over."""

        os.system(
            "cd "
            + str(ATMOS_DATA_PATH)
            + " \n cp "
            + "*"
            + " "
            + str(self.atmos_data_path)
        )
        os.system(
            "cd "
            + str(ATMOS_TMP_PATH)
            + " \n cp "
            + "*"
            + " "
            + str(self.atmos_tmp_path)
        )

    # Iteration 0 is the initial name, itertion Z+ returns
    # a new name. The name alone should be an option to
    # allow renaming to occur.

    # pylint: disable=missing-function-docstring
    def tcam_output(self, path: bool = True) -> str:

        name = (
            "S91"
            + "-hq"
            + str(self.cfg.atm.h_q)
            + "-prcp_land"
            + str(self.cfg.atm.prcp_land)
            + ".nc"
        )
        if path:
            return os.path.join(self.atmos_path, name)
        else:
            return name

    def dq_output(self, path: bool = True) -> str:
        name = "dQ.nc"
        if path:
            return os.path.join(self.atmos_path, name)
        else:
            return name

    def q_output(self, path: bool = True) -> str:
        name = "Q.nc"
        if path:
            return os.path.join(self.atmos_path, name)
        else:
            return name

    def ts_clim(self, it: int, path: bool = True) -> str:
        if it == 0:
            name = self.clim_file("ts", "clim", "", path=False)
            # name = self.clim_name(0, path=False)  # "ts-ECMWF-clim.nc"
        else:
            name = "ts-" + str(it) + "-clim.nc"

        if path:
            return os.path.join(self.atmos_data_path, name)
        else:
            return name

    def ts_clim60(self, it: int, path: bool = True) -> str:
        if it == 0:
            name = self.clim_file("ts", "clim", "60", path=False)
        else:
            name = "ts-" + str(it) + "-clim60.nc"

        if path:
            return os.path.join(self.atmos_data_path, name)
        else:
            return name

    def ts_trend(self, it: int, path: bool = True) -> str:
        if it == 0:
            name = self.clim_file("ts", "trend", "", path=False)
        else:
            name = "ts-" + str(it) + "-trend.nc"
        if path:
            return os.path.join(self.atmos_data_path, name)
        else:
            return name

    def tau_base(self, it: int, path: bool = True) -> str:
        if it == 0:
            name = self.clim_file("tau", "", "", path=False)[:-4]
        else:
            name = "it_" + str(it) + "_tau"
        if path:
            return os.path.join(self.ocean_data_path, name)
        else:
            return name

    def tau_y(self, it: int, path: bool = True) -> str:
        return self.tau_base(it, path=path) + ".y"

    def tau_x(self, it: int, path: bool = True) -> str:
        return self.tau_base(it, path=path) + ".x"

    def tau_clim_base(self, it: int, path: bool = True) -> str:
        if it == 0:
            name = self.clim_file("tau", "clim", path=False)[:-3]
        else:
            name = "it_" + str(it) + "_clim_tau"
        if path:
            return os.path.join(self.ocean_data_path, name)
        else:
            return name

    def tau_clim_y(self, it: int, path: bool = True) -> str:
        return self.tau_clim_base(it, path=path) + ".y"

    def tau_clim_x(self, it: int, path: bool = True) -> str:
        return self.tau_clim_base(it, path=path) + ".x"

    def dq_df(self, it: int, path: bool = True) -> str:
        if it == 0:
            name = "dQdf-sample.nc"
        else:
            name = "it_" + str(it) + "_dq_df.nc"

        if path:
            return os.path.join(self.ocean_data_path, name)
        else:
            return name

    def dq_dt(self, it: int, path: bool = True) -> str:
        if it == 0:
            name = "dQdT-sample.nc"
        else:
            name = "it_" + str(it) + "_dq_dt.nc"

        if path:
            return os.path.join(self.ocean_data_path, name)
        else:
            return name

    def om_run2f_nc(self) -> str:
        return os.path.join(self.ocean_output_path, "om_run2f.nc")

    def ecmwf_sfcwind(self) -> str:
        return os.path.join(self.atmos_data_path, "sfcWind-ECMWF-clim.nc")

    def om_mask(self) -> str:
        return os.path.join(self.ocean_data_path, "om_mask.nc")

    def nino_png(self, it: int) -> str:
        return os.path.join(self.nino_plot_path, "nino_" + str(it) + ".png")

    def nino_nc(self, it: int) -> str:
        return os.path.join(self.nino_data_path, "nino_" + str(it) + ".nc")

    def coupling_video(self, pac: bool = False, mask_land=False) -> str:
        name = "coupling"
        if pac:
            name += "_pac"
        if mask_land:
            name += "_mask"
        return os.path.join(self.gif_path, name + ".gif")

    def prcp_quiver_plot(self) -> str:
        return os.path.join(self.plot_path, "prcp_quiver.png")

    def tuq_trend_plot(self) -> str:
        return os.path.join(self.plot_path, "tuq_trends.png")

    def rep_plot(self, num: str, suffix: str = "") -> str:
        return os.path.join(self.plot_path, "fig_" + str(num) + suffix + ".png")

    def _get_clim_name(self, var_num: int) -> str:
        return self.names[self.cfg.atm.mem[var_num]]

    def clim_file(
        self, var_name: str, typ: str = "clim", appendage: str = "", path: bool = True
    ) -> str:
        name = (
            var_name + "-" + self.input_dict[var_name] + "-" + typ + appendage + ".nc"
        )
        if path:
            return os.path.join(self.atmos_data_path, name)
        else:
            return name

    def clim60_name(self, var_num: int, path: bool = True) -> str:
        return self.clim_file(self.var[var_num], "clim", "60", path=path)

    def clim_name(self, var_num: int, path: bool = True) -> str:
        return self.clim_file(self.var[var_num], "clim", "", path=path)
