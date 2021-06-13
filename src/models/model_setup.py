"""Set up the model, copy the files, get the names."""
import os
from omegaconf import DictConfig
from src.constants import (
    OCEAN_RUN_PATH,
    OCEAN_SRC_PATH,
    OCEAN_DATA_PATH,
    ATMOS_DATA_PATH,
    ATMOS_TMP_PATH,
)


class ModelSetup:
    """Initialise, store, and setup directories for models."""

    def __init__(self, direc: str, cfg: DictConfig) -> None:
        """
        Model directories init.

        Copies directories, relevant files, and creates

        Args:
            direc (str): The main model directory.
            cfg (DictConfig): The config.
        """

        # setup ocean paths
        self.direc = direc
        self.cfg = cfg
        self.ocean_path = os.path.join(direc, "ocean")
        self.ocean_run_path = os.path.join(self.ocean_path, "RUN")
        self.ocean_src_path = os.path.join(self.ocean_path, "SRC")
        self.ocean_data_path = os.path.join(self.ocean_path, "DATA")
        self.ocean_output_path = os.path.join(self.ocean_path, "output")

        # setup atmospheric paths
        self.atmos_path = os.path.join(direc, "atmos")
        self.atmos_data_path = os.path.join(self.atmos_path, "DATA")
        self.atmos_tmp_path = os.path.join(self.atmos_path, "tmp")

        for i in [
            # make ocean paths
            self.ocean_path,
            self.ocean_run_path,
            self.ocean_src_path,
            self.ocean_data_path,
            self.ocean_output_path,
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

    def _init_ocean(self):
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

    def _init_atmos(self):
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

    def ts(self, it: int, path: bool = True) -> str:
        name = int(it) + "_Q.nc"
        if path:
            return os.path.join(self.atmos_data_path, name)
        else:
            return name

    def ts_trend(self, it: int, path: bool = True) -> str:
        name = int(it) + "_Q.nc"
        if path:
            return os.path.join(self.atmos_data_path, name)
        else:
            return name

    def tau_base(self, it: int, path: bool = True) -> str:
        if it == 0:
            name = "tau-ECMWF"
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
