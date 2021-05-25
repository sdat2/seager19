"""set_up.py"""
import os
from src.constants import (
    OCEAN_RUN_PATH,
    OCEAN_SRC_PATH,
    OCEAN_DATA_PATH,
    ATMOS_DATA_PATH,
    ATMOS_TMP_PATH,
)


class ModelSetup:
    """Initialise, store, and setup directories for models."""

    def __init__(self, direc: str):
        """
        Model directories init.

        Args:
            direc (str): The main model directory.
        """

        # setup ocean paths
        self.ocean_path = os.path.join(direc, "ocean")
        self.ocean_run_path = os.path.join(self.ocean_path, "RUN")
        self.ocean_src_path = os.path.join(self.ocean_path, "SRC")
        self.ocean_data_path = os.path.join(self.ocean_path, "DATA")
        self.ocean_output_path = os.path.join(self.ocean_path, "output")

        # setup atmospheric paths
        self.atmos_path = os.path.join(direc, "atmos")
        self.atmos_data_path = os.path.join(self.atmos_path, "DATA")
        self.atmos_tmp_path = os.path.join(self.atmos_path, "tmp")

        # make ocean paths
        os.mkdir(self.ocean_path)
        os.mkdir(self.ocean_src_path)
        os.mkdir(self.ocean_run_path)
        os.mkdir(self.ocean_data_path)
        os.mkdir(self.ocean_output_path)

        # make atmos paths
        os.mkdir(self.atmos_path)
        os.mkdir(self.atmos_data_path)
        os.mkdir(self.atmos_tmp_path)

        # make symlinks in ocean model
        os.symlink(self.ocean_data_path, os.path.join(self.ocean_run_path, "DATA"))
        os.symlink(self.ocean_data_path, os.path.join(self.ocean_src_path, "DATA"))
        os.symlink(self.ocean_output_path, os.path.join(self.ocean_run_path, "output"))
        os.symlink(self.ocean_output_path, os.path.join(self.ocean_src_path, "output"))

        self.init_ocean()
        self.init_atmos()

    def init_ocean(self):
        """initialise the ocean model."""
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

    def init_atmos(self):
        """Creating atmos."""
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
