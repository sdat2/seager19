"""Ocean model."""
import os
import time
import xarray as xr
import logging
from omegaconf import DictConfig
from typeguard import typechecked
import wandb
from src.visualisation.ani import animate_ds
from src.utils import timeit
from src.data_loading.ingrid import linear_qflx_replacement
from src.models.model_setup import ModelSetup


log = logging.getLogger(__name__)


class Ocean:
    """Ocean."""

    def __init__(self, cfg: DictConfig, setup: ModelSetup):
        """Init ocean.

        Args:
            cfg (DictConfig): config.
            setup (ModelSetup): The setup object.
        """
        self.setup = setup
        self.cfg = cfg

    def compile_all(self) -> None:
        """Compile the Fortran/C."""
        os.system("cd " + self.setup.ocean_src_path + " \npwd\nmake all")

    @typechecked
    def run(self, command: str) -> float:
        """Runs a line of bash in the ocean/RUN directory
        and times how long it takes.

        Args:
            command (str): a valid bash command.
            setup (ModelSetup): the model setup containing the file locations.
        """
        rc_prefix = "cd " + self.setup.ocean_run_path + " \n"
        full_command = rc_prefix + command
        ts = time.perf_counter()
        os.system(full_command)
        te = time.perf_counter()
        diff = te - ts
        print(full_command + " %2.5f s\n" % diff)
        return diff

    def edit_run(self):
        """
        Edit the run files to change ocean param.

        Args:
            cfg (DictConfig): [description]
            setup (ModelSetup): The setup script.
        """
        for i in ["om_spin", "om_diag", "om_run2f"]:
            file_name = os.path.join(self.setup.ocean_run_path, i)
            print("editing ", file_name)
            with open(file_name) as read_file:
                string_list = read_file.readlines()
            # read_file.close()
            # print(string_list)
            for j in range(len(string_list)):
                string_list[j] = string_list[j].replace(
                    "+NUMMODE              2",
                    "+NUMMODE              " + str(self.cfg.oc.nummode),
                )

            with open(file_name, "w") as write_file:
                write_file.writelines(string_list)

    @timeit
    @typechecked
    def run_all(self) -> None:
        """Run all the executables.

        Args:
            cfg (DictConfig): the model configs to pass.

        """
        self.edit_run()
        # Run the test to see if it's working.
        self.run("../SRC/" + self.cfg.ocean.tcom_name + " -i om_test")
        log.info("Run.")
        if self.cfg.ocean.spin:
            self.run(
                "../SRC/" + self.cfg.ocean.tcom_name + " -i om_spin -t om_spin.tios"
            )
            self.run("../SRC/" + self.cfg.ocean.tios2cdf_name + " -f output/om_spin")
            self.run("rm -rf output/om_spin.data output/om_spin.indx")
            self.run("cp -f output/om_spin.save output/om_spin.20y.restart")
        if self.cfg.ocean.diag:
            self.run(
                "../SRC/" + self.cfg.ocean.tcom_name + " -i om_diag -t om_diag.tios"
            )
            self.run("../SRC/" + self.cfg.ocean.tios2cdf_name + " -f output/om_diag")
            self.run("rm -rf output/om_diag.data output/om_diag.indx")
            self.run("cp -f output/om_diag.save output/om_diag.2y.restart")
        if self.cfg.ocean.ingrid:
            linear_qflx_replacement(self.setup)
        if self.cfg.ocean.run_through:
            run_time = self.run(
                "../SRC/" + self.cfg.ocean.tcom_name + " -i om_run2f -t om_run2f.tios"
            )
            wandb.log({"ocean_run": run_time})
            self.run("../SRC/" + self.cfg.ocean.tios2cdf_name + " -f output/om_run2f")
            self.run("rm -rf output/om_run2f.data output/om_run2f.indx")

    @timeit
    @typechecked
    def animate_all(self) -> None:
        """Animate the sst into gifs.

        Args:
            cfg (DictConfig): the model configs to pass.

        """
        l_x = list()
        if self.cfg.ocean.diag:
            l_x.append("om_diag")
        if self.cfg.ocean.run_through:
            l_x.append("om_run2f")

        for x in l_x:
            if "qflx" not in x:
                animate_ds(
                    xr.open_dataset(
                        str(os.path.join(self.setup.ocean_output_path, x)) + ".nc",
                        decode_times=False,
                    ),
                    x,
                    self.setup.direc,
                    front_trim=int("om_diag" == x),  # remove first month for om_diag
                )
            else:
                ds = xr.open_dataset(
                    str(os.path.join(self.setup.direc, x)) + ".nc", decode_times=False
                ).rename(
                    {
                        "T": "time",
                        "Y": "y",
                        "X": "x",
                        "Z": "Z",
                    }
                )
                animate_ds(ds, x, self.setup.direc)
