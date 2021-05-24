"""Ocean model."""
import os
import time
import xarray as xr
import wandb
import logging
from omegaconf import DictConfig
from typeguard import typechecked
from src.visualisation.ani import animate_ds
from src.utils import timeit
from src.data_loading.ingrid import linear_qflx_replacement
from src.models.model_setup import ModelSetup


log = logging.getLogger(__name__)


def compile_all(setup: ModelSetup) -> None:
    """Compile the Fortran/C."""
    os.system("cd " + setup.ocean_src_path + " \npwd\nmake all")


@typechecked
def run(command: str, setup: ModelSetup) -> None:
    """Runs a line of bash in the ocean/RUN directory
    and times how long it takes.

    Args:
        command (str): a valid bash command.
    """
    rc_prefix = "cd " + setup.ocean_run_path + " \n"
    full_command = rc_prefix + command
    ts = time.perf_counter()
    os.system(full_command)
    te = time.perf_counter()
    print(full_command + " %2.5f s\n" % ((te - ts)))


def edit_run(cfg: DictConfig, setup: ModelSetup):
    """
    Edit the run files to change ocean param.

    Args:
        cfg (DictConfig): [description]
        setup (ModelSetup): [description]
    """
    for i in ["om_spin", "om_diag", "om_run2f"]:
        file_name = os.path.join(setup.ocean_run_path, i)
        print("editing ", file_name)
        read_file = open(file_name)
        string_list = read_file.readlines()
        read_file.close()
        print(string_list)
        for j in range(len(string_list)):
            string_list[j] = string_list[j].replace(
                "+NUMMODE              2",
                "+NUMMODE              " + str(cfg.oc.nummode),
            )
        write_file = open(file_name, "w")
        write_file.writelines(string_list)


@timeit
@typechecked
def run_all(cfg: DictConfig, setup: ModelSetup) -> None:
    """Run all the executables.

    Args:
        cfg (DictConfig): the model configs to pass.

    """
    edit_run(cfg, setup)
    log.info("Run.")
    if cfg.ocean.spin:
        run("../SRC/" + cfg.ocean.tcom_name + " -i om_spin -t om_spin.tios", setup)
        run("../SRC/" + cfg.ocean.tios2cdf_name + " -f output/om_spin", setup)
        run("rm -rf output/om_spin.data output/om_spin.indx", setup)
        run("cp -f output/om_spin.save output/om_spin.20y.restart", setup)
    if cfg.ocean.diag:
        run("../SRC/" + cfg.ocean.tcom_name + " -i om_diag -t om_diag.tios", setup)
        run("../SRC/" + cfg.ocean.tios2cdf_name + " -f output/om_diag", setup)
        run("rm -rf output/om_diag.data output/om_diag.indx", setup)
        run("cp -f output/om_diag.save output/om_diag.2y.restart", setup)
    if cfg.ocean.ingrid:
        linear_qflx_replacement()
    if cfg.ocean.run_through:
        run("../SRC/" + cfg.ocean.tcom_name + " -i om_run2f -t om_run2f.tios", setup)
        run("../SRC/" + cfg.ocean.tios2cdf_name + " -f output/om_run2f", setup)
        run("rm -rf output/om_run2f.data output/om_run2f.indx", setup)


@timeit
@typechecked
def animate_all(cfg: DictConfig, setup: ModelSetup) -> None:
    """Animate the sst into gifs.

    Args:
        cfg (DictConfig): the model configs to pass.

    """
    l_x = list()
    if cfg.ocean.diag:
        l_x.append("om_diag")
    if cfg.ocean.run_through:
        l_x.append("om_run2f")

    for x in l_x:
        if "qflx" not in x:
            animate_ds(
                xr.open_dataset(
                    str(os.path.join(setup.ocean_output_path, x)) + ".nc",
                    decode_times=False,
                ),
                x,
                wandb.run.dir,
            )
        else:
            ds = xr.open_dataset(
                str(os.path.join(wandb.run.dir, x)) + ".nc", decode_times=False
            ).rename(
                {
                    "T": "time",
                    "Y": "y",
                    "X": "x",
                    "Z": "Z",
                }
            )
            animate_ds(ds, x, wandb.run.dir)
