"""Ocean model."""
import os
import time
import xarray as xr
import wandb
import logging
from omegaconf import DictConfig
from typeguard import typechecked
from src.constants import (
    OCEAN_DATA_PATH,
    OCEAN_OUTPUT_PATH,
    OCEAN_RUN_PATH,
    OCEAN_SRC_PATH,
)
from src.visualisation.ani import animate_ds
from src.utils import timeit
from src.data_loading.ingrid import linear_qflx_replacement


log = logging.getLogger(__name__)


def compile_all() -> None:
    """Compile the Fortran/C."""
    os.system("cd " + str(OCEAN_SRC_PATH) + " \npwd\nmake all")


@typechecked
def run(command: str) -> None:
    """Runs a line of bash in the ocean/RUN directory
    and times how long it takes.

    Args:
        command (str): a valid bash command.
    """
    rc_prefix = "cd " + str(OCEAN_RUN_PATH) + " \n"
    full_command = rc_prefix + command
    ts = time.perf_counter()
    os.system(full_command)
    te = time.perf_counter()
    print(full_command + " %2.5f s\n" % ((te - ts)))


@timeit
@typechecked
def run_all(cfg: DictConfig) -> None:
    """Run all the executables.

    Args:
        cfg (DictConfig): the model configs to pass.

    """
    log.info("Run.")
    if cfg.ocean.spin:
        run("../SRC/" + cfg.ocean.tcom_name + " -i om_spin -t spin.tios")
        run("../SRC/" + cfg.ocean.tios2cdf_name + " -f output/om_spin")
        run("rm -rf output/om_spin.data output/om_spin.indx")
        run("cp -f output/om_spin.save output/om_spin.20y.restart")
    if cfg.ocean.diag:
        run("../SRC/" + cfg.ocean.tcom_name + " -i om_diag -t diag.tios")
        run("../SRC/" + cfg.ocean.tios2cdf_name + " -f output/om_diag")
        run("rm -rf output/om_diag.data output/om_diag.indx")
        run("cp -f output/om_diag.save output/om_diag.2y.restart")
    if cfg.ocean.ingrid:
        linear_qflx_replacement()
    if cfg.ocean.run_through:
        run("../SRC/" + cfg.ocean.tcom_name + " -i om_run2f -t month.tios")
        run("../SRC/" + cfg.ocean.tios2cdf_name + " -f output/om_run2f")
        run("rm -rf output/om_run2f.data output/om_run2f.indx")


@typechecked
def copy_run_rdir(file_name: str) -> None:
    """Copy a file from the ocean/RUN directory.

    Args:
        file_name (str): file name to copy.

    """
    run("cp " + file_name + " " + str(wandb.run.dir))


@typechecked
def copy_output_rdir(file_name: str) -> None:
    """Copy a file from the ocean/output directory.

    Args:
        file_name (str): file name to copy.
    """
    os.system(
        "cd "
        + str(OCEAN_OUTPUT_PATH)
        + " \n cp "
        + file_name
        + " "
        + str(wandb.run.dir)
    )


def copy_data_rdir(file_name: str) -> None:
    """Copy a file from the ocean/DATA directory.

    Args:
        file_name (str): file name to copy.
    """
    os.system(
        "cd " + str(OCEAN_DATA_PATH) + " \n cp " + file_name + " " + str(wandb.run.dir)
    )


@typechecked
def copy_all(cfg: DictConfig) -> None:
    """Copy all relevant files to wandb folder.

    Args:
        cfg (DictConfig): the model configs to pass.

    """
    for x in [
        "om_diag",
        "om_diag.tr",
        "om_diag.log",
        "om_spin",
        "om_spin.tr",
        "om_spin.log",
        "om_run2f",
        "om_run2f.tr",
        "om_run2f.log",
        "month.tios",
        "spin.tios",
        "diag.tios",
    ]:
        copy_run_rdir(x)

    l_x = list()
    if cfg.ocean.spin:
        l_x.append("om_spin.nc")
    if cfg.ocean.diag:
        l_x.append("om_diag.nc")
    if cfg.ocean.run_through:
        l_x.append("om_run2f.nc")
    for x in l_x:
        copy_output_rdir(x)

    for x in [
        "qflx.nc",
        "qflx-0.nc",
    ]:
        copy_data_rdir(x)


@timeit
@typechecked
def animate_all(cfg: DictConfig) -> None:
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
                    str(os.path.join(wandb.run.dir, x)) + ".nc", decode_times=False
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
