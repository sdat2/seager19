"""Ocean model."""
import os
import time
import xarray as xr
import wandb
import logging
from omegaconf import DictConfig
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
def run_all(cfg: DictConfig) -> None:
    """Run all the executables.

    Args:
        cfg (DictConfig): the model configs to pass.

    """
    log.info("Run.")
    run("../SRC/tcom.exe -i om_spin -t spin.tios")
    run("../SRC/tios2cdf.exe -f output/om_spin")
    run("rm -rf output/om_spin.data output/om_spin.indx")
    run("cp -f output/om_spin.save output/om_spin.20y.restart")
    run("../SRC/tcom.exe -i om_diag -t diag.tios")
    run("../SRC/tios2cdf.exe -f output/om_diag")
    run("rm -rf output/om_diag.data output/om_diag.indx")
    run("cp -f output/om_diag.save output/om_diag.2y.restart")
    if cfg.ingrid:
        linear_qflx_replacement()
    run("../SRC/tcom.exe -i om_run2f -t month.tios")
    run("../SRC/tios2cdf.exe -f output/om_run2f")
    run("rm -rf output/om_run2f.data output/om_run2f.indx")


def copy_run_rdir(file_name: str) -> None:
    """Copy a file from the ocean/RUN directory.

    Args:
        file_name (str): file name to copy.
    """
    run("cp " + file_name + " " + str(wandb.run.dir))


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


def copy_all() -> None:
    """Copy all relevant files to wandb folder."""
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

    for x in [
        "om_spin.nc",
        "om_diag.nc",
        "om_run2f.nc",
    ]:
        copy_output_rdir(x)

    for x in [
        "qflx.nc",
        "qflx-0.nc",
    ]:
        copy_data_rdir(x)


@timeit
def animate_all() -> None:
    """Animate the sst into gifs."""
    for x in [
        "om_diag",
        "om_run2f",
        # "qflx",
        # "qflx-0",
    ]:
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
