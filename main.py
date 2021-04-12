"""A file to run model runs from with wandb.

Basically a wrapper for bash commands.

Example:
   Usage of script::
       python3 main.py run_name=test22

"""
import os
import time
import xarray as xr
import wandb
import logging
import hydra
from omegaconf import DictConfig, OmegaConf
from src.constants import (
    FIGURE_PATH,
    OCEAN_DATA_PATH,
    OCEAN_OUTPUT_PATH,
    SRC_PATH,
    OCEAN_RUN_PATH,
    OCEAN_SRC_PATH,
    PROJECT_PATH,
)
from src.visualisation.ani import animate_ds
from src.utils import timeit
from src.data_loading.ingrid import linear_qflx_replacement


log = logging.getLogger(__name__)


def start_wandb(cfg: DictConfig) -> None:
    """Intialises wandb for run."""
    run_dir = os.path.join(PROJECT_PATH, "logs", cfg.run_name)
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)
    wandb.init(
        project=cfg.project,
        entity=cfg.user,
        dir=run_dir,
        save_code=True,
        name=cfg.run_name,
        notes=cfg.notes,
        config=cfg,
    )


def compile_all() -> None:
    """Compile the Fortran/C"""
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
    """run all the executables."""
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
    """Copy a file from the ocean/RUN directory."""
    run("cp " + file_name + " " + str(wandb.run.dir))


def copy_output_rdir(file_name: str) -> None:
    """Copy a file from the ocean/output directory."""
    os.system(
        "cd "
        + str(OCEAN_OUTPUT_PATH)
        + " \n cp "
        + file_name
        + " "
        + str(wandb.run.dir)
    )


def copy_data_rdir(file_name: str) -> None:
    """Copy a file from the ocean/DATA directory."""
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


def print_locations() -> None:
    """Check that the locations are right."""
    file_prefix = os.path.join(wandb.run.dir, wandb.run.name)
    print(file_prefix)
    figure_prefix = os.path.join(FIGURE_PATH, wandb.run.name)
    print(figure_prefix)
    qflx_path = os.path.join(OCEAN_DATA_PATH, "qflx.nc")
    print(qflx_path)
    print(OCEAN_DATA_PATH)
    print(OCEAN_OUTPUT_PATH)


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


@timeit
@hydra.main(config_path=os.path.join(SRC_PATH, "configs"), config_name="config")
def main(cfg: DictConfig):
    """The main function to run the model and animations.

    Args:
        cfg (DictConfig): The hyrda dict config from the wrapper.

    """
    print(OmegaConf.to_yaml(cfg))
    print(cfg)
    """
    start_wandb(cfg)
    compile_all()
    if cfg.run:
        run_all(cfg)
    copy_all()
    if cfg.animate:
        animate_all()
    """


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()
    # 2. Save model inputs and hyperparameters
    # wandb.config
    # config.learning_rate = 0.01
    # 3. Log metrics over time to visualize performance
    # wandb.log({"loss": loss})
    # wandb.log(cfd)
