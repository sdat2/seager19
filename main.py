"""A file to run model runs from with wandb.

Basically a wrapper for bash commands.

Example:
   Usage of script::
       python3 main.py

"""
import os
import time
import xarray as xr
import wandb
import hydra
from omegaconf import DictConfig, OmegaConf
from src.constants import (
    FIGURE_PATH,
    OCEAN_DATA_PATH,
    OCEAN_OUTPUT_PATH,
    SRC_PATH,
    OCEAN_RUN_PATH,
    OCEAN_SRC_PATH,
)
from src.visualisation.ani import animate_xr_da
from src.utils import timeit
from src.data_loading.ingrid import linear_qflx_replacement


def start_wandb() -> None:
    """Intialises wandb for run."""
    wandb.init(
        project="seager19",
        entity="sdat2",
        name="test_11",
        notes="test run, animating files",
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
    if cfg.run:
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


def rdict(index: int) -> dict:
    """Returns renaming dict for xarray.DataArray.

    Args:
        index (int): index on coords

    Returns:
        dict: renaming dict.
    """
    return {
        "T_0" + str(index): "time",
        "Y_0" + str(index): "y",
        "X_0" + str(index): "x",
        "L_0" + str(index): "Z",
    }


@timeit
def animate_all() -> None:
    """Animate the sst into gifs."""
    for x in [
        "om_diag",
        "om_run2f",
        "qflx",
        "qflx-0",
    ]:
        animate_ds(
            xr.open_dataset(str(wandb.run.dir / x) + ".nc", decode_times=False), x
        )


@timeit
def animate_ds(ds: xr.Dataset, file_name: str) -> None:
    """Animate the `xarray.Dataset`."""
    cmap_d = {
        "DYN_PRES": "delta",
        "SST_QFLX": "delta",
        "SST_SST": "sst",
        "qflx": "delta",
        "QFLX": "delta",
        "TDEEP_HTHERM": "sst",
        "TDEEP_TDEEP": "sst",
        "TDEEP_HMODEL": "sst",
    }
    unit_d = {"SST_SST": "$^{\circ}$C"}
    for y in ds.variables:
        y = str(y)
        if "X_" not in y:
            if "Y_" not in y:
                if "L_" not in y:
                    if "T_" not in y or "SST" in y:
                        if "GRID" != y:
                            print(y)
                            da = ds[y]
                            for key in da.coords:
                                num = str(key)[3]
                            da = da.rename(rdict(num))
                            if y in unit_d:
                                da.attrs["units"] = unit_d[y]
                            da.coords["x"].attrs["units"] = "$^{\circ}$E"
                            da.coords["y"].attrs["units"] = "$^{\circ}$N"
                            animate_xr_da(
                                da.where(da != 0.0).isel(Z=0),
                                video_path=os.path.join(
                                    wandb.run.dir, file_name + "_" + y + ".gif"
                                ),
                                vcmap=cmap_d[y],
                            )


@timeit
@hydra.main(config_path=os.path.join(SRC_PATH, "configs"), config_name="config")
def main(cfg: DictConfig):
    """The main function to run the model and animations.

    Args:
        cfg (DictConfig): The hyrda dict config from the wrapper.
    """
    print(OmegaConf.to_yaml(cfg))
    start_wandb()
    compile_all()
    run_all(cfg)
    copy_all()
    animate_all()


if __name__ == "__main__":
    main()
    # 2. Save model inputs and hyperparameters
    # wandb.config
    # config.learning_rate = 0.01
    # 3. Log metrics over time to visualize performance
    # wandb.log({"loss": loss})
    # wandb.log(cfd)
