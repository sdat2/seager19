"""A file to run model runs from with wandb.

Example:
   Usage::
       python3 main.py

"""
# Flexible integration for any Python script
import os
import time
import wandb
from src.constants import FIGURE_PATH, OCEAN_DATA_PATH, OCEAN_OUTPUT_PATH


def start_wandb() -> None:
    wandb.init(
        project="seager19",
        entity="sdat2",
        name="test_5",
        notes="test run, uploading files",
    )

start_wandb()

# 2. Save model inputs and hyperparameters
# wandb.config
# config.learning_rate = 0.01

os.system("echo Test!")
os.system("pwd")


def compile_all() -> None:
    os.system("cd ocean/SRC\npwd\nmake all")

compile_all()


def run(command: str) -> None:
    """Runs a line of bash in the ocean/RUN directory
    and times how long it takes.

    Args:
        command (str): a valid bash command.
    """
    rc_prefix = "cd ocean/RUN\n"
    full_command = rc_prefix + command
    ts = time.perf_counter()
    os.system(full_command)
    te = time.perf_counter()
    print(full_command + " %2.5f s\n" % ((te - ts)))


def run_all() -> None:
    """run all the FORTRAN files."""
    run("../SRC/tcom.exe -i om_spin -t spin.tios")
    run("../SRC/tios2cdf.exe -f output/om_spin")
    run("rm -rf output/om_spin.data output/om_spin.indx")
    run("cp -f output/om_spin.save output/om_spin.20y.restart")
    run("../SRC/tcom.exe -i om_diag -t diag.tios")
    run("../SRC/tios2cdf.exe -f output/om_diag")
    run("rm -rf output/om_diag.data output/om_diag.indx")
    run("cp -f output/om_diag.save output/om_diag.2y.restart")
    run("../SRC/tcom.exe -i om_run2f -t month.tios")
    run("../SRC/tios2cdf.exe -f output/om_run2f")
    run("rm -rf output/om_run2f.data output/om_run2f.indx")


def copy_src_rdir(file_name: str) -> None:
    run("cp " + file_name + " " + str(wandb.run.dir))


def copy_output_rdir(file_name: str) -> None:
    run("cd output\n cp " + file_name + " " + str(wandb.run.dir))

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
    copy_src_rdir(x)

for x in [
    "om_spin.nc",
    "om_diag.nc",
    "om_run2f.nc"
]:
   copy_output_rdir(x)

# Model training here

# 3. Log metrics over time to visualize performance
# wandb.log({"loss": loss})
# wandb.log(cfd)
file_prefix = os.path.join(wandb.run.dir, wandb.run.name)
print(file_prefix)
figure_prefix = os.path.join(FIGURE_PATH, wandb.run.name)
qflx_path = os.path.join(OCEAN_DATA_PATH, "qflx.nc")
print(qflx_path)
print(OCEAN_DATA_PATH)
print(OCEAN_OUTPUT_PATH)
