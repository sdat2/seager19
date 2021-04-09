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


# 1. Start a W&B run
wandb.init(
    project="seager19",
    entity="sdat2",
    name="test_4",
    notes="test run, no files changed",
)


# 2. Save model inputs and hyperparameters
# wandb.config
# config.learning_rate = 0.01

os.system("echo Test!")
os.system("pwd")
rc_prefix = "cd ocean/RUN\n"
os.system("cd ocean/SRC\npwd\nmake all")


def run(command: str):
    full_command = rc_prefix + command
    ts = time.perf_counter()
    os.system(full_command)
    te = time.perf_counter()
    print("{full_command}  %2.5f s\n" % ((te - ts)))


run("../SRC/tcom.exe -i om_spin -t spin.tios")
run("../SRC/tios2cdf.exe -f output/om_spin")
run("rm -rf output/om_spin.data output/om_spin.indx")
run("cp -f output/om_spin.save output/om_spin.20y.restart")
run("../SRC/tcom.exe -i om_diag -t diag.tios")
run("../SRC/tios2cdf.exe -f output/om_diag")
run("rm -rf output/om_diag.data output/om_diag.indx")
run("cp -f output/om_diag.save output/om_diag.2y.restart")
run("../SRC/tcom.exe -i om_run2f -t month.tios")
run("../SRC/tios2cdf.exe -f output/om_run2f ")
run("rm -rf output/om_run2f.data output/om_run2f.indx")

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
