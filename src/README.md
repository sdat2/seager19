# Source code folder

All re-usable source code for the project goes here.

The source folder is structured as follows:

```txt
src
|
├── __init__.py    <- Makes src a Python module
|
├── constants.py   <- Includes project wide constants for easy imports
|
├── main.py        <- The main file to run the model from.
│
├── xr_utils.py     <- Xarray specific utilities for project.
│
├── plot_utils.py    <- Matplotlib specific utilities for project.
|
├── utils.py       <- Miscellaneous utilities for project.
|
├── data_loading   <- Scripts to download or generate data
|
├── preprocessing  <- Scripts to turn raw data into clean data and features for modeling
│
├── models         <- Scripts to train models and then use trained models to make
│                     predictions
│
├── test          <- Scripts for unit tests of your functions
│
├── configs       <- Scripts for unit tests of your functions
│
└── visualisation  <- Scripts to visualise things (animations)
```
