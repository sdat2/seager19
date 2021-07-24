# src/models

```txt
models
|
├── __init__.py       <- Makes src a Python module
|
├── atmos.py          <- Get the data from dropbox through various functions.
|
├── coupling.py       <- Couple the ocean and atmosphere (supervisor for rest of models).
|
├── model_setup.py    <- The file structure class for the class.
|
├── ocean.py          <- The ocean model is run from here through `os.system`.
|
└── poly.py           <- Fit polynomials (with uncertainties).
```
