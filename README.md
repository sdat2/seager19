# seager19

<a href="https://mybinder.org/v2/gh/sdat2/seager19/HEAD?filepath=notebooks%2Fbinder%2Fanalyse_results.ipynb">
<img alt="Binder" src="https://mybinder.org/badge_logo.svg"/></a>
<a href="https://opensource.org/licenses/MIT"><img alt="License: MIT" src=https://img.shields.io/badge/License-MIT-blue.svg></a>
 <a href="https://github.com/psf/black"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>
 <a href=https://www.python.org/downloads/release/python-388/><img src='https://img.shields.io/badge/python-3.8-blue.svg' alt='Python 3.8' /></a>
<a href='https://seager19.readthedocs.io/en/latest/?badge=latest'>
<img src='https://readthedocs.org/projects/seager19/badge/?version=latest' alt='Documentation Status' />
</a>
<a href='https://travis-ci.com/sdat2/seager19'>
    <img src='https://travis-ci.com/sdat2/seager19.svg?branch=main' alt='Build Status' />
</a>
<a href='https://coveralls.io/github/sdat2/seager19?branch=main'><img src='https://coveralls.io/repos/github/sdat2/seager19/badge.svg?branch=main' alt='Coverage Status' /></a>
<a href=https://github.com/sdat2/seager19/actions><img src='https://github.com/sdat2/seager19/actions/workflows/python.yml/badge.svg' alt='python' /></a>
<a href=https://github.com/sdat2/seager19/actions><img src='https://github.com/sdat2/seager19/actions/workflows/fort-c.yml/badge.svg' alt='fort-c' /></a>
<a href="https://zenodo.org/badge/latestdoi/309507184"><img src="https://zenodo.org/badge/309507184.svg" alt="DOI"></a>

Model run results: <https://wandb.ai/sdat2/seager19>

Docker image for gfortran/gcc/cdf/conda: <https://hub.docker.com/repository/docker/sdat2/seager19>

![Coupling over iterations with c_d=2.25e-3 over the Pacific, with the land masked out in green](gifs/coupling_pac_mask.gif)

## Purpose

A repository to contain, analyze, and expand upon the model from:

### Seager et al. 2019, Nature Climate Change, Strengthening Tropical Pacific Zonal Sea Surface Temperature Gradient Consistent with Rising Greenhouse Gases

<https://doi.org/10.1038/s41558-019-0505-x>

### Summary of their paper

- There is a west-east warm-to-cold contrast in the Pacific.

- State-of-the-art GCMs predict weakening contrast with GHG concentration.

- Their simple linear equatorial coupled model agrees
    with observations when forced with ECMWF reanalysis product.

- When forced with CMIP5 multi-model mean relative humidity and
    windspeed their model can reproduce the bias, showing
    that these mean state biases can explain the effect by
     changing the sensitivity of the ocean to warming.

### Citation

The citation for their paper is:

```bibtex
@article{seager2019strengthening,
  title={Strengthening tropical Pacific zonal sea surface temperature gradient consistent with rising greenhouse gases},
  author={Seager, Richard and Cane, Mark and Henderson, Naomi and Lee, Dong-Eun and Abernathey, Ryan and Zhang, Honghai},
  journal={Nature Climate Change},
  volume={9},
  number={7},
  pages={517--522},
  year={2019},
  url={https://doi.org/10.1038/s41558-019-0505-x},
  publisher={Nature Publishing Group}
}
```

### Podcast

The paper is discussed in a podcast available at Deep Convection Season 1 Episode 5:

<a href='https://deep-convection.org/2020/04/13/episode-5-richard-seager/'>
<img src='https://deep-convection.org/wp-content/uploads/2020/02/DC_logo_small_rectangular.png'
alt='Deep Convection: Season 1 Episode 5' width='150' />
</a>

Specifically they discuss the paper in the time-period 33:00-44:30.

### Code

The code and data were taken from this Columbia University website:

<a href='http://kage.ldeo.columbia.edu:81/SOURCES/.LDEO/.ClimateGroup/.PROJECTS/.PublicationsData/.Seager_etal_NCC-2019/'>
<img src='https://worldsummit.ai/wp-content/uploads/sites/4/2020/06/columbia-university-logo-png-columbia-university-crown.jpg', width='150'>
</a>

This website is blocked from some IP addresses, such as the Cambridge VPN (probably my fault),
and so you may need to try a different connection method if it doesn't load.

## Setup

### Get the Github repository

```bash

    git clone https://github.com/sdat2/seager19.git

    cd seager19

```

### If you have root access and can install Linux packages

```bash

    # Enforcing these versions is necessary to be able to use 70s Fortran.
    # If you would rather not change your version, you might be better off
    # Using the docker image instead.

    sudo apt-get install gfortran-4.8

    sudo apt-get install gcc-4.8

    sudo apt-get install --fix-missing libnetcdf-dev libnetcdff-dev

    sudo apt-get install cloc  # brew install cloc for mac
```

### If you want to make the docker environment yourself

```bash

    docker build . -t sdat2/seager19:g4.8

    docker push sdat2/seager19:g4.8
```

### If you need to install the singularity environment

```bash

    # On Jasmin you would need to make a new `tmp` directory in your home directory
    # to be able to install the quite large singularity environment e.g.:

    mkdir /home/users/sithom/tmp

    # and then call the command to create the singularity object

    TMPDIR=/home/users/sithom/tmp SINGULARITY_CACHEDIR=/home/users/sithom/tmp singularity pull docker://sdat2/seager19:g4.8

    # and then to launch the environment

    singularity run seager19_g4.8.sif

```

### Making the environment and testing it works (either in singularity or not)

```bash

    # in case conda isn't activated:

    conda init bash

    source ~/.bash_profile

    # activate the environment

    conda activate ./env/

    # test that the environment and package work before using:

    make test        # Also downloads the data if needed.

```

### Add optional features

```bash

    make jupyter_pro  # adds timing etc.

    make vscode_pro  # sets up nice autodoc, pylint etc.

    make jupyter_dark  # dark mode for jupyter notebooks

    # and to reverse this:
    make jupyter_light

```

### Examples of running the model

```bash

    # default run:

    python src/main.py name=cd_2.25 coup.c_d=2.25e-3

    # Sweep through different levels of Raleigh friction

    # values near the value in the paper
    python src/main.py atm.eps_days=1.05,0.95,0.85,0.75,0.65,0.55 -m

    python src/main.py atm.eps_days=1.7,1.8,1.85 -m 

    python src/main.py atm.eps_days=0.7,0.8 -m 


    # values from other papers quoted in the introduction 
    # to Romps (2014) "Raleigh Damping in the Free Troposphere"
    python src/main.py atm.eps_days=1.25,1.8,2,2.5,3,5,10 -m

    python src/main.py atm.eps_days=2.1,1.75,1.9 -m

    # Seager91, Matsuno66, Yu97, Gill80, Chang82, Sugiyama09, Wu00
    # breaks at 5, 10: Sugiyama09, Wu00. Breaks through tau.
    # looks pretty strange at 3: Chang82.

    python src/main.py atm.k_days=7,8,9,10 -m

    # uncoupled run without syncing:
    python src/main.py name=it_1 coup.iterations=1 coup.c_d=2.25e-3 wandb=false
```

To look at the commands to replicate the paper figures, see `replicate.md`.

## Other handy commands for development of repository

### Make a notebook with helpful magic functions for dark mode

```bash

    make notebook name=your-notebook-name
```

### New python script

```bash

    make py name=src/to-new.py
```

### Update docs of the src directory

```bash

    make autodoc
```

### Check where a variable is referenced in the ocean model

```bash

   grep -R "f1prime" ocean/SRC
```

### Get CMIP5 multimodel means from Columbia

<http://kage.ldeo.columbia.edu:81/SOURCES/.LDEO/.ClimateGroup/.PROJECTS/.IPCC/.CMIP5/.MultiModelMeans/.MMM-v2.3/.historical/.Surface/data.cdf>

<http://kage.ldeo.columbia.edu:81/SOURCES/.LDEO/.ClimateGroup/.PROJECTS/.IPCC/.CMIP5/.MultiModelMeans/.MMM-v2.3/.rcp85/.Surface/data.cdf>

### CMIP6/CMIP5 bias mechanism data

<https://docs.google.com/spreadsheets/d/1QrCLil7uHMRJECOoSL18uk2mvBwMLxqesXbZpUWh3ko/edit?usp=sharing>
