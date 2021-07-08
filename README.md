# seager19

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

Model run results: <https://wandb.ai/sdat2/seager19>

Docker image for gfortran/gcc/cdf/conda: <https://hub.docker.com/repository/docker/sdat2/seager19>

![Coupling over iterations with c_d=2.25e-3 over tha Pacific with the land masked out in green](gifs/coupling_pac_mask.gif)

## Purpose

A repository to contain, analyse, and expand upon the model from:

### Seager et al. 2019, Nature Climate Change, Strengthening Tropical Pacific Zonal Sea Surface Temperature Gradient Consistent with Rising Greenhouse Gases

<https://doi.org/10.1038/s41558-019-0505-x>

### Summary of their paper

- There is a west-east warm-to-cold contrast in the Pacific.

- State of the art GCMs predict weakening contrast with GHG conc.

- In observations it has increased.

- Their simple linear equatorial coupled model agrees with observations when forced with ECMWF reanalysis product.

- When forced with CMIP5 ensemble relative humidity and wind speed their model can reproduce the bias, highlighting
    these fields as the mechanism through which it is created.

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

Specifically the time period 33:00-44:30.

### Code

The code and data were taken from a Columbia University website:

<a href='http://kage.ldeo.columbia.edu:81/SOURCES/.LDEO/.ClimateGroup/.PROJECTS/.PublicationsData/.Seager_etal_NCC-2019/'>
<img src='https://upload.wikimedia.org/wikipedia/en/thumb/f/f1/Columbia_University_shield.svg/1200px-Columbia_University_shield.svg.png', width='150'>
</a>

This website is blocked from some IP addresses, such as the Cambridge VPN (probably my fault),
and so you may need to try a different connection method if it doesn't load.

## Setup

### Get the github repository

```bash

    git clone https://github.com/sdat2/seager19.git

    cd seager19

```

### If you have root access and can install linux packages

```bash

    # Enforcing these versions is necessary to be able to use 70s fortran.
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
    
    # on jasmin you would need to make a new tmp directory in your home directory
    # to be able to install the quite large singularity environment e.g:

    mkdir /home/users/sithom/tmp

    # and then call the command to create the singularity object

    TMPDIR=/home/users/sithom/tmp SINGULARITY_CACHEDIR=/home/users/sithom/tmp singularity pull docker://sdat2/seager19:g4.8

    # and then to launch the envirnonment

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


    # values from other papers quoted in the introduction to Romps (2014) "Raleigh Damping in the Free Troposphere"
    python src/main.py atm.eps_days=1.25,1.8,2,2.5,3,5,10 -m

    python src/main.py atm.eps_days=2.1,1.75,1.9 -m

    # Seager91, Matsuno66, Yu97, Gill80, Chang82, Sugiyama09, Wu00
    # breaks at 5, 10: Sugiyama09, Wu00. Breaks through tau.
    # looks pretty strange at 3: Chang82.

    python src/main.py atm.k_days=7,8,9,10 -m

    python src/main.py atm.k_days=4,5,6 -m

    python src/main.py atm.k_days=16,17,18 -m

    # uncoupled run without syncing:
    python src/main.py name=it_1 coup.iterations=1 coup.c_d=2.25e-3 wandb=false

    # Other runs:

    python src/main.py name=cd_2.0 coup.c_d=2.0e-3

    python src/main.py name=cd_1.0 coup.c_d=1.0e-3

    # testing a different set of param to replicate paper

    python src/main.py archive_dir=/gws/nopw/j04/ai4er/users/sdat2/uc_logs name=pap_2 wandb=false coup.iterations=1 atm.e_frac=0.5 animate=false ocean.ingrid=false

    python src/main.py archive_dir=/gws/nopw/j04/ai4er/users/sdat2/uc_logs name=efrac2 wandb=false coup.iterations=1 atm.e_frac=2 animate=false ocean.ingrid=false

    python src/main.py archive_dir=/gws/nopw/j04/ai4er/users/sdat2/uc_logs name=efrac0.5_fix wandb=false coup.iterations=1 atm.e_frac=0.5 animate=false ocean.ingrid=false

    python src/main.py atm.mem=EEEE atm.vary_cloud_const=true

```

## Other handy commands for development of repo

### Moving old model runs

```bash

    mv -f logs/* /gws/nopw/j04/ai4er/users/sdat2/cd_logs/
        
    mv -f logs/* /gws/nopw/j04/ai4er/users/sdat2/eps_days_logs/

    mv -f k_days_logs /gws/nopw/j04/ai4er/users/sdat2/

```

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

   grep -R "depth" ocean/SRC

   grep -R "aml" ocean/SRC

   grep -R "TDEEP" ocean/SRC

   grep -R "HSFC" ocean/SRC


grep -R "rho" ocean/SRC # 1023 kg m-2


f0=54.6746d+00,f1=-.603459d+00,f2=1.09987d-02,f3=-6.167d-05,
```

### Get CMIP5 multimodel means

<http://kage.ldeo.columbia.edu:81/SOURCES/.LDEO/.ClimateGroup/.PROJECTS/.IPCC/.CMIP5/.dataset/.CMIP5/.MultiModelMeans/.MMM-v2.3/.historical/.Surface/data.cdf>

<http://kage.ldeo.columbia.edu:81/SOURCES/.LDEO/.ClimateGroup/.PROJECTS/.IPCC/.CMIP5/.dataset/.CMIP5/.MultiModelMeans/.MMM-v2.3/.rcp85/.Surface/data.cdf>
