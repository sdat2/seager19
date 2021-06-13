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

![SST output over spin up period](gifs/om_diag_SST_SST.gif)

## Purpose

A repository to contain and analyse the code from:

### Seager et al. 2019, Nature Climate Change, Strengthening Tropical Pacific Zonal Sea Surface Temperature Gradient Consistent with Rising Greenhouse Gases

<https://doi.org/10.1038/s41558-019-0505-x>

### Summary of paper

- There is a west-east warm-to-cold contrast in the Pacific.

- GCMs predict weakening contrast with GHG conc.

- In observations it has increased.

- Their simple linear model agrees with observations.

### Citation

The citation for this paper is:

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

The paper is discussed in a podcast available at:

<a href='https://deep-convection.org/2020/04/13/episode-5-richard-seager/'>
<img src='https://deep-convection.org/wp-content/uploads/2020/02/DC_logo_small_rectangular.png'
alt='Deep Convection: Season 1 Episode 5' width='150' />
</a>

## Code

The code and data was taken from:

<a href='http://kage.ldeo.columbia.edu:81/SOURCES/.LDEO/.ClimateGroup/.PROJECTS/.PublicationsData/.Seager_etal_NCC-2019/'>
<img src='https://upload.wikimedia.org/wikipedia/en/thumb/f/f1/Columbia_University_shield.svg/1200px-Columbia_University_shield.svg.png', width='150'>
</a>

## Videos of outputs

### ocean-model/RUN/output/

- SST_SST full in om_run2f: <https://youtu.be/JA97IWPmwxs>
- DYN_PRES in om_run2f: <https://youtu.be/5oRMWWAK1sM>
- TDEEP_HMODEL in om_run2f: <https://youtu.be/n25l6uYWEzY>
- TDEEP_HTHERM in om_run2f: <https://youtu.be/ikOo6VTXfkg>
- TDEEP_TDEEP in om_run2f: <https://youtu.be/BSRyTuESzLA>

## Setup

```bash
    git clone https://github.com/sdat2/seager19.git

    cd seager19

    make env

    conda activate ./env

    python3 src/data_loading/download.py

    sudo apt-get install gfortran

    sudo apt-get install gcc

    sudo apt-get install cloc

    cloc --report-file=docs/lang.txt $(git ls-files)

    make test        # Also downloads the data.

    make jupyter_pro

    make report

    docker build . -t sdat2/seager19:init

    docker push sdat2/seager19:init

    TMPDIR=/home/users/sithom/tmp SINGULARITY_CACHEDIR=/home/users/sithom/tmp singularity pull docker://sdat2/seager19:g4.8

    singularity run seager19_g4.8.sif

    conda init bash

    source ~/.bash_profile

    conda activate ./env/
    
```

### compilers in singularity container

```bash
    # x86_64-linux-gnu-gcc-4.8
    # x86_64-linux-gnu-gfortran-4.8
```
