"""setup.py allows the installation of the project by pip."""
from typing import List
from setuptools import find_packages, setup

REQUIRED: List[str] = ["sithom"]

setup(
    name="src",
    version="0.0.1",
    author="Simon Thomas",
    author_email="sdat2@cam.ac.uk",
    description="ENSOTrend Project",
    url="https://github.com/sdat2/seager19",
    install_requires=REQUIRED,
    packages=find_packages(),
)
