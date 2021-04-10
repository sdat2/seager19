"""setup.py allows the installation of the project by pip."""
from setuptools import find_packages, setup

setup(
    name="src",
    version="0.0.1",
    author="Simon Thomas",
    author_email="sdat2@cam.ac.uk",
    description="ENSO MRes project 2020",
    url="https://github.com/sdat2/seager19",
    packages=find_packages(),
    # test_suite="src.tests.test_all.suite",
    # setup_requires=["pytest-runner"],
    # tests_require=["pytest"],
)
