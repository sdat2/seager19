# !/bin/bash
sphinx-apidoc -f -o . ..
rm setup.rst
rm modules.rst
rm main.rst
rm src.test.rst
