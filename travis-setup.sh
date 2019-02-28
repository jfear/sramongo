#!/bin/bash
set -euo pipefail

# Download and set up conda
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/anaconda
export PATH=~/anaconda/bin:$PATH

# Add channels in the specified order.
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
conda config --add channels jfear

conda install -y --file requirements.txt

# we also want conda-build to build a conda package
conda install -y conda-build

~/anaconda/bin/python setup.py install
