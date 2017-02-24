#!/bin/bash
set -euo pipefail
set -x

# Download and set up conda
curl -O https://fastdl.mongodb.org/linux/mongodb-linux-x86_64-ubuntu1204-3.4.1.tgz
tar -zxf mongodb-linux-x86_64-ubuntu1204-3.4.1.tgz && mv mongodb-linux-x86_64-ubuntu1204-3.4.1  ~/mongo

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

conda install -y python=3.5
conda install -y --file requirements.txt

~/anaconda/bin/python setup.py install
