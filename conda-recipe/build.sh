#!/bin/bash
$PYTHON setup.py install

# copy scripts over
mv ./bin/* "$PREFIX/bin/"
