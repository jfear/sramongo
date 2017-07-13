#!/bin/bash

# Build the dist files and upload to pypi
python setup.py sdist
python setup.py bdist_wheel

if [[ $TRAVIS_PULL_REQUEST = "false" && $TRAVIS_TAG ]]; then
    pip install twine -y
    twine upload -u $PYPI_USERNAME -p $PYPI_PASSWORD dist/*
fi
