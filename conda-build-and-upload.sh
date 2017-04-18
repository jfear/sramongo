#!/bin/bash

# Always build with conda to test any build issues, but only upload if we're on
# master branch and it's not a pull request.
export SRAMONGO_VERSION=$(./_version.py version)
export SRAMONGO_BUILD=$(./_version.py build)
conda build conda-recipe

if [[ $TRAVIS_BRANCH = "master" && $TRAVIS_PULL_REQUEST = "false" ]]; then
  conda install anaconda-client -y
  anaconda \
    -t $ANACONDA_TOKEN \
    upload \
    -u jfear \
    $(conda build --output conda-recipe)
fi
