#!/bin/bash

# Always build with conda to test any build issues, but only upload if we're on
# master branch and it's not a pull request.
export SRAMONGO_BUILD=$(git rev-list --count `git describe --tags --abbrev=0`..HEAD --count)
conda build conda-recipe

if [[ $TRAVIS_BRANCH = "master" && $TRAVIS_PULL_REQUEST = "false" ]]; then
  conda install anaconda-client -y
  anaconda \
    -t $ANACONDA_TOKEN \
    upload \
    -u jfear \
    $(conda build --output conda-recipe)
fi
