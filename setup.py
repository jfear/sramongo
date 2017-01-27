#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = [i.strip() for i in open('requirements.txt').readlines()]

setup(
    name='sramongo',
    version='0.0.1',
    description="A package to download metadata from SRA/Biosample/Geo and dump into a mongo database.",
    author="Justin Fear",
    author_email='justin.fear@nih.gov',
    url='https://github.com/jfear/sra2mongo',
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    entry_points={
        'console_scripts':
        [
            'sra2mongo = sramongo.sra2mongo:main',
        ],
    },
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)

