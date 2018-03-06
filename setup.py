#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

requirements = [i.strip() for i in open('requirements.txt').readlines()]

setup(
    name='sramongo',
    version='0.0.3',
    description="A package to download metadata from SRA/Biosample/Geo and dump into a mongo database.",
    author="Justin M Fear",
    author_email='justin.m.fear@gmail.com',
    url='https://github.com/jfear/sramongo',
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

