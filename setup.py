#!/bin/env python
# -*- coding: utf-8 -*-
"""

    This file is part of norton_beer.

    :copyright: Copyright 2023 by the shs team, see AUTHORS

    All rights reserved.

"""


from __future__ import print_function

import os
import subprocess as sp
from setuptools import setup, find_packages


DISTNAME = 'norton_beer'
DESCRIPTION = ('Norton-Beer apodization with analytical Fourier transformation.')
MAINTAINER = 'Konstantin Ntokas'
MAINTAINER_EMAIL = 'k.ntokas@fz-juelich.de'
VERSION = '0.1.0.dev0'


def git_info():
    try:
        git_rev = sp.Popen(['git', 'describe', '--always', '--dirty'], stdout=sp.PIPE).communicate()[0].decode().strip()
    except Exception as ex:
        print("%s %s", type(ex), ex)
        print("no git revision number found!")
        git_rev = "???"
    try:
        git_count = sp.Popen(['git', 'rev-list', 'HEAD', '--count'], stdout=sp.PIPE).communicate()[0].decode().strip()
        git_count_str = str(int(git_count)).zfill(5)
    except Exception as ex:
        print("%s %s", type(ex), ex)
        print("no git version number found!")
        git_count_str = "???"
    return git_rev, git_count_str


def write_version_py(filename='norton_beer/version.py'):
    git_rev, git_count_str = git_info()
    version_string = \
        '# THIS FILE IS GENERATED FROM THE NORTON_BEER SETUP.PY\n' \
        'VERSION = "{}"\n' \
        'REVISION = "{}"\n' \
        'COMMIT = "{}"\n'.format(VERSION, git_rev, git_count_str)
    with open(os.path.join(os.path.dirname(__file__), filename), 'w') as vfile:
        vfile.write(version_string)


def _main():

    write_version_py()

    setup(
        name=DISTNAME,
        description=DESCRIPTION,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        version=VERSION,
        # test_suite="py.test",
        setup_requires=["pytest", "pytest-runner", "pytest-flake8", "pytest-cov"],
        install_requires=["numpy"],
        classifiers=[
            'Development Status :: 1 - planning',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: BSD License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Operating System :: Unix',
            'Operating System :: MacOS',
        ],
        packages=find_packages(),
        zip_safe=False,
    )


if __name__ == "__main__":
    _main()
