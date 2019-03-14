import os
import shutil

from setuptools import find_packages, setup

import versioneer

NAME = "mykit"
DESCRIPTION = 'Utilites for manipulating inputs and outputs of various electronic structure calculator'
AUTHOR = 'Minye Zhang'
EMAIL = 'stevezhang@pku.edu.cn'
URL = 'https://github.com/minyez/' + NAME
REQUIRES_PYTHON = '>=3.6.0'
BINDIR = "tools"
REQUIRED = [
    'numpy', 'scipy', 'spglib',
]
LICENSE = "LICENSE"

# add executables in tools/ to scripts
SCRIPTS = []
tools = os.path.join(os.path.dirname(__file__), BINDIR)
for fn in os.listdir(tools):
    fpath = os.path.join(tools, fn)
    if os.access(fpath, os.X_OK):
        SCRIPTS.append(os.path.join(BINDIR, fn))

# start setup function
setup(
    name            = NAME,
    description     = DESCRIPTION,
    author          = AUTHOR,
    author_email    = EMAIL,
    url             = URL,
    packages        = [NAME],
    license         = LICENSE,
    version         = versioneer.get_version(),
    install_requires = REQUIRED,
    include_package_data = True,
    package_data = {
        NAME: ['vasp/metadata/*', 'core/metadata/*'],
    },
    python_requires = REQUIRES_PYTHON,
    #scripts         = SCRIPTS,
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy'
    ],
    cmdclass=versioneer.get_cmdclass(),
)

shutil.rmtree('__pycache__')
