import os
import shutil

from setuptools import find_packages, setup

NAME = "mykit"
DESCRIPTION = 'Utilites for manipulating inputs and outputs of various electronic structure calculator'
AUTHOR = 'Minye Zhang'
EMAIL = 'stevezhang@pku.edu.cn'
URL = 'https://github.com/minyez/' + NAME
REQUIRES_PYTHON = '>=3.6.0'

setup(
    name            = NAME,
    description     = DESCRIPTION,
    author          = AUTHOR,
    author_email    = EMAIL,
    url             = URL,
    # packages        = find_packages(exclude=("tests", "doc")),
    packages        = [NAME],
    license         = "LICENSE",
    version         = "0.0.1",
    include_package_data = True,
    package_data = {
        NAME: ['vasp/metadata/*', 'core/metadata/*'],
    },
    python_requires = REQUIRES_PYTHON,
)

# shutil.rmtree("build")
# shutil.rmtree("dist")
# shutil.rmtree(NAME+".egg-info")
# shutil.rmtree(os.path.join(NAME, NAME+".egg-info"))
