<img src="https://github.com/Ajax23/PoreMS/blob/master/docsrc/pics/logo_text_sub.svg" width="60%">

--------------------------------------

[![PyPI Version](https://img.shields.io/badge/PyPI-0.2.1-orange)](https://pypi.org/project/porems/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/Ajax23/PoreMS/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4525195.svg)](https://doi.org/10.5281/zenodo.4525195)
[![Build Status](https://travis-ci.com/Ajax23/PoreMS.svg?branch=master)](https://travis-ci.com/Ajax23/PoreMS)
[![codecov](https://codecov.io/gh/Ajax23/PoreMS/branch/master/graph/badge.svg)](https://codecov.io/gh/Ajax23/PoreMS)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/Ajax23/PoreMS.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Ajax23/PoreMS/context:python)

## Documentation

Online documentation is available at [ajax23.github.io/PoreMS](https://ajax23.github.io/PoreMS/).

<img src="https://github.com/Ajax23/PoreMS/blob/master/docsrc/pics/pore.svg" width="60%">

The docs include an example for generating [molecules](https://ajax23.github.io/PoreMS/molecule.html) and [pores](https://ajax23.github.io/PoreMS/pore.html), and an [API reference](https://ajax23.github.io/PoreMS/api.html). Visit [process](https://ajax23.github.io/PoreMS/process.html) for an overview of the programs operating principle.

An examplary [workflow](https://ajax23.github.io/PoreMS/workflow.html) has been provided for using the PoreMS package to create a pore system and run molecular dynamics simulation using [Gromacs](http://www.gromacs.org/).

## Dependencies

PoreMS supports Python 3.5+.

Installation requires [numpy](https://numpy.org/), [pandas](https://pandas.pydata.org/) and [matplotlib](https://matplotlib.org/).


## Installation

The latest stable release (and older versions) can be installed from PyPI:

    pip install porems

You may instead want to use the development version from Github:

    pip install git+https://github.com/ajax23/porems.git#egg=porems

    pip install git+https://github.com/ajax23/porems.git@develop#egg=porems

Or download the repository and install in the top directory via:

    pip install .


## Testing

To test porems, run the test in the test directory.


## Development

PoreMS development takes place on Github: [www.github.com/Ajax23/PoreMS](https://github.com/Ajax23/PoreMS)

Please submit any reproducible bugs you encounter to the [issue tracker](https://github.com/Ajax23/PoreMS/issues).


## How to Cite PoreMS

When citing PoreMS please use the following: **Kraus et al., Molecular Simulation, 2021, DOI: [10.1080/08927022.2020.1871478](https://doi.org/10.1080/08927022.2020.1871478)**

Additionaly, to assure reproducability of the generated pore systems, please cite the **Zenodo DOI** corresponding to the used PoreMS version. (Current DOI is listed in the badges.)

## Published Work
* Kraus et al., 2021. PoreMS: a software tool for generating silica pore models with user-defined surface functionalisation and pore dimensions. Molecular Simulation, pp.1-11, doi:[10.1080/08927022.2020.1871478](https://doi.org/10.1080/08927022.2020.1871478).
  - Data-Repository: doi:[10.18419/darus-1170](https://doi.org/10.18419/darus-1170)
* Kobayashi et al., 2020. Confined Ru‐catalysts in a Two‐phase Heptane/Ionic Liquid Solution: Modeling Aspects. ChemCatChem, 13(2), pp.739-746, doi:[10.1002/cctc.202001596](https://doi.org/10.1002/cctc.202001596).
  - Data-Repository: doi:[10.18419/darus-1138](https://doi.org/10.18419/darus-1138)
* Ziegler et al., 2019. Olefin Metathesis in Confined Geometries: A Biomimetic Approach toward Selective Macrocyclization. Journal of the American Chemical Society, 141(48), pp.19014-19022, doi:[10.1021/jacs.9b08776](https://doi.org/10.1021/jacs.9b08776).
  - Data-Repository: doi:[10.18419/darus-477](https://doi.org/10.18419/darus-477)
