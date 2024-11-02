<img src="https://github.com/Ajax23/PoreMS/blob/master/docsrc/pics/logo_text_sub.svg" width="60%">

--------------------------------------

[![PyPI Version](https://img.shields.io/badge/PyPI-0.3.0-orange)](https://pypi.org/project/porems/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/Ajax23/PoreMS/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14028652.svg)](https://doi.org/10.5281/zenodo.14028652)
[![Build Status](https://github.com/Ajax23/PoreMS/actions/workflows/workflow.yml/badge.svg)](https://github.com/Ajax23/PoreMS/actions/workflows/workflow.yml)
[![codecov](https://codecov.io/gh/Ajax23/PoreMS/branch/master/graph/badge.svg)](https://codecov.io/gh/Ajax23/PoreMS)

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
* Nguyen et al., 2024. Effects of Surfaces and Confinement on Formic Acid Dehydrogenation Catalyzed by an Immobilized Ru–H Complex: Insights from Molecular Simulation and Neutron Scattering. ACS Catalysis, doi:[doi.org/10.1021/acscatal.4c02626](https://doi.org/10.1021/acscatal.4c02626)
  - Data-Repository: doi:[]()
* Kraus et al., 2023. Axial Diffusion in Liquid-Saturated Cylindrical Silica Pore Models. The Journal of Physical Chemistry C, doi:[10.1021/acs.jpcc.3c01974](https://doi.org/10.1021/acs.jpcc.3c01974).
  - Data-Repository: doi:[10.18419/darus-3067](https://doi.org/10.18419/darus-3067)
* Kraus and Hansen, 2022. An atomistic view on the uptake of aromatic compounds by cyclodextrin immobilized on mesoporous silica. Adsorption, doi:[10.1007/s10450-022-00356-w](https://doi.org/10.1007/s10450-022-00356-w).
  - Data-Repository: doi:[10.18419/darus-2154](https://doi.org/10.18419/darus-2154)
* Kraus et al., 2021. PoreMS: a software tool for generating silica pore models with user-defined surface functionalisation and pore dimensions. Molecular Simulation, 47(4), pp.306-316, doi:[10.1080/08927022.2020.1871478](https://doi.org/10.1080/08927022.2020.1871478).
  - Data-Repository: doi:[10.18419/darus-1170](https://doi.org/10.18419/darus-1170)
* Ziegler et al., 2021. Confinement Effects for Efficient Macrocyclization Reactions with Supported Cationic Molybdenum Imido Alkylidene N-Heterocyclic Carbene Complexes. ACS Catalysis, 11(18), pp. 11570-11578, doi:[10.1021/acscatal.1c03057](https://doi.org/10.1021/acscatal.1c03057)
  - Data-Repository: doi:[10.18419/darus-1752](https://doi.org/10.18419/darus-1752)
* Kobayashi et al., 2021. Confined Ru-catalysts in a Two-phase Heptane/Ionic Liquid Solution: Modeling Aspects. ChemCatChem, 13(2), pp.739-746, doi:[10.1002/cctc.202001596](https://doi.org/10.1002/cctc.202001596).
  - Data-Repository: doi:[10.18419/darus-1138](https://doi.org/10.18419/darus-1138)
* Ziegler et al., 2019. Olefin Metathesis in Confined Geometries: A Biomimetic Approach toward Selective Macrocyclization. Journal of the American Chemical Society, 141(48), pp.19014-19022, doi:[10.1021/jacs.9b08776](https://doi.org/10.1021/jacs.9b08776).
  - Data-Repository: doi:[10.18419/darus-477](https://doi.org/10.18419/darus-477)
