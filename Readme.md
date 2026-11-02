# AbRaCD (AntiBody Random Coodinate Descent)

## Overview

AbRaCD (Antibody Random Coordinate Descent) is a sampling-based method for predicting antibody complementarity-determining region (CDR) loop structures. It generates and ranks multiple candidate conformations using antibody-specific structural constraints, then refines the selected models. The [AbRaCD web server](https://abracd.iqf.csic.es) provides access to the method.

This repository contains the scripts, data, and instructions needed to run AbRaCD locally and reproduce the results reported in the associated paper. It is intended to make the evaluation workflow straightforward to run and inspect.

> Martín-Hernández, I., Chacón, P., & López-Blanco, J. R. (manuscript in preparation). *AbRaCD: Fast sampling-based prediction of antibody CDR loops*. Bioinformatics Advances. [Publication details to be added when final.]




## Requirements

- **Python 3.10 or later**, with `mpi4py`
```bash
python -m pip install mpi4py
```

- **Perl 5.32 or later**
- **MPI**, recommended
- **PyRosetta**, required. [Download PyRosetta](https://www.pyrosetta.org/downloads)

Tested PyRosetta releases:

- `PyRosetta4.Release.python310.ubuntu.release-347`
- `PyRosetta4.Release.python310.ubuntu 2025.50 (weekly)`
- `PyRosetta4.Debug.python312.linux 2024.39 (weekly)`

Other Python and PyRosetta versions may work, but have not been tested. Compatibility is more likely with commonly used versions; very old versions may cause issues.

## Installation

The repository includes Python functions that must be installed before running the refinement workflow. Installation instructions will be added once the setup has been finalized.

```bash
python -m pip install -e 
```
  
### RCD

AbRaCD is integrated into [RCD](https://github.com/chaconlab/RCD.git). Follow the installation instructions in the RCD repository.




## Quick start

A minimal example showing how to run AbRaCD on one input structure.

## Reproducing the paper results

Instructions for reproducing the benchmark results reported in the paper.

## Repository contents

A brief overview of the files and datasets included in this repository.

## Troubleshooting

Solutions to common problems encountered when installing or running AbRaCD.

## Status

This project is currently under development.

