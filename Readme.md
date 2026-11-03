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


It very important to check that the KORP and dunbrack binarie files are present with in the RCD installation.


## Quick start

A minimal example showing how to run AbRaCD on one input structure.

Its higly recomended to start with a structure with only the atoms of the H-L chain wanted. Otherwise, the presence of other chains or the antigen can cause malfuntioning of the program.

The Antibody can be in chothia or raw format.  

### Chothia run

The start and end point of the CDR is automatically detected:

```bash
bash chothia_onescript_complete  {file} {CDR}
```
Example:
```bash
bash chothia_onescript_complete  data/Benchmark/RosettaBenchmark/chothia/1nlb.pdb H3
```

Outputs will be generated in the results folder with the pdb name. As the critalografic structure iths knowed, also the RMSD will be availables


If you want to test some diferent secuence that the cristalografic one:
```bash
bash chothia_onescript_complete  {file} {CDR}   {sequence}
```
Example:

Example:
```bash
bash chothia_onescript_complete  data/Benchmark/RosettaBenchmark/chothia/1nlb.pdb  H3  ARFLRLQYFDV
```
Outputs will be generated in the results folder with the pdb name. In this ocasion, a structure with the best ranked decoy and a list  ranked decoys will be given. Its up of the user select wich ones wants.


`Posible script of use of insertpdb for select th ranked that we want? or all completed when alternative secuence is given`


### Raw run

The start and end point should be given

Be sure that the structure have all of the residues of each chain are linearly numbered (no numeric gaps before or after th target loop)


```bash
bash raw_onescript_complete {file} {CDR} {start} {end}
```
Example:

```bash
bash raw_onescript_complete  data/Benchmark/RosettaBenchmark/Processed/1nlb.pdb H3  97  107
```


```bash
bash raw_onescript_complete {file} {CDR} {start} {end}   {sequence}
```
Example:

```bash
bash raw_onescript_complete data/Benchmark/RosettaBenchmark/Processed/1nlb.pdb H3  97  107  ARFLRLQYFDV
```

## Run a group of cases together 




## Reproducing the paper results

Instructions for reproducing the benchmark results reported in the paper.

## Repository contents

A brief overview of the files and datasets included in this repository.

## Troubleshooting

Solutions to common problems encountered when installing or running AbRaCD.

## Status

This project is currently under development.

