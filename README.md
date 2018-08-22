# Multi-level HIV model

This repository contains the source code and scripts for analysis for a
multi-level HIV model. The genome-level model is based on Kauffmans NK model,
the individual-level model is a hybrid between a deterministic ODE model
and a stochastic mutation and immune response scheme.
The population-level model uses Frasers HIV transmission and progression rates
in combination with a contact network based on Leungs binding-site model.

## compilation (Linux)

just run the command

`$ make`

from the project folder. Be sure to install the libraries `gsl` and `pthread`.
You will need the development versions of these libraries (usually called `<libname>-dev`).
The code is written in C++ standard C++11.

## running a simulation (Linux)

execute the command

`$ ./mlev-hiv-sim [-m <mode>] [-s <seed>] [-i <id>] [-n <popsize>] [-b <burnin>] [-l <length>]`

where <seed> must be a non-negative integer, and <id> is a string used for naming data files.
Valid modes are `plev` for the simulation of an epidemic, `cohort` for a cohort,
and `ilev` for a simulation of one individual with highly detailed output.
Using the flag -h will print a help message
data will be written to files in the `data` folder.

## output files

When the simulation has ended, the following files should be in the data/ folder:

- `stats-file-<id>.xml` -- file on-the-fly statistics about hosts and MHC alleles
- `herit-file-<id>.xml` -- file with data that can be used to calculate heritability (e.g. in R)
- `ancs-file-<id>.xml` -- file containing the entire ancestor trace
- `crosssectional-sample-file-<id>.xml` -- at regular (long) time intervals,
  a large set of hosts are written to this file
- `longitudinal-sample-file-<id>.xml` -- at regular (short) time intervals,
  a random host is written to this file
- `immresp-file-<id>.xml` -- file containing all immune responses
- TODO ditch the `file` in the filenames, ...


## jupyter notebooks for interpreting the data

A couple of jupyter notebooks can be found in the folder `scripts` for plotting
figures. You'll need [Jupyter](http://jupyter.org) to use them.

- `ilev-analysis.ipynb` -- for plotting individual-level simulations (`-m ilev`)
- `cohort-analysis.ipynb` -- for plotting cohort data (`-m cohort`)
- `plev-analysis.ipynb` -- for plotting population-level simulations (`-m plev`)

Only notebooks with cleared output are tracked in the repository (with extension `.ipynb.cln`).
In order to edit and run the notebooks with jupyter, first run

`$ make -f nbs.makefile`

When you are ready to commit changes to the notebooks, first run

`$ make -f cln.makefile`

## plotting the contact network

the p-level simulation outputs `graphviz` strings in the `DOT` language.
These can be transformed into a figure with the `contact-networks.ipynb`
notebook.
