# Snakemake workflow: CLimate Impact of AeroSOls 

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for analysing CMIP6 model output. This workflow is designed work with CMIP6
model data stored in DKRZ data format. 

## Installation

Before you can run the workflow the a conda environment need to be built first:

```
conda env create -f dustysnake -p ./dustysnake
```
Then activate the environment: 

```
conda activate ./dustysnake

``` 
## Usage

The workflow can be configured in the `config/config.yaml` file. Here which
experiment and activities to analyze can be defined. Current version has been primarily developed for analyzing AerChemMIP type of experiments.  

* Uses the `lookup_*.yaml` files are used to find the  CMIP6 files, to avoid to extensive globbing search. 
  To make it faster when working on mounted file systems. 

* Extending the workflow can be done trough adding additional rules. Following the standard snakemake format.

### Running the workflow:
Default is to run the `all` rule defined in the `Snakefile`:
```
snakemake -j2 
``` 

The `-j` argument specify how many cores the workflow will use. 
  
To print which rules are defined in the workflow use:
```
snakemake -l

calc_ERF_surf
calculate_ERF_TOA
calculate_SW_ERF
calculate_ERF_TOA_LW
calc_cloud_radiative_effect
calc_direct_radiative_effect
calc_absorption
calc_clim_PI_control
calc_experiment_climalogies
calc_feedback
plot_ERFs
plot_atm_abs
plot_global_avaraged_ERFs
plot_ERFaci
plot_change_cdnc
plot_change_lwp
plot_change_clivi
plot_change_clt
plot_change_ts
plot_change_tas
plot_change_prs
plot_emidust
plot_feedback_decomposed
plot_depdust
generate_table
all
clim_od550dust_aerocom_CMIP6
clim_od550aer_aerocom_CMIP6

``` 
To generate the decomposition tables use:
```
snakemake -j2 generate_table
```

## Mounting of CMIP6 storage server
This workflow relies on the having access to the Betzy and Nird storage system. Most of the 
CMIP6 data archived in Norway is stored at the Betzy, while the NorESM output is archived on the 
NIRD server. However it should also be possible to make the workflow work on any storage systems also as the data is archived following the DKRZ data format. 