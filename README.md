# Replication repository

This repository contains replication files the working paper "When Respondents Don't Care Anymore: Identifying the Onset of Careless Responding".

Below, we give a detailed explanation how all individual results from the paper are produced with the `R` scripts in this repository. We used the following machine to obtain all results: 

- A working station with an AMD EPYC 9654 96-Core Processor with 192 threads running Ubuntu 24.04 LTS.

This machine was using `R` version 4.4.2 and our package `carelessonset` version 0.1.0. This specific version can be installed with the following commands:

```R
install.packages("remotes")
remotes::install_github("my-anonymous-account/carelessonset")
```

If you already have package `remotes` installed, you can skip the first line.

Package `carelessonset` expects that you have `tensorflow`, `keras`, as well as the language Python (Version 3) installed. If not, then run:
```R
## install Python (version 3.10.12 here)
# install.packages("reticuate") # if package 'reticulate' isn't installed
reticulate::install_python(version = '3.10.12')

## Windows users should install version 3.10.11 because version 3.10.12 isn't
## available on Windows. To do so, uncomment the next line:
# reticulate::install_python(version = '3.10.11')

## tensorflow
install.packages("tensorflow")
tensorflow::install_tensorflow()

## keras
install.packages("keras")
keras::install_keras()
```
Alternatively, calling the function `carelessonset()` in package `carelessonset` will guide you through the installation process.


Please also note that for the simulations in the paper to run, you need package `simstudy` installed. You can install it via:

```R
install.packages("simstudy")
```

**Important note:** On some Mac platforms, we encountered errors in package `simstudy` when generating data. The issue may affect our scripts for the simulations and for the illustrative figures with simulated data. However, it seems to be an issue with package `simstudy` and not with our code. We did not encounter any issues with our replication scripts on Linux and Windows platforms.


## Illustrative figures

Figures 2--4 are illustrations for CODERS' results for simulated careless respondents. These figures are produced by the script `example_plots.R` in folder `illustration`.


## Simulations

Folders prefixed by `simulations_` contain replication files for the simulations.

- Raw results are always saved in subfolders named `results`. 
- Processed results ready to be plotted are always saved in subfolders named `analyzed`. 
- Figures are always saved in subfolders named `plots`.


### Simulations from the main text and with additional sample sizes

The folder `simulations_main` contains the materials for the simulations featured in the main text. It furthermore contains materials for additional results of traditional methods in this simulation design, as well as materials for the same design but with additional sample sizes (both are found in Appendix C).

Running the scripts therein---in the following steps---reproduces our figures.

1. Perform simulations and generate results:
   - `j05sim_main_careless.R`: for nonzero prevalence of carelessness
   - `j05sim_main_nocareless.R`: for when there is no carelessness

To produce Figures 6--9, then run the following scripts:

2. `j05_gatherresults_CODERS.R`: analyze and gather results for CODERS
3. `plotmaker_main.R`: read in the analyzed results and produce the figures

To produce Figures C1 and C2 in Appendix C (traditional methods), run the following scripts:

2. Analyze results for the traditional methods:
   - `j05_analyze_careless_benchmarks.R`: for nonzero prevalence of carelessness
   - `j05_analyze_nocareless_benchmarks.R`: for when there is no carelessness
3. Read in the analyzed results and produce the figures:
   - `plotmaker_benchmarks.R`: Figure C1
   - `plotmaker_benchmarks_subgroups.R`: Figure C2

To produce Figures C7--C10 in Appendix C (additional sample sizes), run the following scripts:

2. `j05_gatherresults_vary-n_CODERS.R`: analyze and gather results for CODERS
3. `plotmaker_main_vary-n.R`: read in the analyzed results and produce the figures


### Simulations on decreasing survey length (Appendix C)

The folder `simulations_limitation/shortsurvey` contains the materials for the simulations in Appendix C on decreasing survey length.

To produce Figures C3--C6 in Appendix C, run the following scripts in this order:

1. Perform simulations and generate results:
   - `j05sim_shortsurvey_careless.R`: for nonzero prevalence of carelessness
   - `j05sim_shortsurvey_nocareless.R`: for when there is no carelessness
2. `j05_shortsurvey_gatherresults_CODERS.R`: analyze and gather results for CODERS
3. `plotmaker_shortsurvey.R`: read in the analyzed results and produce the figures 


### Simulations on carelessness around midpoint (Appendix C)

The folder `simulations_limitation/middling` contains the materials for the simulation in Appendix C on careless respondents randomly choosing from the middle answer categories. 

Running the scripts therein---in the following steps---reproduces our figures.

1. `j05sim_limit_onlymiddling.R`: perform simulations and generate results

To produce Figures C11 and C12, then run the following scripts:

2. `j05sim_analyze_limit_onlymiddling_CODERS.R`: analyze results for CODERS
3. Produce figures:
   - `plotmaker_onlymiddling_CODERS_careful.R`: read in the analyzed results and produce Figure C11
   - `plotmaker_onlymiddling_CODERS_divergence.R`: read in the analyzed results and produce Figure C12

To produce Figure C13, run the following scripts:

2. `j05_analyze_limit_onlymiddling_benchmarks.R`: analyze results for the traditional methods
3. `plotmaker_onlymiddling_benchmarks.R`: read in the analyzed results and produce Figure C13
   

### Simulations on temporary carelessness (Appendix C)

The folder `simulations_additional/multiple_changepoints` contains the materials for the simulations in Appendix C on tmporary carelessness. 

Running the scripts therein---in the following steps---reproduces our figures.

1. `j05sim_multipleCP_careless.R`: perform simulations and generate results

To produce Figures C14 and C15, then run the following scripts:

2. `j05_gatherresults_multipleCP_CODERS.R`: analyze and gather results for CODERS
3. `plotmaker_multipleCP.R`: read in the analyzed results and produce Figure C14
C15

To produce Figures C16 and C17, run the following scripts:

2. `j05_analyze_multipleCP_benchmarks.R`: analyze results for the traditional methods
3. Produce figures:
   - `plotmaker_benchmarks_multipleCP.R`: read in the analyzed results and produce Figure C16
   - `plotmaker_benchmarks_multipleCP_subgroups.R`: read in the analyzed results and produce Figure C17


## Empirical application

The folder `application` contains replication files for the empirical application, with the following folder structure:

- `data`: folder containing the raw data and preprocessing scripts
- `results`: folder containing results in `.RData` files
- `plots`: figures featured in the main text and appendix
- `plots_individual_series`: plots of CODERS series for all flagged respondents (not featured in the paper)
- `plotfuns.R`: contains helper functions used for making plots

To reproduce our analyses, run the following scripts in this order:

- `0_antonym+reliability.R`: get results for antonym and reliability
- `1_apply-CODERS.R`: get results for CODERS
- `2_analyze-CODERS.R`: analyze CODERS results (stores information on all flagged respondents in `results.RData`)
- `3_calculate_cormat+figure5`: produce Figure 5 (correlation matrix and marginal response probabilities of cleaned data used in simulations)
- `4_plots-for-all-respondents`: produce plots of CODERS series for all flagged respondents (stored in subfolder `plots_individual_series`)

To produce the figures and tables, then run the following scripts:

- `table1.R`: produce Table 1
- `figure10.R`: produce Figure 10
- `table2.R`: produce Table 2
- `table3_figure11+D1.R`: produce Table 3, Figure 11 (both in main text), and Figure D1 (in Appendix D)
