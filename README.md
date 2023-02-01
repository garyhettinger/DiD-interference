# DiD-interference
This repository implements code from research on difference-in-differences methodology in the presence of interference (https://arxiv.org/abs/2301.06697).

## Example

An example simulated dataset is provided in `example_simulated_data.RData` from the case when both treatment and outcome are simulated from observed covariates W. 

An example script calling a simulation and running a bootstrapped doubly robust analysis for the Average Treatment Effects on the Treated (ATT) and Neighboring Control (ATN) are provided in `did_example_run.R`. The example is run with only 10 bootstrap replicates (nboots) since performance is not optimized in this repository, but 500 replicates is more preferable.

## Functions

Functions to fit outcome and propensity score models as well as use these models to derive estimates are found in `did_estimate_functions.R`. Functions to generate the simulations and to call the simulation functions with proper parameters are found in `did_sim_generation.R` and `did_get_simulated_data.R`, respectively.
