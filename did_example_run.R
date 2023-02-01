source("did_get_simulated_data.R")
source("did_estimate_functions.R")

################################################################################
# This file provides an example script
################################################################################

# Get simulated data

sim_data = get_simulated_dataset(seed=1, ps_true=T, outcome_true=T, n=1000, nperiods=13, 
                              period_effects=F)

# Run doubly robust method

estims = get_estimates(data=sim_data$data, trt_col_name="city", trts=3, neighbors=2, 
                       trt_ctls=1, neighbor_ctls=4, nboots=10, period_col="period",
                       method="dr", ps_vars=c("W1", "W2", "W3", "W4"), 
                       or_vars=c("W1", "W2", "W3", "W4"), ps_ia_terms=c(), 
                       or_ia_terms=c())
estims$att
estims$atn
