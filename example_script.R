library(dplyr)
library(stats)
library(ggplot2)
source("call_funcs_v0224.R")
source("sims_generate_funcs_v0224.R")
source("bootstrap_funcs_v0224.R")

# Simulate data 
n_periods = 12
sim_data = gen_sim(seed=1, n_trt=250, n_ngh=250, n_periods=n_periods)
head(sim_data)
saveRDS(sim_data, "example_simulated_data.RData")

# Estimate effects according to model assumptions
or_covs = ps_covs = paste0("X", 1:4) # misspecified would be paste0("W", 1:4)
point_ests = get_est_info(data=sim_data, n_periods=n_periods, or_covs=or_covs, ps_covs=ps_covs, 
                          is_boot=F, wt_col="boot_wts", wt_norm0=T, y1_col="Y1")
boot_ests = collect_boots(data=sim_data, n_periods=n_periods, nboots=10, 
                    or_covs=or_covs, ps_covs=ps_covs, wt_norm0=T, block_col="Zip")
