source("did_sim_generation.R")


################################################################################
# This file calls the sim generation functions with associated parameters
################################################################################


### get_outcome_params_by_period: creates outcome model parameters according to number of periods (2, 4, or 13)
### PARAM @n_time: number of time periods before tax (n_time=2,4,or 13)
### RETURN: list of n_time: int, pre_lambdas: n_var x n_time period-specific pre-tax effect parameters, 
#### post_lambdas: n_var x n_time period-specific post-tax effect parameters, gamma_p: vector of period-specific effects
get_outcome_params_by_period = function(n_time, period_effects) {
  pre_lambdas_base = c(0.5, 0.25, 0.25, 0.75)*3
  tau_zp = rbind(rep(0,n_time), rep(1,n_time), rep(-2,n_time), rep(0,n_time)) # can change for period-specific effects
  if (n_time == 2) {
    period_lambda_weights = c(0.95, 1.05)
    gamma_p = c(0, 0.5)
    period_effect_weights = c(0.9, 1.1)
  }
  # Based on seasonal ideas
  else if (n_time == 4) {
    period_lambda_weights = c(0.85, 1, 1.1, 0.95)
    gamma_p = c(0, 0.5, 1, 0.5)
    period_effect_weights = c(0.8, 1.05, 1.2, 0.95)
  }
  # Based on granular seasonal ideas
  else if (n_time == 13) {
    period_lambda_weights = c(0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.15, 1, 0.95, 0.9)
    gamma_p = c(0, 0.2, 0.3, 0.5, 0.6, 0.7, 0.9, 1.0, 1.1, 0.9, 0.6, 0.4, 0.3)
    period_effect_weights = c(0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15, 1.2, 1.15, 1, 0.95, 0.9)
  }
  pre_lambdas = c()
  for (i in 1:n_time) {
    pre_lambdas = cbind(pre_lambdas, period_lambda_weights[i]*pre_lambdas_base)
  }
  post_lambdas = 2 * pre_lambdas
  if (period_effects) tau_zp = tau_zp * period_effect_weights
  return(list(n_time=n_time, pre_lambdas=pre_lambdas, post_lambdas=post_lambdas, 
              gamma_p=gamma_p, tau_zp=tau_zp))
  
}

### get_full_data_set: compiles parameters in order to call simulation generation function 
### PARAM @n: int for number of total stores in study
### PARAM @ps_true: boolean saying if correctly specified the ps model
### PARAM @outcome_true: boolean saying if correctly specified the outcome model
### PARAM @nperiods: number of periods in simulation
### RETURN: simulated data object
get_full_data_set = function(n, ps_true, outcome_true, nperiods, period_effects) {
  simulate_A_from_W = if (ps_true) T else F
  simulate_Y_from_W = if (outcome_true) T else F
  params = get_outcome_params_by_period(n_time=nperiods, period_effects=period_effects)
  beta_a = rbind(c(0.1,0.2,0.5,0.7), c(0.1, 0.6, -0.4, 1.1), c(0.25, 0.8, -0.1, 0.3), c(0.2,-0.5, 0.4, -0.8))
  beta_0_a = c(-0.05, 0.1, 0.2, 0.4) #if (n <= 200) c(-0.3, 0.1, 1.3, 1.3) else c(-0.45, -0.2, 1.1, 1.1)
  alpha_a = c(1,2,6,0)
  n_zips = NULL #as.integer(c(n/(4*5), n/(4*5), n/(2*5))) # 1:1:2 store ratio, 5 stores per zip-code
  D_vars = c(2,3,4)
  D_var_signs = c(1,1,-1)#/3
  sim_data = get_simulated_data(n=n, n_zips=n_zips, n_time=params$n_time, n_stat_bin=1, n_stat_cont=2, n_tv=1,
                                beta_0_a=beta_0_a, beta_a=beta_a, 
                                pre_lambdas=params$pre_lambdas, 
                                post_lambdas=params$post_lambdas, 
                                alpha_a=alpha_a, gamma_p=params$gamma_p, gamma_y=-0.5, tau_zp=params$tau_zp,
                                D_vars=D_vars, D_var_signs=D_var_signs, simulate_A_from_W=simulate_A_from_W,
                                simulate_Y_from_W=simulate_Y_from_W, n_att=as.integer(n*0.43))
  if ((min(table(sim_data$data$city)) < 10*nperiods) | 
      (min(aggregate(W1 ~ city, FUN=function(x) length(unique(x)), 
                     data=unique(sim_data$data[,c("ID","W1","city")]))$V1) < 2)) return(NULL) # if simulation doesn't work for group, skip
  return(sim_data)
}

### get_simulated_data: gets simulated data, estimators for simulated data, and returns a list for individual simulation
### PARAM @seed: int for reproducible simulations
### PARAM @ps_true: boolean saying if correctly specified the ps model
### PARAM @outcome_true: boolean saying if correctly specified the outcome model
### PARAM @n: int for number of total stores in study
### PARAM @nperiods: number of periods in simulation
### PARAM @period_effects: boolean T if policy effects change over time
### RETURN: list of att effects, variances, lower, and upper CI bounds
get_simulated_dataset = function(seed, ps_true, outcome_true, n, nperiods, period_effects) {
  set.seed(seed)
  sim_data = NULL
  while(is.null(sim_data)) {
    sim_data = get_full_data_set(n=n, ps_true=ps_true, outcome_true=outcome_true, nperiods=nperiods, period_effects=period_effects)
    if (is.null(sim_data)) print("Bad City Assignments")
  }
  return(sim_data)
}


