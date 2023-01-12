library(MASS)
library(nnet)
library(lme4)

################################################################################
# This file provides a simulated data set based on a number of user-specified options
################################################################################

### simulate_tv_X: function to simulate time_varying X variables
### PARAM @baseline_tv_X: matrix with one column per tv X variables
### PARAM @n_time: number of pre tax periods
### PARAM @n_stat: number of static continuous and binary variables
### PARAM @n_tv: number of time-varying continuous variables
### RETURN: data.frame with columns for time-varying X variables
simulate_tv_X = function(baseline_tv_X, n_time, n_stat) {
  time_multiplier = 0.1
  n = nrow(baseline_tv_X)
  n_tv = ncol(baseline_tv_X)
  tv_X = matrix(nrow=n_time*n, ncol=n_tv)
  
  # For each tv variable, set it so it is increasing on average over time
  for (i in 1:n_tv) {
    spec_tv_X = prev_spec_tv_X = baseline_tv_X[,i]
    for (j in 2:n_time) {
      prev_spec_tv_X = sapply(prev_spec_tv_X, function(x) rnorm(1, mean=(x + 0.25), sd=0.25))
      spec_tv_X = c(spec_tv_X, prev_spec_tv_X)
    }
    tv_X[,i] = spec_tv_X
  }
  tv_X = cbind(data.frame(ID=rep(1:n, n_time), period=rep(1:n_time, each=n)), data.frame(tv_X))
  names(tv_X) = c("ID", "period", paste0("X", (n_stat+1):(n_stat+n_tv)))
  return(tv_X)
}

### simulate_X: function to simulate X variables
### PARAM @n: total number of stores in study
### PARAM @n_time: number of periods
### PARAM @n_stat_bin: number of static binary variables
### PARAM @n_stat_cont: number of static continuous variables
### PARAM @n_tv: number of time-varying continuous variables
### RETURN: data.frame with columns for X variables
simulate_X = function(n, n_time, n_stat_bin, n_stat_cont, n_tv) {
  stat_stat_corr = 0 #0.4
  stat_tv_corr = 0 #0.3
  tv_tv_corr = 0 #0.5
  # bin_threshs = rep(0, n_stat_bin)
  n_stat = n_stat_bin + n_stat_cont
  n_var = n_stat_bin + n_stat_cont + n_tv
  Sigma = matrix(nrow=n_var, ncol=n_var)
  Sigma[1:n_stat, 1:n_stat] = stat_stat_corr
  Sigma[(n_stat+1):n_var,1:n_stat] = stat_tv_corr
  Sigma[1:n_stat, (n_stat+1):n_var] = stat_tv_corr
  Sigma[(n_stat+1):n_var,(n_stat+1):n_var] = tv_tv_corr
  diag(Sigma) = rep(1, n_var)
  X = mvrnorm(n, mu=rep(0,n_var), Sigma=Sigma)
  # for (i in 1:n_stat_bin) {
  #   X[,i] = as.integer(X[,i] >= bin_threshs[i])
  # }
  tv_X = simulate_tv_X(matrix(X[,(n_stat+1):n_var], ncol=n_tv), n_time, n_stat)
  stat_X = cbind(data.frame(ID=1:n), data.frame(X[,1:n_stat]))
  X = merge(stat_X, tv_X, by="ID")
  return(X)
}

### simulate_A: function to simulate treatment
#### 1:non-neighboring control (e.g., Baltimore), 
#### 2:neighboring control (e.g. Border counties)
#### 3:treatment (e.g. Philadelphia)
#### 4:non-neighboring control for neighbor (e.g. Non-border counties)
### PARAM @n: total number of stores in the study
### PARAM @data: data.frame with columns for X variables
### PARAM @beta_0_a: treatment specific intercept terms in multinomial model
### PARAM @beta_a: treatment specific covariate coefficient terms in multinomial model
### PARAM @n_var: total number of confounders
### PARAM @simulate_from_W: boolean that is true if observed data constitutes DGM
### RETURN: data.frame with columns for X variables
simulate_A_multinomial = function(n, data, beta_0_a, beta_a, n_var, simulate_from_W) {
  cov_prefix = if (simulate_from_W) "W" else "X"
  vars = paste0(cov_prefix,1:n_var)
  A = matrix(nrow=n, ncol=length(beta_0_a))
  tmp = matrix(beta_0_a, nrow=n, ncol=length(beta_0_a), byrow=T)
  baseline_data = aggregate(data[,vars], by=list(ID=data$ID), FUN=mean)
  for (k in 1:n_var) {
    tmp = tmp + baseline_data[[paste0(cov_prefix,k)]] %*% t(beta_a[k,])
  }
  probs = exp(tmp)
  for(i in 1:n){
    A[i,] = rmultinom(n=1, size=1, prob=probs[i,])
  }
  A = apply(A, 1, function(x) which(x == 1))
  A = data.frame(ID=baseline_data$ID, city=A)
  return(A)
}

simulate_A_binomial = function(n, data, beta_0_a, beta_a, n_var, simulate_from_W, labels) {
  cov_prefix = if (simulate_from_W) "W" else "X"
  vars = paste0(cov_prefix,1:n_var)
  A = rep(NA, n)
  tmp = rep(beta_0_a[2]-beta_0_a[1])
  baseline_data = aggregate(data[,vars], by=list(ID=data$ID), FUN=mean)
  for (k in 1:n_var) {
    tmp = tmp + baseline_data[[paste0(cov_prefix,k)]] * (beta_a[k,2] - beta_a[k,1])
  }
  probs = exp(tmp) / (1 + exp(tmp))
  for(i in 1:n){
    A[i] = rbinom(n=1, size=1, prob=probs[i])
  }
  A = ifelse(A == 1, labels[2], labels[1])
  A = data.frame(ID=baseline_data$ID, city=A)
  return(A)
}

simulate_A = function(n, data, beta_0_a, beta_a, n_var, simulate_from_W, n_att) {
  if (n_att > 0) {
    A_att = simulate_A_binomial(n_att, data[data$ID <= n_att,], beta_0_a[c(1,3)], 
                                beta_a[,c(1,3)], n_var, simulate_from_W, c(1,3))
    A_atn = simulate_A_binomial(n-n_att, data[data$ID > n_att,], beta_0_a[c(2,4)], 
                                beta_a[,c(2,4)], n_var, simulate_from_W, c(2,4))
    return(rbind(A_att, A_atn))
  }
  return(simulate_A_multinomial(n, data, beta_0_a, beta_a, n_var, simulate_from_W))
}

### simulate_D: function to simulate distance variable from W
### PARAM @n: total number of stores in the study
### PARAM @data: data.frame with columns for W variables
### PARAM @vars: variables affecting distance
### PARAM @var_signs: sign of variable coefficient in distance model
### RETURN: data.frame with columns for ID and D variables
simulate_D = function(n, data, vars, var_signs, simulate_from_W) {
  baseline_data = data[data$period == 1,]
  D_pred = 0
  cov_prefix = if (simulate_from_W) "W" else "X"
  
  # e.g. D* ~ N(plogis(W1 + W2 - W3)), D = trunc(D*, top=0.9, bottom=0.1)
  for (i in 1:length(vars)) {
    var = vars[i]
    var.sign = var_signs[i]
    D_pred = D_pred + var.sign * baseline_data[,paste0(cov_prefix, var)]
  }
  D = ifelse(baseline_data$city == 2, pmin(0.9, pmax(0.1, rnorm(n=n, mean=plogis(D_pred), sd=rep(0.05,n)))), rep(0, n))^2
  D = data.frame(ID=baseline_data$ID, D=D)
  return(D)
}

### simulate_Y: function to simulate outcome data
### PARAM @n: number of total stores in the study
### PARAM @data: data.frame with columns for time and covariate variables
### PARAM @n_time: number of periods
### PARAM @n_var: number of total covariates
### PARAM @pre_lambdas: coefficients of X/W in linear outcome model in pre-tax period (dim: n_var x n_time)
### PARAM @post_lambdas: coefficients of X/W in linear outcome model in post-tax period (dim: n_var x n_time)
### PARAM @alpha_a: Unobserved, time-invariant effects specific to treatment group (length: 3=n_trt)
### PARAM @gamma_p: period-specific time effects (length: n_time)
### PARAM @gamma_y: time-specific effect that changes pre and post tax (length: 1)
### PARAM @tau_z: (tau(1)=0, tau(2) = ATN, tau(3) = ATT, length: 3=n_trt)
### PARAM @simulate_from_W: boolean representing if observed covariates are true DGM or not
### RETURN: data.frame with added columns for PreTaxTarget and PostTaxTarget
simulate_Y = function(n, data, n_time, n_var, pre_lambdas, post_lambdas, alpha_a, 
                      gamma_p, gamma_y, tau_zp, simulate_from_W) {
  cov_prefix = if (simulate_from_W) "W" else "X"
  X0 = 20
  alpha_i = data.frame(ID=1:n, alpha_i=rnorm(n=n, mean=0, sd=1)) # individual specific intercepts
  data = merge(data, alpha_i, by="ID")
  for (period in 1:n_time) {
    sub_data = data[data$period == period,]
    pre_tmp = X0 + sub_data$alpha_i + alpha_a[sub_data$city] + gamma_p[period]
    post_tmp = pre_tmp + gamma_y + tau_zp[cbind(sub_data$city, sub_data$period)]*(1-sub_data$D)
    for (i in 1:n_var) {
      pre_tmp = pre_tmp + sub_data[[paste0(cov_prefix,i)]] * pre_lambdas[i,period] 
      post_tmp = post_tmp + sub_data[[paste0(cov_prefix,i)]] * post_lambdas[i,period] 
    }
    data[data$period == period, "PreTaxTarget"] = pre_tmp + rnorm(nrow(sub_data), mean=0, sd=0.5)
    data[data$period == period, "PostTaxTarget"] = post_tmp + rnorm(nrow(sub_data), mean=0, sd=0.5)
  }
  return(data)
}

### get_true_effects: function to collect true ATT and ATN from simulation parameters and data
### PARAM @data: data.frame with columns for X variables (now unused since calculate an effective tau_z)
### PARAM @tau_z: tau parameters in outcome model
### RETURN: list with true ATT and ATN for simulation
get_true_effects = function(data, tau_zp) {
  true.ATT = mean(tau_zp[3,])
  true.ATN = mean(tau_zp[2,]) #mean(tau_z[2] * (1-data[data$A == 2, "D"]))
  return(list(ATT=true.ATT, ATN=true.ATN))
}

### get_W: function to transform X into observed data
### PARAM @data: data.frame with columns for X variables
### PARAM @n_stat_bin: number of static binary variables
### PARAM @n_stat_cont: number of static continuous variables
### PARAM @n_tv: number of time-varying continuous variables
### PARAM @n_zips: vector of number of zip-codes (and thus unique static continuous values per variable) 
#### in each region
### RETURN: data.frame with columns for W variables
get_W = function(data, n_stat_bin, n_stat_cont, n_tv, n_zips) {
  obs_data = data
  n_stat = n_stat_bin + n_stat_cont
  n_var = n_stat_bin + n_stat_cont + n_tv
  bin_label = paste0("X", n_stat_bin) 
  cont_label_1 = paste0("X", n_stat_bin+1)
  cont_label_2 = paste0("X", n_stat)
  tv_label = paste0("X", n_var)
  w_bin_label = paste0("W", n_stat_bin) 
  w_cont_label_1 = paste0("W", n_stat_bin+1)
  w_cont_label_2 = paste0("W", n_stat)
  w_tv_label = paste0("W", n_var)
  
  # W1* = (0.6 + X1 * X3  25)^3
  obs_data[[w_bin_label]] = (0.6 + data[[bin_label]] * data[[cont_label_2]] / 25)^3
  
  # W2* = 10 + X2 / (1 + exp(X3))
  obs_data[[w_cont_label_1]] = 10 + data[[cont_label_1]] / (1 + exp(data[[cont_label_2]]))
  
  # W3* = exp(0.5*X3)
  obs_data[[w_cont_label_2]] = exp(0.5*data[[cont_label_2]])
  
  # W4* = (20 + X2 + X4)^2
  obs_data[[w_tv_label]] = (20+data[[cont_label_1]] + data[[tv_label]])^2
  
  # Only one unique value per zip code: choose closest
  # Note: Not currently functional since cyclic relationship between covariates and treatment
  if (!is.null(n_zips)) {
    cities = unique(obs_data$city)
    for (i in 1:length(cities)) {
      city = cities[i]
      n_zip = n_zips[i]
      W2_k_means = kmeans(x=obs_data[obs_data$city == city, w_cont_label_1], centers=n_zip)
      obs_data[obs_data$city == city, w_cont_label_1] = W2_k_means$centers[W2_k_means$cluster]
      W3_k_means = kmeans(x=obs_data[obs_data$city == city, w_cont_label_2], centers=n_zip)
      obs_data[obs_data$city == city, w_cont_label_2] = W3_k_means$centers[W3_k_means$cluster]
    }
  }
  
  # Wj* = (Wj* - mean(Wj*))/sd(Wj*)
  obs_data[[w_bin_label]] = scale(obs_data[[w_bin_label]])
  #obs_data[[w_bin_label]] = as.integer(obs_data[[w_bin_label]] > 0) # binarize W1
  obs_data[[w_cont_label_1]] = scale(obs_data[[w_cont_label_1]])
  obs_data[[w_cont_label_2]] = scale(obs_data[[w_cont_label_2]])
  obs_data[[w_tv_label]] = scale(obs_data[[w_tv_label]])
  return(obs_data)
}

### get_simulated_data: calls simulation function and compiles data into dataframe
### PARAM @n: number of units in study
### PARAM @n_zips: number of zip codes in each region (length: n_trt)
### PARAM @n_time: number of periods before tax (= number of periods after)
### PARAM @n_stat_bin: number of static binary variables to simulate
### PARAM @n_stat_bin: number of static continuous variables to simulate
### PARAM @n_stat_bin: number of time-varying continuous variables to simulate
### PARAM @beta_0_a: intercept terms for logistic reg ps model (length: 3, 3=n_trt)
### PARAM @beta_a: coefficients of X/W for logistic reg ps model (dim: n_var x 3, 3=n_trt)
### PARAM @pre_lambdas: coefficients of X/W in linear outcome model in pre-tax period (dim: n_var x n_time)
### PARAM @post_lambdas: coefficients of X/W in linear outcome model in post-tax period (dim: n_var x n_time)
### PARAM @alpha_a: Unobserved, time-invariant effects specific to treatment group (length: 3=n_trt)
### PARAM @gamma_p: period-specific time effects (length: n_time)
### PARAM @gamma_y: time-specific effect that changes pre and post tax (length: 1)
### PARAM @tau_zp: (tau(1,)=rep(0,nperiods) for Balt, tau(2,) = period-specific ATN, 
#### tau(3,) = period-specific ATT, tau(4,)=rep(0,nperiods) for non-border,(dim=n_trt x nperiods)
### PARAM @D_vars: vector of variables that affect distance 
### PARAM @D_var_signs: vector of sign for direction which variables affect distance (positive if associated with higher distance)
### PARAM @simulate_A_from_W: if T, simulate treatment from observed data, F is misspecified
### PARAM @simulate_Y_from_W: if T, simulate outcome from observed data, F is misspecified
### RETURN: list of the simulated data frame as well as the true effects
get_simulated_data = function(n, n_zips, n_time, n_stat_bin, n_stat_cont, n_tv, beta_0_a, beta_a, 
                              pre_lambdas, post_lambdas, alpha_a, gamma_p, gamma_y, tau_zp,
                              D_vars, D_var_signs, simulate_A_from_W, simulate_Y_from_W, n_att) {
  n_var = n_stat_bin + n_stat_cont + n_tv
  X = simulate_X(n=n, n_time=n_time, n_stat_bin=n_stat_bin, n_stat_cont=n_stat_cont, n_tv=n_tv)
  W = get_W(data=X, n_stat_bin=n_stat_bin, n_stat_cont=n_stat_cont, n_tv=n_tv, n_zips=n_zips)
  data = merge(X, W[,!(colnames(W) %in% paste0("X",1:n_var))], by=c("ID", "period"))
  A = simulate_A(n=n, data=data, beta_0_a=beta_0_a, beta_a=beta_a, n_var=n_var, simulate_from_W=simulate_A_from_W, n_att=n_att)
  data = merge(data, A, by="ID")
  D = simulate_D(n=n, data=data, vars=D_vars, var_signs=D_var_signs, simulate_from_W=simulate_A_from_W)
  data = merge(data, D, by="ID")
  effective_tau_zp = rbind(tau_zp[1,], tau_zp[2,] / mean(1-data[data$city == 2, "D"]), tau_zp[3,], tau_zp[4,])
  data = simulate_Y(n=n, data=data, n_time=n_time, n_var=n_var, pre_lambdas=pre_lambdas, post_lambdas=post_lambdas, 
                    alpha_a=alpha_a, gamma_p=gamma_p, gamma_y=gamma_y, tau_zp=effective_tau_zp, simulate_from_W=simulate_Y_from_W)
  true_effects = get_true_effects(data=data, tau_zp=tau_zp)
  return(list(data=data, effects=true_effects))
}



