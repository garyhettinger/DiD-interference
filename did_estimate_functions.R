################################################################################
# This file implements the models to estimate ATT and ATN with different DiD estimators
################################################################################


### get_bootstrap_sample: uses a stratified bootstrap sampling approach to do a bootstrap within each treatment group 
### PARAM @data: data.frame with city variable to do stratified bootstrap across
### RETURN: stratified bootstrap sample data.frame of input data.frame
get_bootstrap_sample = function(data) {
  new_data = data.frame()
  new_ID_counter = max(data$ID) + 1 # so can make repeated IDs distinct
  for (city in unique(data$city)) {
    trt_IDs = unique(data[data$city == city, "ID"]) # get IDs within city
    sample_ids = sample(trt_IDs, length(trt_IDs), replace=T) # sample among those city IDs
    for (id in sample_ids) {
      sub_data = data[data$ID == id,]
      sub_data$ID = new_ID_counter
      # Add the data from each ID for each time it is included in bootstrap sample
      new_data = rbind(new_data, sub_data) 
      new_ID_counter = new_ID_counter + 1
    }
  }
  return(new_data)
}

### set_treatment_categories_binomial: takes sub portion of data frame according to defined treatment/control groups
#### and labels trt as 0/1 accordingly
### PARAM @data: data.frame with city variable to use for treatments
### PARAM @trt_cols: city values that are considered in treatment group
### PARAM @control_cols: city values considered as control region
### RETURN: data.frame with labels for trt and only stores from cities in treated/control groups
set_treatment_categories_binomial = function(data, trt_col_name, trts, ctls) {
  sub_data = data[data[[trt_col_name]] %in% c(trts, ctls),]
  sub_data$trt = ifelse(sub_data[[trt_col_name]] %in% trts, rep(1, nrow(sub_data)), rep(0, nrow(sub_data)))
  return(sub_data)
}

### get_stacked_data: takes data with one row per m-time with columns for both Pre and Post Tax Target
### and creates a stacked df with separate rows for Pre and Post Tax outcomes
### PARAM @data: data.frame with PreTaxTarget and PostTaxTarget columns
### RETURN: data.frame with double the rows and a single Target column
get_stacked_data = function(data) {
  data_0 = data[,colnames(data) != "PostTaxTarget"]
  colnames(data_0) = ifelse(colnames(data_0) == "PreTaxTarget", rep("Target", ncol(data_0)), colnames(data_0))
  data_0$Year = 0
  data_1 = data[,colnames(data) != "PreTaxTarget"]
  colnames(data_1) = ifelse(colnames(data_1) == "PostTaxTarget", rep("Target", ncol(data_1)), colnames(data_1))
  data_1$Year = 1
  return(rbind(data_0, data_1))
}

### unstack_data: reverses get_stacked data by creating a dataframe with separate Pre and PostTax outcome columns
### PARAM @stacked_data: data.frame with separate rows for Pre and PostTax outcomes
### RETURN: data.frame with half the rows and a single Target (and Fitted Target if applicable) column
unstack_data = function(stacked_data) {
  data_0 = stacked_data[1:(nrow(stacked_data)/2),]
  colnames(data_0) = ifelse(colnames(data_0) %in% c("Target", "FittedTarget"), sapply(colnames(data_0), function(x) paste0("PreTax",x)),
                            colnames(data_0))
  
  data_1 = stacked_data[(nrow(stacked_data)/2 + 1):nrow(stacked_data),c("Target", "FittedTarget")]
  colnames(data_1) = ifelse(colnames(data_1) %in% c("Target", "FittedTarget"), sapply(colnames(data_1), function(x) paste0("PostTax",x)),
                            colnames(data_1))
  return(cbind(data_0, data_1))
}


### fit_ps_model_binomial: fits a logistic regression propensity score model to data and adds columns for trt probabilities
### PARAM @data: data.frame with columns for (ID, trt, vars, PreTaxTarget, PostTaxTarget)
### PARAM @variables: variables in PS model
### PARAM @ia_terms: vector of interaction terms for PS model
### RETURN: data.frame with added columns for Control_Prob and Treat_Prob
fit_ps_model_binomial = function(data, variables, ia_terms=c()) {
  form = as.formula(paste0("trt ~ ", paste0(c(variables, ia_terms), collapse = " + "))) # no interactions
  ps_data = aggregate(cbind(data[,variables]), by=data[,c("ID","trt")], FUN=mean)
  model = glm(form, family=binomial, data=ps_data)
  probs = data.frame(ID=ps_data$ID, Control_Prob = 1-model$fitted.values, Treat_Prob = model$fitted.values)
  data = merge(data, probs, by="ID")
  return(data)
}

### fit_outcome_model: fits outcome model on difference in pre and post tax target
### PARAM @data: data.frame with columns for (ID, trt, vars, PreTaxTarget, PostTaxTarget)
### PARAM @variables: confounder variables in model
### PARAM @ia_terms: vector of interaction terms for OR model
### RETURN: data.frame with FittedTarget outcome variable
fit_outcome_model = function(data, variables, ia_terms=c()) {
  # Fit on control data
  form = as.formula(paste0("Target ~ ", paste0(c(variables, ia_terms), collapse=" + ")))
  outcome_model_data = data[data$trt == 0,]
  outcome_model_data$Target = outcome_model_data$PostTaxTarget - outcome_model_data$PreTaxTarget
  model = lm(form, data=outcome_model_data)
  
  # Estimate on full data
  data$FittedTarget = predict(model, newdata=data)
  return(data)
}


### get_twfe_estimate: fits twfe model on outcome and interprets coefficient as ATT/ATN
### PARAM @data: unstacked data.frame with columns for (ID, trt, vars, PreTaxTarget, PostTaxTarget)
### PARAM @variables: confounder variables in model
### PARAM @use_period_interactions: boolean if want to address time-varying covariate effects
### RETURN: point estimate for effect
get_twfe_estimate = function(data, variables, use_period_interactions=F) {
  # Fit outcome model with additional Year and trt based variables using full data
  stacked_data=get_stacked_data(data)
  stacked_weights = if (is.null(weights)) NULL else c(weights, weights)
  twfe_vars = c("Year", "trt", variables)
  if (use_period_interactions) {
    twfe_vars = c(twfe_vars, paste0("Year:", variables))
  }
  twfe_vars = c(twfe_vars, "Year:trt")
  form = as.formula(paste0("Target ~ ", paste0(twfe_vars, collapse=" + "), " + (1 | ID)"))
  model = lmer(form, data=stacked_data, weights=stacked_weights)
  effect = summary(model)$coefficients["Year:trt", "Estimate"]
  return(effect)
}

### get_or_estimate: implements Heckman outcome estimator based on an outcome model 
### PARAM @data: data.frame with columns for (ID, trt, vars, PreTaxTarget, PostTaxTarget)
### PARAM @vars: variables to include in model
### PARAM @ia_terms: additional interaction terms to include in OR model
### RETURN: point estimate for effect
get_or_estimate = function(data, vars, ia_terms=c()) {
  data = fit_outcome_model(data, vars, ia_terms=ia_terms)
  trt_data = data[data$trt == 1,]
  trt_diff = trt_data[,"PostTaxTarget"] - trt_data[,"PreTaxTarget"]
  control_diff = trt_data[,"FittedTarget"]
  effect = mean(trt_diff - control_diff)
  return(effect)
}

### get_ipw_estimate: implements Abadie outcome estimator based on fitted propensity scores
### PARAM @data: data.frame with columns for (ID, trt, vars, PreTaxTarget, PostTaxTarget)
### PARAM @vars: vector of column names for variables in PS model
### PARAM @ia_terms: vector of interaction terms to include in PS model
### RETURN: point estimate of effect
get_ipw_estimate = function(data, vars, ia_terms=c()) {
  data = fit_ps_model_binomial(data, vars, ia_terms=ia_terms)
  treat.A = as.integer(data$trt == 1)
  treat.A.hat = data$Treat_Prob
  control.A = as.integer(data$trt == 0)
  control.A.hat = data$Control_Prob
  w1 = treat.A / mean(treat.A)
  w0 = (treat.A.hat * control.A / control.A.hat) / mean(treat.A)
  effect = mean((w1-w0) * (data$PostTaxTarget - data$PreTaxTarget))
  return(effect)
}

### get_dr_estimate: implements Santanna & Zhao/Li & Li DR estimator based on fitted propensity scores and outcome estimates
### PARAM @data: data.frame with columns for (ID, trt, vars, PreTaxTarget, PostTaxTarget)
### PARAM @ps_vars: vector of variables for PS model
### PARAM @or_vars: vector of variables for OR model
### PARAM @ps_ia_terms: vector of interaction terms for PS model
### PARAM @or_ia_terms: vector of interaction terms for OR model
### RETURN: point estimate of effect
get_dr_estimate = function(data, ps_vars, or_vars, ps_ia_terms=c(), or_ia_terms=c()) {
  data = fit_outcome_model(data, or_vars, ia_terms=or_ia_terms)
  data = fit_ps_model_binomial(data, ps_vars, ia_terms=ps_ia_terms)
  treat.A = as.integer(data$trt == 1)
  treat.A.hat = data$Treat_Prob
  control.A = as.integer(data$trt == 0)
  control.A.hat = data$Control_Prob
  w1 = treat.A / mean(treat.A)
  w0 = (treat.A.hat * control.A / control.A.hat) / mean(treat.A)
  delta.Y = data$PostTaxTarget - data$PreTaxTarget
  delta.mu = data$FittedTarget
  effect = mean((w1-w0)*(delta.Y-delta.mu))
  return(effect)
}

get_estimates_single = function(full_data, trt_col_name, trts, ctls, period_col=c("period"), method="dr", ps_vars=NULL, 
                         or_vars=NULL, ps_ia_terms=c(), or_ia_terms=c()) {
  data = set_treatment_categories_binomial(full_data, trt_col_name, trts=trts, ctls=ctls)
  periods = sort(unique(data[[period_col]]))
  effs = c()
  for (p in periods) {
    p_data = data[data[[period_col]] == p,]
    if (method == "dr") eff = get_dr_estimate(p_data, ps_vars, or_vars, ps_ia_terms, or_ia_terms)
    else if (method == "ipw") eff = get_ipw_estimate(p_data, vars, ps_ia_terms)
    else if (method == "or") eff = get_or_estimate(p_data, vars, or_ia_terms)
    else if (method == "twfe_cond") eff = get_twfe_estimate(p_data, variables, use_period_interactions=T)
    else eff = get_twfe_estimate(p_data, variables, use_period_interactions=F)
    effs = c(effs, eff)
  }
  df = data.frame(period=periods, effect=effs)
  df = rbind(df, list(period="FULL", effect=mean(effs)))
  return(df)
}

get_bootstrap_estimates_single = function(data, trt_col_name, trts, ctls, nboots=500, 
                                          period_col="period", method="dr", ps_vars=NULL, 
                                          or_vars=NULL, ps_ia_terms=c(), or_ia_terms=c()) {
  boot_df = data.frame()
  for (i in 1:nboots) {
    set.seed(i)
    boot_data = get_bootstrap_sample(data)
    boot_res = get_estimates_single(boot_data, trt_col_name, trts, ctls, period_col, method, 
                             ps_vars, or_vars, ps_ia_terms, or_ia_terms)
    names(boot_res) = c("period", paste0("effect_",i))
    if (nrow(boot_df) == 0) boot_df = boot_res
    else boot_df = merge(boot_df, boot_res, by="period")
  }
  results = data.frame()
  for (p in boot_df$period) {
    estims = unlist(boot_df[boot_df$period == p,2:ncol(boot_df)])
    eff = mean(estims)
    pctiles = stats::quantile(estims, c(0.025, 0.975))
    res = list(period=p, effect=eff, lower=pctiles[1], upper=pctiles[2])
    results = rbind(results, res)
  }
  return(results)
}

get_estimates = function(data, trt_col_name, trts, neighbors, trt_ctls, neighbor_ctls,
                         nboots=500, period_col="period", method="dr", ps_vars=NULL, 
                         or_vars=NULL, ps_ia_terms=c(), or_ia_terms=c()) {
  att_estims = get_bootstrap_estimates_single(data, trt_col_name, trts, trt_ctls, nboots, 
                                              period_col, method, ps_vars, 
                                              or_vars, ps_ia_terms, or_ia_terms)
  atn_estims = get_bootstrap_estimates_single(data, trt_col_name, neighbors, neighbor_ctls, nboots, 
                                              period_col, method, ps_vars, 
                                              or_vars, ps_ia_terms, or_ia_terms)
  return(list(att=att_estims, atn=atn_estims))
}
