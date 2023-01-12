library(DRDID)
library(lme4)
library(stats)

################################################################################
# This file implements the models to estimate ATT and ATN with different DiD estimators
################################################################################

weighted_mean = function(x, w) {
  if (is.null(w)) return(mean(x))
  return(weighted.mean(x,w))
}

### get_bootstrap_sample: uses a stratified bootstrap sampling approach to do a bootstrap within each treatment group 
### PARAM @data: data.frame with city variable to do stratified bootstrap across
### RETURN: stratified bootstrap sample data.frame of input data.frame
get_bootstrap_sample = function(data, vars) {
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
  for (var in vars) {
    for (city in unique(new_data$city)) {
      if (length(unique(new_data[new_data$city == city, var])) == 1) new_data = get_bootstrap_sample(data, vars)
    }
  }
  return(new_data)
}

get_bootstrap_sample_cluster = function(data, vars) {
  new_data = data.frame()
  new_ID_counter = max(data$ID) + 1 # so can make repeated IDs distinct
  for (cluster in unique(data$Cluster)) {
    trt_IDs = unique(data[data$Cluster == cluster, "ID"]) # get IDs within city
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
set_treatment_categories_binomial = function(data, trt_cols, control_cols) {
  sub_data = data[data$city %in% c(trt_cols, control_cols),]
  sub_data$trt = ifelse(sub_data$city %in% trt_cols, rep(1, nrow(sub_data)), rep(0, nrow(sub_data)))
  return(sub_data)
}

get_stacked_data = function(data) {
  data_0 = data[,colnames(data) != "PostTaxTarget"]
  colnames(data_0) = ifelse(colnames(data_0) == "PreTaxTarget", rep("Target", ncol(data_0)), colnames(data_0))
  data_0$Year = 0
  data_1 = data[,colnames(data) != "PreTaxTarget"]
  colnames(data_1) = ifelse(colnames(data_1) == "PostTaxTarget", rep("Target", ncol(data_1)), colnames(data_1))
  data_1$Year = 1
  return(rbind(data_0, data_1))
}

unstack_data = function(stacked_data) {
  data_0 = stacked_data[1:(nrow(stacked_data)/2),]
  colnames(data_0) = ifelse(colnames(data_0) %in% c("Target", "FittedTarget"), sapply(colnames(data_0), function(x) paste0("PreTax",x)),
                           colnames(data_0))
  
  data_1 = stacked_data[(nrow(stacked_data)/2 + 1):nrow(stacked_data),c("Target", "FittedTarget")]
  colnames(data_1) = ifelse(colnames(data_1) %in% c("Target", "FittedTarget"), sapply(colnames(data_1), function(x) paste0("PostTax",x)),
                            colnames(data_1))
  return(cbind(data_0, data_1))
}


### fit_ps_model_binomial: fits a propensity score model to data and adds columns for trt probabilities
### PARAM @data: data.frame with columns for (ID, trt, (variables))
### PARAM @variables: variables in PS model
### PARAM @super: boolean that is T if use superlearner PS model (currently elastic net) compared to multinomial logistic reg
### RETURN: data.frame with added columns for Control_Prob and Treat_Prob
fit_ps_model_binomial = function(data, variables, weights) {
  form = as.formula(paste0("trt ~ ", paste0(variables, collapse = " + "))) # no interactions
  ps_data = aggregate(cbind(data[,variables]), by=data[,c("ID","trt")], FUN=mean)
  if (!is.null(weights)) weights = aggregate(weights, by=data[,c("ID","trt")], FUN=mean)$weights
  model = glm(form, family=binomial, data=ps_data, weights=weights)
  probs = data.frame(ID=ps_data$ID, Control_Prob = 1-model$fitted.values, Treat_Prob = model$fitted.values)
  data = merge(data, probs, by="ID")
  return(data)
}

### fit_outcome_model: fits outcome model on difference in pre and post tax target
### PARAM @data: data.frame with columns for (PreTaxTarget, PostTaxTarget)
### PARAM @variables: confounder variables in model
### RETURN: data.frame with FittedTarget outcome variable
fit_outcome_model = function(data, variables, weights) {
  # Fit on control data
  form = as.formula(paste0("Target ~ ", paste0(variables, collapse=" + ")))
  outcome_model_data = data[data$trt == 0,]
  outcome_model_data$Target = outcome_model_data$PostTaxTarget - outcome_model_data$PreTaxTarget
  outcome_model_weights = if (is.null(weights)) NULL else weights[data$trt == 0]
  model = lm(form, data=outcome_model_data, weights=outcome_model_weights)
  
  # Estimate on full data
  data$FittedTarget = predict(model, newdata=data)
  return(data)
}

fit_full_outcome_model = function(data, variables, weights) {
  # Fit on control data
  form = as.formula(paste0("Target ~ ", paste0(c(variables, "season"), collapse=" + "), " + (1 | ID)"))
  outcome_model_data = data[data$trt == 0,]
  outcome_model_data$Target = outcome_model_data$PostTaxTarget - outcome_model_data$PreTaxTarget
  outcome_model_weights = if (is.null(weights)) NULL else weights[data$trt == 0]
  model = lmer(form, data=outcome_model_data, weights=outcome_model_weights)
  
  # Estimate on full data
  data$FittedTarget = predict(model, newdata=data, allow.new.levels=T)
  return(data)
}

fit_full_both_outcome_model = function(data, variables, weights) {
  stacked_data = get_stacked_data(data)
  # Fit on entire
  form = as.formula(paste0("Target ~ ", paste0(c(variables, "season"), collapse=" + "), " + ", 
                           paste0(paste0(variables, ":Year"), collapse=" + "), " + Year*trt + (1 | ID)"))
  stacked_weights = if (is.null(weights)) NULL else c(weights, weights)
  model = lmer(form, data=stacked_data, weights=stacked_weights)
  
  # Estimate on full data
  pred_data = stacked_data
  pred_data$trt = 0
  stacked_data$FittedTarget = predict(model, newdata=pred_data, allow.new.levels=T)
  data = unstack_data(stacked_data)
  data$FittedTarget = data$PostTaxFittedTarget - data$PreTaxFittedTarget
  return(data)
}

fit_zip_outcome_model = function(data, variables, weights) {
  # Fit on control data
  form = as.formula(paste0("Target ~ ", paste0(variables, collapse=" + "), " + (1 | zip)"))
  outcome_model_data = data[data$trt == 0,]
  outcome_model_data$Target = outcome_model_data$PostTaxTarget - outcome_model_data$PreTaxTarget
  outcome_model_weights = if (is.null(weights)) NULL else weights[data$trt == 0]
  model = lmer(form, data=outcome_model_data, weights=outcome_model_weights)
  
  # Estimate on full data
  data$FittedTarget = predict(model, newdata=data, allow.new.levels=T)
  return(data)
}

### get_twfe_estimator: fits twfe model on outcome and interprets coefficient as ATT/ATN
### PARAM @stacked_data: data.frame with column for Target
### PARAM @variables: confounder variables in model
### RETURN: list of effect and variance estimates
get_manual_twfe_estimate = function(data, variables, use_interactions=F, weights=NULL) {
  # Fit outcome model with additional Year and trt based variables using full data
  stacked_data=get_stacked_data(data)
  stacked_weights = if (is.null(weights)) NULL else c(weights, weights)
  twfe_vars = c("Year", "trt", variables)
  if (use_interactions) {
    twfe_vars = c(twfe_vars, paste0("Year:", variables))
  }
  twfe_vars = c(twfe_vars, "Year:trt")
  form = as.formula(paste0("Target ~ ", paste0(twfe_vars, collapse=" + "), " + (1 | ID)"))
  model = lmer(form, data=stacked_data, weights=stacked_weights)
  effect = summary(model)$coefficients["Year:trt", "Estimate"]
  var =(summary(model)$coefficients["Year:trt", "Std. Error"])^2
  return(list(effect=effect, var=var))
}

### get_outcome_variance: calculates residual variance and variance of mean
get_outcome_model_projection_variance_unstacked = function(data, vars, trt_diff, control_diff) {
  trt_data = data[data$trt == 1,]
  control_data = data[data$trt == 0,]
  res_sigma_sq = sum(((control_data[,"PostTaxTarget"] - control_data[,"PreTaxTarget"]) - 
                        (control_data[,"FittedTarget"]))^2)/(nrow(control_data)-length(vars)-1)
  X = as.matrix(cbind(data.frame(W0=rep(1, nrow(control_data))), control_data[,vars]))
  X.star = as.matrix(cbind(data.frame(W0=rep(1, nrow(trt_data))), trt_data[,vars]))
  vm = tryCatch({ res_sigma_sq * (X.star %*% solve(t(X) %*% X) %*% t(X.star)) },
                error=function(cond) { return(NA) })
  v1 = sum(vm) / length(trt_diff)^2
  v2 = var(trt_diff - control_diff) / length(trt_diff)
  v = v1 + v2
  return(v)
}

### get_outcome_estimator_unstacked: implements Heckman outcome estimator based on outcome projections
### PARAM @data: data.frame with columns for (trt, PreTaxTarget, PostTaxTarget, FittedTarget)
### RETURN: list of effect and variance estimates
get_manual_outcome_estimate = function(data, vars, weights=NULL) {
  data = fit_outcome_model(data, vars, weights)
  trt_data = data[data$trt == 1,]
  trt_weights = weights[data$trt == 1]
  trt_diff = trt_data[,"PostTaxTarget"] - trt_data[,"PreTaxTarget"]
  control_diff = trt_data[,"FittedTarget"]
  effect = weighted_mean(trt_diff - control_diff, trt_weights)
  var = get_outcome_model_projection_variance_unstacked(data=data, vars=vars, 
                                                        trt_diff=trt_diff, control_diff=control_diff)
  return(list(effect=effect, var=var))
}

get_manual_outcome_estimate_sens = function(data, vars, weights=NULL) {
  data = fit_full_outcome_model(data, vars, weights)
  trt_data = data[data$trt == 1,]
  trt_weights = weights[data$trt == 1]
  trt_diff = trt_data[,"PostTaxTarget"] - trt_data[,"PreTaxTarget"]
  control_diff = trt_data[,"FittedTarget"]
  effect = weighted_mean(trt_diff - control_diff, trt_weights)
  var = get_outcome_model_projection_variance_unstacked(data=data, vars=vars, 
                                                        trt_diff=trt_diff, control_diff=control_diff)
  return(list(effect=effect, var=var))
}

get_manual_outcome_estimate_sens2 = function(data, vars, weights=NULL) {
  data = fit_zip_outcome_model(data, vars, weights)
  trt_data = data[data$trt == 1,]
  trt_weights = weights[data$trt == 1]
  trt_diff = trt_data[,"PostTaxTarget"] - trt_data[,"PreTaxTarget"]
  control_diff = trt_data[,"FittedTarget"]
  effect = weighted_mean(trt_diff - control_diff, trt_weights)
  var = get_outcome_model_projection_variance_unstacked(data=data, vars=vars, 
                                                        trt_diff=trt_diff, control_diff=control_diff)
  return(list(effect=effect, var=var))
}

get_manual_outcome_estimate_sens3 = function(data, vars, weights=NULL) {
  data = fit_full_both_outcome_model(data, vars, weights)
  trt_data = data[data$trt == 1,]
  trt_weights = weights[data$trt == 1]
  trt_diff = trt_data[,"PostTaxTarget"] - trt_data[,"PreTaxTarget"]
  control_diff = trt_data[,"FittedTarget"]
  effect = weighted_mean(trt_diff - control_diff, trt_weights)
  var = get_outcome_model_projection_variance_unstacked(data=data, vars=vars, 
                                                        trt_diff=trt_diff, control_diff=control_diff)
  return(list(effect=effect, var=var))
}

### get_ipw_estimator: implements Abadie outcome estimator based on fitted propensity scores
### PARAM @data: data.frame with columns for (trt, Treat_Prob, Control_Prob)
### PARAM @normalize_weights: boolean T if use Sant'anna denominator for w0, F if use Li and Li
### RETURN: list of effect and variance estimates
get_manual_ipw_estimate = function(data, vars, stabilize=F, weights=NULL) {
  data = fit_ps_model_binomial(data, vars, weights)
  treat.A = as.integer(data$trt == 1)
  treat.A.hat = data$Treat_Prob
  control.A = as.integer(data$trt == 0)
  control.A.hat = data$Control_Prob
  w1 = treat.A / mean(treat.A)
  w0 = (treat.A.hat * control.A / control.A.hat) / mean(treat.A)
  if (stabilize) {
    w1 = w1 / mean(w1)
    w0 = w0 / mean(w0)
  }
  effect = weighted_mean((w1-w0) * (data$PostTaxTarget - data$PreTaxTarget), weights)
  #var = var((w1-w0) * (data$PostTaxTarget - data$PreTaxTarget))/length(treat.A) #- w1*effect)/length(treat.A)
  var = mean(((w1-w0)*(data$PostTaxTarget - data$PreTaxTarget) - effect)^2)
  return(list(effect=effect, var=var))
}

### get_dr_estimator_stacked: implements Santanna/Li DR estimator based on fitted propensity scores and outcome estimates
### PARAM @data: data.frame with columns for (trt, Treat_Prob, Control_Prob, PreTaxTarget, PostTaxTarget, FittedTarget)
### PARAM @normalize_weights: boolean T if use Sant'anna denominator for w0, F if use Li and Li
### RETURN: list of effect and variance estimates
get_manual_dr_estimate = function(data, ps_vars, or_vars, weights=NULL) {
  data = fit_outcome_model(data, or_vars, weights)
  data = fit_ps_model_binomial(data, ps_vars, weights)
  treat.A = as.integer(data$trt == 1)
  treat.A.hat = data$Treat_Prob
  control.A = as.integer(data$trt == 0)
  control.A.hat = data$Control_Prob
  w1 = treat.A / mean(treat.A)
  w0 = (treat.A.hat * control.A / control.A.hat) / mean(treat.A)
  delta.Y = data$PostTaxTarget - data$PreTaxTarget
  delta.mu = data$FittedTarget
  effect = weighted_mean((w1-w0)*(delta.Y-delta.mu), weights)
  #var = var((w1-w0)*(delta.Y-delta.mu)-w1*effect)/length(treat.A)
  var = mean(((w1-w0)*(delta.Y-delta.mu)-effect)^2)/length(treat.A)
  return(list(effect=effect, var=var))
}

# TODO
get_manual_dr_estimate_sens = function(data, ps_vars, or_vars, weights=NULL) {
  data = fit_full_outcome_model(data, or_vars, weights)
  data = fit_ps_model_binomial(data, ps_vars, weights)
  treat.A = as.integer(data$trt == 1)
  treat.A.hat = data$Treat_Prob
  control.A = as.integer(data$trt == 0)
  control.A.hat = data$Control_Prob
  w1 = treat.A / mean(treat.A)
  w0 = (treat.A.hat * control.A / control.A.hat) / mean(treat.A)
  delta.Y = data$PostTaxTarget - data$PreTaxTarget
  delta.mu = data$FittedTarget
  effect = weighted_mean((w1-w0)*(delta.Y-delta.mu), weights)
  #var = var((w1-w0)*(delta.Y-delta.mu)-w1*effect)/length(treat.A)
  var = mean(((w1-w0)*(delta.Y-delta.mu)-effect)^2)/length(treat.A)
  return(list(effect=effect, var=var))
}

get_manual_dr_estimate_sens2 = function(data, ps_vars, or_vars, weights=NULL) {
  data = fit_zip_outcome_model(data, or_vars, weights)
  data = fit_ps_model_binomial(data, ps_vars, weights)
  treat.A = as.integer(data$trt == 1)
  treat.A.hat = data$Treat_Prob
  control.A = as.integer(data$trt == 0)
  control.A.hat = data$Control_Prob
  w1 = treat.A / mean(treat.A)
  w0 = (treat.A.hat * control.A / control.A.hat) / mean(treat.A)
  delta.Y = data$PostTaxTarget - data$PreTaxTarget
  delta.mu = data$FittedTarget
  effect = weighted_mean((w1-w0)*(delta.Y-delta.mu), weights)
  #var = var((w1-w0)*(delta.Y-delta.mu)-w1*effect)/length(treat.A)
  var = mean(((w1-w0)*(delta.Y-delta.mu)-effect)^2)/length(treat.A)
  return(list(effect=effect, var=var))
}

get_manual_dr_estimate_sens3 = function(data, ps_vars, or_vars, weights=NULL) {
  data = fit_full_both_outcome_model(data, or_vars, weights)
  data = fit_ps_model_binomial(data, ps_vars, weights)
  treat.A = as.integer(data$trt == 1)
  treat.A.hat = data$Treat_Prob
  control.A = as.integer(data$trt == 0)
  control.A.hat = data$Control_Prob
  w1 = treat.A / mean(treat.A)
  w0 = (treat.A.hat * control.A / control.A.hat) / mean(treat.A)
  delta.Y = data$PostTaxTarget - data$PreTaxTarget
  delta.mu = data$FittedTarget
  effect = weighted_mean((w1-w0)*(delta.Y-delta.mu), weights)
  #var = var((w1-w0)*(delta.Y-delta.mu)-w1*effect)/length(treat.A)
  var = mean(((w1-w0)*(delta.Y-delta.mu)-effect)^2)/length(treat.A)
  return(list(effect=effect, var=var))
}

get_manual_dr_estimate_cluster = function(data, cluster, ps_vars, or_vars, weights=NULL) {
  data = fit_outcome_model(data, or_vars, weights)
  data = fit_ps_model_binomial(data, ps_vars, weights)
  data = data[(data$Cluster == cluster) | (data$trt == 0),]
  treat.A = as.integer(data$trt == 1)
  treat.A.hat = data$Treat_Prob
  control.A = as.integer(data$trt == 0)
  control.A.hat = data$Control_Prob
  w1 = treat.A / mean(treat.A)
  w0 = (treat.A.hat * control.A / control.A.hat) / mean(treat.A)
  delta.Y = data$PostTaxTarget - data$PreTaxTarget
  delta.mu = data$FittedTarget
  effect = weighted_mean((w1-w0)*(delta.Y-delta.mu), weights)
  #var = var((w1-w0)*(delta.Y-delta.mu)-w1*effect)/length(treat.A)
  var = mean(((w1-w0)*(delta.Y-delta.mu)-effect)^2)/length(treat.A)
  return(list(effect=effect, var=var))
}

get_drdid_twfe_estimate = function(data, vars, use_interactions=F, weights=NULL) {
  y1 = data$PostTaxTarget
  y0 = data$PreTaxTarget
  D = data$trt
  X = data[,vars]
  if (use_interactions) {
    twfe = twfe_did_panel(y1=y1, y0=y0, D=D, covariates=X, boot=F, i.weights=weights)
  }
  else {
    twfe = twfe_did_panel(y1=y1, y0=y0, D=D, covariates=NULL, boot=F, i.weights=weights)
  }
  return(list(effect=twfe$ATT, var=twfe$se^2))
}

get_drdid_outcome_estimate = function(data, vars, weights=NULL) {
  y1 = data$PostTaxTarget
  y0 = data$PreTaxTarget
  D = data$trt
  X = data[,vars]
  outcome = reg_did_panel(y1=y1, y0=y0, D=D, covariates=X, boot=F, i.weights=weights)
  return(list(effect=outcome$ATT, var=outcome$se^2))
}

get_drdid_ipw_estimate = function(data, vars, stabilize=F, weights=NULL) {
  y1 = data$PostTaxTarget
  y0 = data$PreTaxTarget
  D = data$trt
  X = data[,vars]
  if (stabilize) {
    ipw = std_ipw_did_panel(y1=y1, y0=y0, D=D, covariates=X, boot=F, i.weights=weights)
  }
  else {
    ipw = ipw_did_panel(y1=y1, y0=y0, D=D, covariates=X, boot=F, i.weights=weights)
  }
  return(list(effect=ipw$ATT, var=ipw$se^2))
}

get_drdid_dr_estimate = function(data, vars, weights=NULL) {
  y1 = data$PostTaxTarget
  y0 = data$PreTaxTarget
  D = data$trt
  X = data[,vars]
  dr = drdid_panel(y1=y1, y0=y0, D=D, covariates=X, boot=F, i.weights=weights)
  return(list(effect=dr$ATT, var=dr$se^2))
}

