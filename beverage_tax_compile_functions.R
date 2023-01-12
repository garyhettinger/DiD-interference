source("beverage_tax_estimate_functions.R")
library(bayesboot)

collect_estimates_df = function(data, ps_vars, or_vars, manual=T, cluster=NA, weights=NULL) {
  if (manual) {
    twfe = get_manual_twfe_estimate(data, or_vars[or_vars != "PreTaxTarget"], use_interactions=F, weights=weights)
    twfe_cond = get_manual_twfe_estimate(data, or_vars[or_vars != "PreTaxTarget"], use_interactions=T, weights=weights)
    outcome = get_manual_outcome_estimate(data, or_vars, weights=weights)
    ipw = get_manual_ipw_estimate(data, ps_vars, stabilize=F, weights=weights)
    ipw_std = get_manual_ipw_estimate(data, ps_vars, stabilize=T, weights=weights)
    dr = get_manual_dr_estimate(data, ps_vars, or_vars, weights=weights)
    if (!is.na(cluster)) dr = get_manual_dr_estimate_cluster(data, cluster, ps_vars, or_vars, weights=weights)
  }
  else {
    twfe = get_drdid_twfe_estimate(data, or_vars[or_vars != "PreTaxTarget"], use_interactions=F, weights=weights)
    twfe_cond = get_drdid_twfe_estimate(data, or_vars[or_vars != "PreTaxTarget"], use_interactions=T, weights=weights)
    outcome = get_drdid_outcome_estimate(data, or_vars, weights=weights)
    ipw = get_drdid_ipw_estimate(data, ps_vars, stabilize=F, weights=weights)
    ipw_std = get_drdid_ipw_estimate(data, ps_vars, stabilize=T, weights=weights)
    dr = get_drdid_dr_estimate(data, unique(c(ps_vars, or_vars)), weights=weights)
  }
  
  effects = data.frame(twfe = twfe$effect, twfe_cond = twfe_cond$effect,
                       outcome = outcome$effect, 
                       ipw=ipw$effect, ipw_std=ipw_std$effect, 
                       dr = dr$effect)
  variances = data.frame(twfe = twfe$var, twfe_cond = twfe_cond$var,
                         outcome = outcome$var, 
                         ipw=ipw$var, ipw_std=ipw_std$var, 
                         dr = dr$var)
  return(list(effects=effects, variances=variances))
}

transform_to_ratio_estimate = function(data, additive_effects) {
  ratio_effects = additive_effects
  y11 = mean(data[data$trt == 1, "PostTaxTarget"])
  for (name in names(additive_effects)) {
    ratio_effects[[name]] = y11 / (y11 - additive_effects[[name]])
  }
  return(ratio_effects)
}

aggregate_time_estimates = function(period_estims, periods) {
  effects = variances = 0
  for (period in periods) {
    period_label = paste0("period_", period)
    effects = effects + period_estims[[period_label]]$effects / length(periods)
    variances = variances + period_estims[[period_label]]$variances / (length(periods))^2
  }
  return(list(effects=effects, variances=variances))
}

get_parametric_lowers_uppers = function(estims) {
  lowers = estims$effects - 1.96 * sqrt(estims$variances)
  uppers = estims$effects + 1.96 * sqrt(estims$variances)
  return(list(lowers=lowers, uppers=uppers))
}

get_period_specific_estimates = function(full_data, ps_vars, or_vars, manual=T, cluster=NA, weights=NULL) {
  # Fit full twfe model across periods
  period_estims = list()
  periods = unique(full_data$period)
  for (period in periods) {
    period_data = full_data[full_data$period == period,]
    period_weights = weights[full_data$period == period]
    period_label = paste0("period_", period)
    period_estims[[period_label]] = collect_estimates_df(period_data, ps_vars, or_vars, manual=manual, cluster=cluster, weights=period_weights)
  }
  return(period_estims)
}

get_period_specific_ratio_estimates = function(full_data, period_estims) {
  # Fit full twfe model across periods
  period_ratio_estims = list()
  periods = unique(full_data$period)
  for (period in periods) {
    period_data = full_data[full_data$period == period,]
    period_label = paste0("period_", period)
    period_ratio_estims[[period_label]] = list()
    period_ratio_estims[[period_label]]$effects = transform_to_ratio_estimate(period_data, period_estims[[period_label]]$effects)
  }
  return(period_ratio_estims)
}

aggregate_to_season = function(estims) {
  season_period_map = list(Winter=1:3, Spring=4:6, Summer=7:9, Fall=10:13)
  season_estims = list()
  for (season in names(season_period_map)) {
    season_label = paste0("season_", season)
    periods_for_season = season_period_map[[season]]
    estims[[season_label]] = aggregate_time_estimates(estims, periods_for_season)
  }
  return(estims)
}

aggregate_to_year = function(estims, full_data, ps_vars, or_vars, ratio=F, get_sens=F, weights=NULL) {
  estims$year = aggregate_time_estimates(estims, unique(full_data$period))
  if (get_sens) {
    estims$year$effects$outcome_sens = get_manual_outcome_estimate_sens(full_data, or_vars, weights=weights)$effect
    estims$year$effects$ipw_sens = get_manual_ipw_estimate(full_data, ps_vars, weights=weights)$effect
    estims$year$effects$outcome_sens2 = get_manual_outcome_estimate_sens2(full_data, or_vars, weights=weights)$effect
    estims$year$effects$outcome_sens3 = get_manual_outcome_estimate_sens3(full_data, or_vars[or_vars != "PreTaxTarget"], weights=weights)$effect
    estims$year$effects$dr_sens = get_manual_dr_estimate_sens(full_data, ps_vars, or_vars, weights=weights)$effect
    estims$year$effects$dr_sens2 = get_manual_dr_estimate_sens2(full_data, ps_vars, or_vars, weights=weights)$effect
    estims$year$effects$dr_sens3 = get_manual_dr_estimate_sens3(full_data, ps_vars, or_vars[or_vars != "PreTaxTarget"], weights=weights)$effect
    if (ratio) estims$year$effects[,c("outcome_sens", "ipw_sens", "dr_sens", "outcome_sens2", "dr_sens2","outcome_sens3", "dr_sens3")] = transform_to_ratio_estimate(full_data, 
                     estims$year$effects[,c("outcome_sens", "ipw_sens", "dr_sens", "outcome_sens2", "dr_sens2", "outcome_sens3", "dr_sens3")])
  }
  return(estims)
}

get_estimates = function(data, ps_vars, or_vars, get_season=T, manual=T, ratio=F, get_sens=F, cluster=NA, weights=NULL) {
  estims = get_period_specific_estimates(data, ps_vars, or_vars, manual=manual, weights=weights, cluster=cluster)
  if (ratio) estims = get_period_specific_ratio_estimates(data, estims)
  if (get_season) estims = aggregate_to_season(estims)
  estims = aggregate_to_year(estims, data, ps_vars, or_vars, ratio, get_sens, weights=weights)
  return(estims)
}

run_parametric_estimates = function(data, ps_vars, or_vars, get_season, manual, ratio, get_sens, trt_cols, control_cols) {
  full_data = set_treatment_categories_binomial(data=data, trt_cols=trt_cols, control_cols=control_cols)
  estims = get_estimates(full_data, ps_vars, or_vars, get_season=get_season, manual=manual, ratio=ratio, get_sens=get_sens)
  return(estims)
}


get_freq_bootstrap_ind = function(data, ps_vars, or_vars, get_season, manual, ratio, trt_cols, control_cols, get_sens, cluster=NA) {
  data = if (!is.na(cluster)) get_bootstrap_sample_cluster(data, unique(c(ps_vars, or_vars))) else get_bootstrap_sample(data, unique(c(ps_vars, or_vars)))
  full_data = set_treatment_categories_binomial(data=data, trt_cols=trt_cols, control_cols=control_cols)
  estims = get_estimates(full_data, ps_vars, or_vars, get_season=get_season, manual=manual, ratio=ratio, get_sens=get_sens, cluster=cluster)
  return(estims)
}

run_freq_bootstrap = function(data, ps_vars, or_vars, nboots, get_season, manual, ratio, get_sens, trt_cols, control_cols, cluster=NA, truth=NA, ncores=5) {
  boot_estims = parallel::mclapply(1:nboots, function(x) 
    get_freq_bootstrap_ind(data, ps_vars, or_vars, get_season, manual, ratio, trt_cols, control_cols, get_sens, cluster), mc.cores=ncores)
  
  # Untangling return from mclapply to have one row per bootstrap
  boot_stack = stack(unlist(boot_estims))
  boot_df = data.frame(sim=1:nboots)
  for (col in unique(boot_stack$ind)) {
    boot_df[[col]] = boot_stack[boot_stack$ind == col,"values"]
  }
  
  # Aggregate bootstrap samples 
  boot_results = sapply(colnames(boot_df), function(x) c(mean(boot_df[[x]], na.rm=T), var(boot_df[[x]], na.rm=T),
                                                         quantile(boot_df[[x]], probs=c(0.025, 0.975), na.rm=T)))
  
  # Store bootstrap effects, variances, and CIs
  results = list()
  for (col in colnames(boot_results[,2:ncol(boot_results)])) {
    subcols = strsplit(col, split=".", fixed=T)[[1]]
    if (subcols[2] != "effects") next
    if (!(subcols[1] %in% names(results))) results[[subcols[1]]] = list(effects=data.frame(truth=truth), variances=data.frame(truth=truth), 
                                                                        lowers=data.frame(truth=truth), uppers=data.frame(truth=truth))
    results[[subcols[1]]]$effects[[subcols[3]]] = boot_results[1,col]
    results[[subcols[1]]]$variances[[subcols[3]]] = boot_results[2,col]
    results[[subcols[1]]]$lowers[[subcols[3]]] = boot_results[3,col]
    results[[subcols[1]]]$uppers[[subcols[3]]] = boot_results[4,col]
  }
  return(results)
}

get_bb_data = function(bb_ids, orig_data) {
  bb_ids$BBID = 1:nrow(bb_ids)
  data = merge(bb_ids, orig_data, by="ID")
  data$ID = data$BBID
  return(data)
}

run_bayes_bootstrap = function(data, ps_vars, or_vars, nboots, get_season, manual, ratio, get_sens, trt_cols, control_cols, truth=NA) {
  full_data = set_treatment_categories_binomial(data=data, trt_cols=trt_cols, control_cols=control_cols)
  bb_func_call = function (ids) unlist(get_estimates(get_bb_data(ids, full_data), ps_vars, or_vars, get_season=get_season, manual=manual, ratio=ratio, get_sens=get_sens))
  boot_df = bayesboot(data=data.frame(ID=unique(full_data$ID)), statistic=bb_func_call,
                         R=nboots, R2=4000, use.weights=F)
  # Aggregate bootstrap samples 
  boot_results = sapply(colnames(boot_df), function(x) c(mean(boot_df[[x]], na.rm=T), var(boot_df[[x]], na.rm=T),
                                                         quantile(boot_df[[x]], probs=c(0.025, 0.975), na.rm=T)))
  # Store bootstrap effects, variances, and CIs
  results = list()
  for (col in colnames(boot_results[,2:ncol(boot_results)])) {
    subcols = strsplit(col, split=".", fixed=T)[[1]]
    if (subcols[2] != "effects") next
    if (!(subcols[1] %in% names(results))) results[[subcols[1]]] = list(effects=data.frame(truth=truth), variances=data.frame(truth=truth), 
                                                                        lowers=data.frame(truth=truth), uppers=data.frame(truth=truth))
    results[[subcols[1]]]$effects[[subcols[3]]] = boot_results[1,col]
    results[[subcols[1]]]$variances[[subcols[3]]] = boot_results[2,col]
    results[[subcols[1]]]$lowers[[subcols[3]]] = boot_results[3,col]
    results[[subcols[1]]]$uppers[[subcols[3]]] = boot_results[4,col]
  }
  return(results)
}

get_bb_data_wt = function(ids, bb_wts, orig_data) {
  bb_ids = data.frame(ID=ids$ID, weights=bb_wts)
  orig_data$row = 1:nrow(orig_data)
  merged_data = merge(bb_ids, orig_data, by="ID")
  weights = merged_data[order(merged_data$row),"weights"]
  return(weights)
}

run_bayes_bootstrap_wt = function(data, ps_vars, or_vars, nboots, get_season, manual, ratio, get_sens, trt_cols, control_cols, truth=NA) {
  full_data = set_treatment_categories_binomial(data=data, trt_cols=trt_cols, control_cols=control_cols)
  
  bb_func_call = function (ids, weights) unlist(get_estimates(full_data, ps_vars, or_vars, get_season=get_season, manual=manual, ratio=ratio, get_sens=get_sens, weights=get_bb_data_wt(ids, weights, full_data)))
  boot_df = bayesboot(data=data.frame(ID=unique(full_data$ID)), statistic=bb_func_call,
                      R=nboots, use.weights=T)
  # Aggregate bootstrap samples 
  boot_results = sapply(colnames(boot_df), function(x) c(mean(boot_df[[x]], na.rm=T), var(boot_df[[x]], na.rm=T),
                                                         quantile(boot_df[[x]], probs=c(0.025, 0.975), na.rm=T)))
  # Store bootstrap effects, variances, and CIs
  results = list()
  for (col in colnames(boot_results[,2:ncol(boot_results)])) {
    subcols = strsplit(col, split=".", fixed=T)[[1]]
    if (subcols[2] != "effects") next
    if (!(subcols[1] %in% names(results))) results[[subcols[1]]] = list(effects=data.frame(truth=truth), variances=data.frame(truth=truth), 
                                                                        lowers=data.frame(truth=truth), uppers=data.frame(truth=truth))
    results[[subcols[1]]]$effects[[subcols[3]]] = boot_results[1,col]
    results[[subcols[1]]]$variances[[subcols[3]]] = boot_results[2,col]
    results[[subcols[1]]]$lowers[[subcols[3]]] = boot_results[3,col]
    results[[subcols[1]]]$uppers[[subcols[3]]] = boot_results[4,col]
  }
  return(results)
}

get_pretrends_data = function(data) {
  data_0 = data[data$period == 1, c("ID", "PreTaxTarget")]
  data_1 = data[data$period != 1, colnames(data) != "PostTaxTarget"]
  names(data_1) = ifelse(names(data_1) == "PreTaxTarget", "PostTaxTarget", names(data_1))
  pretrends_data = merge(data_0, data_1, by=c("ID"))
  pretrends_data$period = pretrends_data$period - 1
  return(pretrends_data)
}

run_parametric_pretrends = function(data, ps_vars, or_vars, trt_cols, control_cols) {
  pretrends_data = get_pretrends_data(data)
  estims = run_parametric_estimates(pretrends_data, ps_vars, or_vars, get_season=F, manual=F, ratio=F, get_sens=F, trt_cols=trt_cols, control_cols=control_cols)
  return(estims)
}

run_freq_bootstrap_pretrends = function(data, ps_vars, or_vars, nboots, trt_cols, control_cols, truth=NA, ncores=5) {
  pretrends_data = get_pretrends_data(data)
  estims = run_freq_bootstrap(pretrends_data, ps_vars, or_vars, nboots, get_season=F, manual=F, ratio=F, get_sens=F, trt_cols=trt_cols, control_cols=control_cols)
  return(estims)
}

run_bayes_bootstrap_wt_pretrends = function(data, ps_vars, or_vars, nboots, trt_cols, control_cols, truth=NA) {
  pretrends_data = get_pretrends_data(data)
  estims = run_bayes_bootstrap_wt(data, ps_vars, or_vars, nboots, get_season=F, manual=F, ratio=F, get_sens=F, trt_cols=trt_cols, control_cols=control_cols)
  return(estims)
}
