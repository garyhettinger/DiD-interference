library(MASS)
get_ps0_info = function(data, ps_covs, trt_col="A", 
                        is_boot=F, wt_col="boot_wts", avg_ps_X=T, full_period_data=NULL) {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  ps0_form = as.formula(paste0(trt_col, " ~ ", paste0(ps_covs, collapse=" + ")))
  
  if (avg_ps_X) {
    full_period_data$boot_wts = sapply(full_period_data$ID, function(x) data[data$ID == x, "boot_wts"])
    ps_data = full_period_data %>% group_by_at(c("ID", trt_col)) %>% 
      summarise_at(c(ps_covs, "boot_wts"), mean) %>% ungroup %>% data.frame()
    ps0_model = glm(ps0_form, data=ps_data, family=binomial(link="logit"), weights=boot_wts)
    ps_map = data.frame(ID=ps_data$ID, ps=predict(ps0_model, newdata=ps_data, type="response"))
    piA_mean = sapply(data$ID, function(x) ps_map[ps_map$ID == x, "ps"])
  }
  else {
    ps0_model = glm(ps0_form, data=data, family=binomial(link="logit"), weights=boot_wts)
    piA_mean = predict(ps0_model, newdata=data, type="response")
  }
  return(list(ps0_form=ps0_form, ps0_model=ps0_model, ps=piA_mean))
}

get_ps0_wts = function(data, ps, trt_col="A", normalize=T, trim=F) {
  wts = ifelse(data[[trt_col]] == 1, 1, ps/(1-ps))
  if (trim) wts = WeightIt::trim(w=wts, at=0.95, lower=T, treat=data[[trt_col]])
  if (normalize) wts[data[[trt_col]]==0] = wts[data[[trt_col]]==0]/mean(wts[data[[trt_col]]==0])
  return(wts)
}

get_or0_info = function(data, or_covs, outcome_col="deltaY", trt_col="A",
                        is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  or0_form = as.formula(paste0(outcome_col, " ~ ", paste0(or_covs, collapse=" + ")))
  or0_model = glm(or0_form, data=data[data[[trt_col]]==0,], family=gaussian(link="identity"), weights=boot_wts)
  mu0 = predict(or0_model, newdata=data, type="response")
  return(list(or0_form=or0_form, or0_model=or0_model, predMu0=mu0))
}

get_theta0_est_twfe = function(data, or_covs, outcome_col0="Y0", outcome_col1="Y1", 
                               trt_col="A", is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  form = as.formula(paste0("Y ~ ", paste0(c(or_covs, "Time", 
                                            trt_col, paste0(trt_col, ":Time"), 
                                            paste0(or_covs, ":Time")), collapse=" + ")))
  stacked_data1 = data %>% dplyr::select(-!!sym(outcome_col0)) %>% rename(Y=outcome_col1) %>% mutate(Time=1)
  stacked_data0 = data %>% dplyr::select(-!!sym(outcome_col1)) %>% rename(Y=outcome_col0) %>% mutate(Time=0)
  stacked_data = rbind(stacked_data1, stacked_data0)
  m1 = lm(form, data=stacked_data, weights=boot_wts)
  est = m1$coefficients[[paste0("Time:", trt_col)]]
  return(est)
}

get_theta0_est_naive = function(data, outcome_col="deltaY", trt_col="A", is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  ctl_data = data[data[[trt_col]]==0,]
  return(stats::weighted.mean(x=ctl_data[[outcome_col]], w=ctl_data$boot_wts))
}

get_theta0_est_or = function(data, predMu0, outcome_col="deltaY", trt_col="A", is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  return(stats::weighted.mean(x=predMu0[data[[trt_col]]==1], w=data$boot_wts[data[[trt_col]]==1]))
}

get_theta0_est_ipw = function(data, wts, outcome_col="deltaY", trt_col="A", is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  return(stats::weighted.mean(x=wts[data[[trt_col]]==0]*data[[outcome_col]][data[[trt_col]]==0], w=data$boot_wts[data[[trt_col]]==0]))
}

get_theta0_est_dr = function(data, wts, predMu0, outcome_col="deltaY", trt_col="A", is_boot=F, wt_col="boot_wts") {
  data$boot_wts = if (!is_boot) rep(1, nrow(data)) else data[[wt_col]]
  return(stats::weighted.mean(x=ifelse(data[[trt_col]]==1, wts*predMu0, wts*(data[[outcome_col]] - predMu0)), w=data$boot_wts)/
           stats::weighted.mean(x=data[[trt_col]]==1, w=data$boot_wts))
}

get_theta0_ests = function(data, or_covs, ps_covs, is_boot=F, trt_col="A", 
                           outcome_col="deltaY", outcome_col0="Y0", outcome_col1="Y1",
                           wt_col="boot_wts", wt_norm=T, avg_ps_X=T, full_period_data=NULL) {
  ps0_info = get_ps0_info(data=data, ps_covs=ps_covs, trt_col=trt_col, is_boot=is_boot, 
                          wt_col=wt_col, avg_ps_X=avg_ps_X, full_period_data=full_period_data)
  or0_info = get_or0_info(data=data, or_covs=or_covs, outcome_col=outcome_col, 
                          trt_col=trt_col, is_boot=is_boot, wt_col=wt_col)
  ps0_wts = get_ps0_wts(data=data, ps=ps0_info$ps, trt_col=trt_col, normalize=wt_norm)
  
  twfe_est = get_theta0_est_twfe(data=data, or_covs=or_covs, outcome_col0=outcome_col0, 
                                 outcome_col1=outcome_col1, trt_col=trt_col, is_boot=is_boot, 
                                 wt_col=wt_col) 
  naive_est = get_theta0_est_naive(data=data, outcome_col=outcome_col, trt_col=trt_col, 
                                   is_boot=is_boot, wt_col=wt_col)
  or_est = get_theta0_est_or(data=data, predMu0=or0_info$predMu0, outcome_col=outcome_col, 
                             trt_col=trt_col, is_boot=is_boot, wt_col=wt_col)
  ipw_est = get_theta0_est_ipw(data=data, wts=ps0_wts, outcome_col=outcome_col, 
                               trt_col=trt_col, is_boot=is_boot, wt_col=wt_col)
  dr_est = get_theta0_est_dr(data=data, wts=ps0_wts, predMu0=or0_info$predMu0, 
                             outcome_col=outcome_col, trt_col=trt_col, is_boot=is_boot, 
                             wt_col=wt_col)
  Ymean = stats::weighted.mean(x=data[data[[trt_col]]==1, outcome_col], w=data[data[[trt_col]]==1, wt_col])
  ests = list(TWFE_ADJ=twfe_est, DIM=(Ymean-naive_est), 
              OR=(Ymean-or_est), IPW=(Ymean-ipw_est), DR=(Ymean-dr_est))
  
  return(list(ests=ests))
}

get_theta0_ests_w_mult = function(data, or_covs, ps_covs, is_boot=F, trt_col="A", 
                                  outcome_col="deltaY", outcome_col0="Y0", outcome_col1="Y1",
                                  wt_col="boot_wts", wt_norm=T, avg_ps_X=T, full_period_data=NULL) {
  ps0_info = get_ps0_info(data=data, ps_covs=ps_covs, trt_col=trt_col, is_boot=is_boot, 
                          wt_col=wt_col, avg_ps_X=avg_ps_X, full_period_data=full_period_data)
  or0_info = get_or0_info(data=data, or_covs=or_covs, outcome_col=outcome_col, 
                          trt_col=trt_col, is_boot=is_boot, wt_col=wt_col)
  ps0_wts = get_ps0_wts(data=data, ps=ps0_info$ps, trt_col=trt_col, normalize=wt_norm)
  
  twfe_est = get_theta0_est_twfe(data=data, or_covs=or_covs, outcome_col0=outcome_col0, 
                                 outcome_col1=outcome_col1, trt_col=trt_col, is_boot=is_boot, 
                                 wt_col=wt_col) 
  naive_est = get_theta0_est_naive(data=data, outcome_col=outcome_col, trt_col=trt_col, 
                                   is_boot=is_boot, wt_col=wt_col)
  or_est = get_theta0_est_or(data=data, predMu0=or0_info$predMu0, outcome_col=outcome_col, 
                             trt_col=trt_col, is_boot=is_boot, wt_col=wt_col)
  ipw_est = get_theta0_est_ipw(data=data, wts=ps0_wts, outcome_col=outcome_col, 
                               trt_col=trt_col, is_boot=is_boot, wt_col=wt_col)
  dr_est = get_theta0_est_dr(data=data, wts=ps0_wts, predMu0=or0_info$predMu0, 
                             outcome_col=outcome_col, trt_col=trt_col, is_boot=is_boot, 
                             wt_col=wt_col)
  Ymean = stats::weighted.mean(x=data[data[[trt_col]]==1, outcome_col], w=data[data[[trt_col]]==1, wt_col])
  Ymean0 = stats::weighted.mean(x=data[data[[trt_col]]==1, outcome_col0], w=data[data[[trt_col]]==1, wt_col])
  Ymean1 = stats::weighted.mean(x=data[data[[trt_col]]==1, outcome_col1], w=data[data[[trt_col]]==1, wt_col])
  ests_add = list(TWFE_ADJ=twfe_est, DIM=(Ymean-naive_est), 
                  OR=(Ymean-or_est), IPW=(Ymean-ipw_est), DR=(Ymean-dr_est))
  ests_mult = list(TWFE_ADJ=twfe_est, DIM=(Ymean1/(Ymean0 + naive_est)), 
                   OR=(Ymean1/(Ymean0 + or_est)), IPW=(Ymean1/(Ymean0 + ipw_est)), DR=(Ymean1/(Ymean0 + dr_est)))
  return(list(ests_add=ests_add, ests_mult=ests_mult))
}