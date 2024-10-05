library(dplyr)
library(stats)
library(MASS)
source("ctl_component_funcs_v0224.R")

flatten_list_res = function(ests_list) {
  ests_df = NULL
  for (mdl in names(ests_list)) {
    mdl_ests = data.frame(t(ests_list[[mdl]]))
    names(mdl_ests) = paste0(mdl, 1:ncol(mdl_ests))
    ests_df = if (is.null(ests_df)) mdl_ests else cbind(ests_df, mdl_ests)
  }
  return(ests_df)
}

get_est_info = function(data, n_periods, or_covs, ps_covs, is_boot=F, 
                        wt_col="boot_wts", wt_norm0=T, y1_col="Y1") {
  spec_res = list()
  for (m in 1:n_periods) {
    est_ctl_info_att = get_theta0_ests(data=data[(data$period == m) & (data$P == 1),], or_covs=or_covs, ps_covs=ps_covs,
                                   is_boot=is_boot, trt_col="A", wt_col=wt_col, wt_norm=wt_norm0,
                                   avg_ps_X=T, full_period_data=data[data$P==1,])$ests
    est_ctl_info_atn = get_theta0_ests(data=data[(data$period == m) & (data$P == 0),], or_covs=or_covs, ps_covs=ps_covs,
                                       is_boot=is_boot, trt_col="hA",  wt_col=wt_col, wt_norm=wt_norm0,
                                       avg_ps_X=T, full_period_data=data[data$P==0,])$ests
    est_info = list(theta_att=flatten_list_res(est_ctl_info_att), theta_atn=flatten_list_res(est_ctl_info_atn))
    spec_res = add_info(res=spec_res, info=est_info)
  }
  avg_est_info = list()
  for (est in names(spec_res)) {
    avg_df = colMeans(spec_res[[est]])
    avg_est_info[[est]] = list()
    for (model in names(est_info[[est]])) {
      avg_est_info[[est]][[model]] = avg_df[[model]]
    }
  }
  return(avg_est_info)
}

add_info = function(res, info) {
  dfs = list(theta_att=data.frame(info$theta_att), 
             theta_atn=data.frame(info$theta_atn))
  if (is.null(res)) res = dfs
  else for (lbl in names(dfs)) res[[lbl]] = rbind(res[[lbl]], dfs[[lbl]])
  return(res)
}

