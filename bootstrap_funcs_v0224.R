get_boot_wts = function(seed, data, block_col=NULL) {
  #-----------------------------------------------------------------------------
  set.seed(seed)
  if (is.null(block_col)) wts <- rexp(n=nrow(data), rate=1) #c(MCMCpack::rdirichlet(n=1, alpha=rep(1,nrow(data))))
  else {
    blocks = data.frame(block=unique(data[[block_col]]))
    blocks$wts = rexp(n=nrow(blocks), rate=1) #c(MCMCpack::rdirichlet(n=1, alpha=rep(1,nrow(blocks))))
    wts = c(sapply(data[[block_col]], function(x) blocks$wts[blocks$block == x]))
  }
  # normalize per treatment group (e.g., stratified bootstrap)
  wts[(data$A == 1) & (data$P == 1)] <- wts[(data$A == 1) & (data$P == 1)] * 
    sum((data$A == 1) & (data$P == 1)) / sum(wts[(data$A == 1) & (data$P == 1)])
  wts[(data$A == 0) & (data$P == 1)] <- wts[(data$A == 0) & (data$P == 1)] * 
    sum((data$A == 0) & (data$P == 1)) / sum(wts[(data$A == 0) & (data$P == 1)])
  wts[(data$hA ==  1) & (data$P == 0)] <- wts[(data$hA ==  1) & (data$P == 0)] * 
    sum((data$hA ==  1) & (data$P == 0)) / sum(wts[(data$hA ==  1) & (data$P == 0)])
  wts[(data$hA ==  0) & (data$P == 0)] <- wts[(data$hA ==  0) & (data$P == 0)] * 
    sum((data$hA ==  0) & (data$P == 0)) / sum(wts[(data$hA ==  0) & (data$P == 0)])
  #-----------------------------------------------------------------------------
  return(wts)
}

collect_boots = function(data, n_periods, nboots, or_covs, ps_covs, wt_norm0=T, block_col=NULL) {
  boot_res = NULL
  for (i in 1:nboots) {
    uniq_data = unique(data[,c("ID", "A", "hA", "P", "Zip")])
    uniq_data$boot_wts = get_boot_wts(seed=i, data=uniq_data, block_col=block_col)
    data = merge(data[,names(data) != "boot_wts"], uniq_data, by=c("ID", "A", "hA", "P", "Zip"))
    est_info = get_est_info(data=data, n_periods=n_periods, or_covs=or_covs, 
                            ps_covs=ps_covs, is_boot=T, wt_col="boot_wts", wt_norm0=wt_norm0, y1_col="Y1")
    est_dfs = list(theta_att=data.frame(est_info$theta_att), theta_atn=data.frame(est_info$theta_atn))
    if (is.null(boot_res)) boot_res = est_dfs
    else for (est in names(est_dfs)) boot_res[[est]] = rbind(boot_res[[est]], est_dfs[[est]])
  }
  agg_res = lapply(boot_res, function(x) data.frame(cbind(t(apply(x, 2, function(y) mean(y, na.rm=T))), 
                                                          t(apply(x, 2, function(y) sd(y, na.rm=T))),
                                                          t(apply(x, 2, function(y) quantile(y, probs=0.025, na.rm=T))),
                                                          t(apply(x, 2, function(y) quantile(y, probs=0.975, na.rm=T))))))
  for (est in names(boot_res)) names(agg_res[[est]]) = c(paste0(names(boot_res[[est]]), ".mean"), paste0(names(boot_res[[est]]), ".sd"), 
                                                         paste0(names(boot_res[[est]]), ".lower"), paste0(names(boot_res[[est]]), ".upper"))
  return(agg_res)
}

