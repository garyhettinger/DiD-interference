library(MASS)
library(dplyr)
library(FNN)
library(ggplot2)

expit = function(x) {return(exp(x)/(1+exp(x)))}

# enforces heterogeneity of confounder X4 across m-time
period_hetero_X = function(M) {
  if (M %in% c(2, 4)) return(seq(1, 1.5, length.out=M))
  else if (M == 12) return(rep(seq(1, 1.5, length.out=4), 3)) 
  warning("Bad n_periods")
  return(NULL)
}

# eta_m^(1) (enforces heterogeneity of treatment effects across m-time)
period_hetero_eff = function(M) {
  if (M %in% c(2, 4)) return((sin(pi*(1:M)/4)+2)/mean(sin(pi*(1:M)/4)+2))
  else if (M == 12) return(rep((sin(pi*(1:M)/4)+2)/mean(sin(pi*(1:M)/4)+2),3))
  warning("Bad n_periods")
  return(NULL)
}

# eta_m^(0) (enforces heterogeneity of confounder effects across m-time)
period_hetero_conf = function(M) {
  if (M %in% c(2, 4)) return((cos(pi*(1:M)/4)+2)/mean(cos(pi*(1:M)/4)+2))
  else if (M == 12) return(rep((cos(pi*(1:M)/4)+2)/mean(cos(pi*(1:M)/4)+2),3))
  warning("Bad n_periods")
  return(NULL)
}

# mu_{m} among treated (Y1 outcome functions among C=1)
sim_trt_func = function(X, D, eff_scale=1, conf_scale=1) {
  mu00.true = conf_scale*(-0.43 + (X %*% c(-1.07, -0.95, 0.51, 0.70)))
  mu10.true = mu00.true + eff_scale*(-3.25*(1-D))
  return(list(mu10=mu10.true, mu00=mu00.true))
}

# mu_{m} among neighboring control (Y1 outcome functions among C=0)
sim_ngh_func = function(X, D, eff_scale=1, conf_scale=1) {
  mu00.true = conf_scale*(-1.21 + (X %*% c(-0.33, 0.48, -0.36, -0.33)))
  mu01.true = mu00.true + eff_scale*(0.3 + 0.75*as.integer(D < 0.3) + 0.5*as.integer(D < 0.65))
  return(list(mu01=mu01.true, mu00=mu00.true))
}


gen_sim = function(seed, n_trt, n_ngh, n_periods, for_truth=F) {
  set.seed(seed)
  n = n_trt + n_ngh
  
  # Separated X so can fix size of comparisons
  Sigma.trt = rbind(c(1, -0.35, -0.21, 0),
                c(-0.35, 1, 0.77, 0.50),
                c(-0.21, 0.77, 1, 0.45),
                c(0, 0.50, 0.45, 1))
  X.trt = mvrnorm(n=n_trt, mu=rep(0,4), Sigma=diag(4))#Sigma.trt)
  Sigma.ngh = rbind(c(1, -0.38, -0.38, -0.20),
                    c(-0.38, 1, 0.77, 0.42),
                    c(-0.38, 0.77, 1, 0.47),
                    c(-0.20, 0.42, 0.47, 1))
  X.ngh = mvrnorm(n=n_ngh, mu=rep(0,4), Sigma=diag(4))#Sigma.ngh)
  
  # Add time-varying covariate
  X.hetero = period_hetero_X(M=n_periods)
  for (j in 2:n_periods) {
    X.trt = rbind(X.trt, cbind(tail(X.trt[,1:3], n_trt), rnorm(n=n_trt, mean=tail(X.hetero[j]*X.trt[,4], n_trt), sd=0.1)))
    X.ngh = rbind(X.ngh, cbind(tail(X.ngh[,1:3], n_ngh), rnorm(n=n_ngh, mean=tail(X.hetero[j]*X.ngh[,4], n_ngh), sd=0.1)))
  }
  X = rbind(X.trt, X.ngh)
  
  # A and D based on averages per unit since not time-varying
  trt_ids = rep(1:n_trt, n_periods)
  X.trt.mean = as.matrix(aggregate(X.trt ~ trt_ids, FUN=mean)[,-1]) # check if this keeps order
  ngh_ids = rep(1:n_ngh, n_periods)
  X.ngh.mean = as.matrix(aggregate(X.ngh ~ ngh_ids, FUN=mean)[,-1]) # check if this keeps order
  
  # Simulate intervention group indicator within groups
  A.mean.trt = expit(cbind(1,X.trt.mean) %*% c(-0.11, 0.40, -0.66, 0.39, -0.28))
  A.trt.head = rbinom(n=n_trt, size=1, prob=A.mean.trt)
  A.trt = rep(A.trt.head, n_periods)
  A.mean.ngh = expit(cbind(1,X.ngh.mean) %*% c(-0.27, 0.33, 0.43, -0.26, -0.31))
  A.ngh.head = rbinom(n=n_ngh, size=1, prob=A.mean.ngh)
  A.ngh = rep(A.ngh.head, n_periods)
  
  # Simulate distances within groups
  D.mean.trt = cbind(1,X.trt.mean) %*% c(0.2, 0.3, 0, -0.4, 0.3) #c(0, 0.3, 0.35, -0.4, 0.3)
  D.trt.head <- expit(rnorm(n=n_trt, mean=D.mean.trt, sd=0.5))
  D.trt = rep(D.trt.head, n_periods)
  D.mean.ngh = cbind(1,X.ngh.mean) %*% c(0.5, 0, 0.45, 0, -0.25) #c(0.5, 0.4, 0.45, -0.3, -0.25)
  D.ngh.head <- expit(rnorm(n=n_ngh, mean=D.mean.ngh, sd=0.5))
  D.ngh = rep(D.ngh.head, n_periods)
  
  if (for_truth) {
    df = data.frame(ID = c(rep(1:n_trt, n_periods), rep((n_trt+1):n, n_periods)), 
                    period=c(rep(1:n_periods, each=n_trt), rep(1:n_periods, each=n_ngh)), 
                    P=c(rep(1,n_trt*n_periods), rep(0, n_ngh*n_periods)),
                    D=c(D.trt, D.ngh), A=c(A.trt, rep(0, n_ngh*n_periods)), hA=c(rep(0, n_trt*n_periods), A.ngh),
                    X1=X[,1], X2=X[,2], X3=X[,3], X4=X[,4])
    return(df[order(df$period, df$ID),])
  }
  
  # Set zip code clusters
  D.trt.norm.jitt = rnorm(n=n_trt, mean=(D.trt.head - mean(D.trt.head))/
    sd(D.trt.head), sd=0.2)
  X.trt.jitt = apply(X.trt, 2, function(x) rnorm(n=length(x), mean=x, sd=0.2))
  clust.trt.mat = cbind(4*D.trt.norm.jitt[A.trt.head == 1], X.trt.jitt[1:n_trt,][A.trt.head == 1,])
  D.ngh.norm.jitt = rnorm(n=n_ngh, mean=(D.ngh.head - mean(D.ngh.head))/
    sd(D.ngh.head), sd=0.2)
  X.ngh.jitt = apply(X.ngh, 2, function(x) rnorm(n=length(x), mean=x, sd=0.2))
  clust.ngh.mat = cbind(4*D.ngh.norm.jitt[A.ngh.head == 1], X.ngh.jitt[1:n_ngh,][A.ngh.head == 1,])
  n.zip.trt.1 = n_trt / 5 
  n.zip.trt.0 = n_trt / 5 
  n.zip.ngh.1 = n_ngh / 5 
  n.zip.ngh.0 = n_ngh / 5 
  zip.trt.1.map = split(D.trt.norm.jitt[A.trt.head == 1], ggplot2::cut_number(x=D.trt.norm.jitt[A.trt.head == 1], n=n.zip.trt.1, labels=FALSE))
  zip.trt.0.map = split(D.trt.norm.jitt[A.trt.head == 0], ggplot2::cut_number(x=D.trt.norm.jitt[A.trt.head == 0], n=n.zip.trt.0, labels=FALSE))
  zip.ngh.1.map = split(D.ngh.norm.jitt[A.ngh.head == 1], ggplot2::cut_number(x=D.ngh.norm.jitt[A.ngh.head == 1], n=n.zip.ngh.1, labels=FALSE))
  zip.ngh.0.map = split(D.ngh.norm.jitt[A.ngh.head == 0], ggplot2::cut_number(x=D.ngh.norm.jitt[A.ngh.head == 0], n=n.zip.ngh.0, labels=FALSE))
  Zip.trt.head = rep(NA, n_trt)
  for (idx in names(zip.trt.1.map)) Zip.trt.head = ifelse((A.trt.head == 1) & (D.trt.norm.jitt %in% zip.trt.1.map[[idx]]), 
                                                          as.integer(idx), Zip.trt.head)
  for (idx in names(zip.trt.0.map)) Zip.trt.head = ifelse((A.trt.head == 0) & (D.trt.norm.jitt %in% zip.trt.0.map[[idx]]), 
                                                          as.integer(idx)+n.zip.trt.1, Zip.trt.head)
  Zip.ngh.head = rep(NA, n_ngh)
  for (idx in names(zip.ngh.1.map)) Zip.ngh.head = ifelse((A.ngh.head == 1) & (D.ngh.norm.jitt %in% zip.ngh.1.map[[idx]]), 
                                                          as.integer(idx)+n.zip.trt.1+n.zip.trt.0, Zip.ngh.head)
  for (idx in names(zip.ngh.0.map)) Zip.ngh.head = ifelse((A.ngh.head == 0) & (D.ngh.norm.jitt %in% zip.ngh.0.map[[idx]]), 
                                                          as.integer(idx)+n.zip.trt.1+n.zip.trt.0+n.zip.ngh.1, Zip.ngh.head)
  
  Zip.trt = rep(Zip.trt.head, n_periods)
  Zip.ngh = rep(Zip.ngh.head, n_periods)
  
  # Set spillover mapping
  D.kmeans.trt.head = rep(NA, n_trt)
  D.kmeans.trt.head[A.trt.head == 1] = Zip.trt.head[A.trt.head==1]
  D.kmeans.trt = rep(D.kmeans.trt.head, n_periods)
  D.kmeans.ngh.head = rep(NA, n_ngh)
  D.kmeans.ngh.head[A.ngh.head == 1] = Zip.ngh.head[A.ngh.head==1] - (n.zip.trt.1 + n.zip.trt.0)
  D.kmeans.ngh = rep(D.kmeans.ngh.head, n_periods)
  
  # Set full data.frame now
  full_df = data.frame(ID = c(rep(1:n_trt, n_periods), rep((n_trt+1):n, n_periods)), 
                  period=c(rep(1:n_periods, each=n_trt), rep(1:n_periods, each=n_ngh)), 
                  Zip=c(Zip.trt, Zip.ngh), D=c(D.trt, D.ngh), 
                  P=c(rep(1, n_trt*n_periods), rep(0, n_ngh*n_periods)), 
                  A=c(A.trt, rep(0, n_ngh*n_periods)), 
                  hA=c(rep(0, n_trt*n_periods), A.ngh),
                  fD=c(D.kmeans.trt, D.kmeans.ngh), X1=X[,1], X2=X[,2], X3=X[,3], X4=X[,4])
  full_df = full_df[order(full_df$period, full_df$ID),]
  
  # Set within-group outcome correlation random intercepts
  Y0.re = rnorm(n)
  delt.re = rnorm(n)
  conf_hetero = period_hetero_conf(M=n_periods)
  eff_hetero = period_hetero_eff(M=n_periods)
  
  # Simulate each period
  for (m in 1:n_periods) {
    df = full_df[full_df$period == m,]
    X = as.matrix(df[,paste0("X", 1:4)])
    
    # Set trt and ngh covariance: correlated within zip, even after adjusting for covariates
    Zip.mat = matrix(as.integer(matrix(rep(df$Zip, n), nrow=n, byrow=F) == matrix(rep(df$Zip, n), nrow=n, byrow=T)), nrow=n) * (1-diag(n))
    fD.vec = as.integer(matrix(rep(df$fD, n), nrow=n, byrow=F) == matrix(rep(df$fD, n), nrow=n, byrow=T))
    fD.mat = matrix(ifelse(is.na(fD.vec), 0, fD.vec), nrow=n) * (1-diag(n))
    trt.mat = matrix(as.integer(matrix(rep(df$P == 1, n), nrow=n, byrow=F) | matrix(rep(df$P == 1, n), nrow=n, byrow=T)), nrow=n)
    ngh.mat = matrix(as.integer(matrix(rep(df$P == 0, n), nrow=n, byrow=F) | matrix(rep(df$P == 0, n), nrow=n, byrow=T)), nrow=n)
    
    sigsq.y0.trt = 4
    rho.y0.trt = 0.2
    sigsq.y0.ngh = 2
    rho.y0.ngh = 0.2
    Sigma.Y0 = sigsq.y0.trt*(diag(n)*trt.mat + rho.y0.trt*trt.mat*Zip.mat) + 
      sigsq.y0.ngh*(diag(n)*ngh.mat + rho.y0.ngh*ngh.mat*Zip.mat)
    
    sigsq.delt.trt = 1
    rho.delt.trt = 0.2
    sigsq.delt.ngh = 1
    rho.delt.ngh = 0.2
    rho.delt.fD = -0.35
    Sigma.delt = sigsq.delt.trt*(diag(n)*trt.mat + rho.delt.trt*trt.mat*Zip.mat) + 
      sigsq.delt.ngh*(diag(n)*ngh.mat + rho.delt.ngh*ngh.mat*Zip.mat) +
      rho.delt.trt*sigsq.delt.ngh*fD.mat
    
    # Set Y0
    Y0.mean = ifelse(df$P == 1, ifelse(df$A == 1, conf_hetero[m]*4.63, 0) + 
                       cbind(1,X) %*% (conf_hetero[m]*c(10.48, 1, -3.52, 0.61, 1.28)),
                     ifelse(df$hA == 1, conf_hetero[m]*1.75, 0) + 
                       cbind(1,X) %*% (conf_hetero[m]*c(11.08, 1, -1.57, -1.28, -0.92)))
    Y0 = mvrnorm(n=1, mu=Y0.mean, Sigma=Sigma.Y0) + Y0.re
    full_df[full_df$period==m,"Y0"] = Y0
    
    # Set Y1 with heterogeneous by distance for deltaY for trt
    sim_trts = sim_trt_func(X=X, D=df$D, 
                            eff_scale=eff_hetero[m], conf_scale=conf_hetero[m])
    sim_nghs = sim_ngh_func(X=X, D=df$D, 
                            eff_scale=eff_hetero[m], conf_scale=conf_hetero[m])
    deltaY.mean = ifelse(df$P == 1, ifelse(df$A == 1, sim_trts$mu10, sim_trts$mu00),
                         ifelse(df$hA == 1, sim_nghs$mu01, sim_nghs$mu00))
    deltaY = mvrnorm(n=1, mu=deltaY.mean, Sigma=Sigma.delt) + delt.re
    full_df[full_df$period==m,"deltaY"] = deltaY
    full_df[full_df$period==m,"Y1"] = Y0 + deltaY
  }
  W = cbind( exp(full_df$X1/2) , 10+full_df$X2/(1+exp(full_df$X1)) , 
             (full_df$X1*full_df$X3/25 + 0.6)^3 , (full_df$X2+full_df$X4+20)^2)
  W = apply(W, 2, function(x) (x-mean(x))/sd(x))
  full_df$W1 = W[,1]
  full_df$W2 = W[,2]
  full_df$W3 = W[,3]
  full_df$W4 = W[,4]
  full_df$boot_wts = 1
  return(full_df)
}

get_true_params = function(n_periods, n=1000000) {
  sim_data = gen_sim(seed=1, n_trt=n, n_ngh=n, n_periods=n_periods, for_truth=T)
  
  eff_scale = period_hetero_eff(M=n_periods)
  conf_scale = period_hetero_conf(M=n_periods)
  eff_scales = eff_scale[sim_data$period]
  conf_scales = conf_scale[sim_data$period]
  
  trt_data = sim_data[sim_data$A == 1,]
  trt_X = as.matrix(trt_data[,paste0("X", 1:4)])
  sim_trts = sim_trt_func(X=trt_X, D=trt_data$D, eff_scale=eff_scales[sim_data$A==1], 
                          conf_scale=conf_scales[sim_data$A==1])
  att = mean(sim_trts$mu10 - sim_trts$mu00)
  
  ngh_data = sim_data[sim_data$hA == 1,]
  ngh_X = as.matrix(ngh_data[,paste0("X", 1:4)])
  sim_nghs = sim_ngh_func(X=ngh_X, D=ngh_data$D, eff_scale=eff_scales[sim_data$hA==1], 
                          conf_scale=conf_scales[sim_data$hA==1])
  atn = mean(sim_nghs$mu01 - sim_nghs$mu00)
  return(list(att=att, atn=atn, mu10_att=mean(sim_trts$mu10), mu00_att=mean(sim_trts$mu00), 
              mu01_atn=mean(sim_nghs$mu01), mu00_atn=mean(sim_nghs$mu00)))
}

