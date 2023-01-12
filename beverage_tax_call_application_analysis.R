source("beverage_tax_compile_functions.R")
library(rlist)
s_vars_att_or = c("WtPPUperOz", "PercentWhite", "Classification", "PreTaxTarget")
s_vars_att_ps = c("WtPPUperOz", "PercentWhite", "AverageHouseValue", "Classification")
s_vars_atn_or = c("WtPPUperOz", "IncomePerHousehold", "Classification", "PreTaxTarget")
s_vars_atn_ps = c("WtPPUperOz", "PercentBlack", "AverageHouseValue", "IncomePerHousehold", "Classification")

p_vars_att_or = c("WtPPUperOz", "PercentWhite", "AverageHouseValue", "IncomePerHousehold", "PreTaxTarget")
p_vars_att_ps = c("WtPPUperOz", "PercentWhite", "PercentBlack")
p_vars_atn_or = c("WtPPUperOz", "PercentBlack", "AverageHouseValue", "PreTaxTarget")
p_vars_atn_ps = c("WtPPUperOz", "PercentWhite", "AverageHouseValue")

#################### ATT #######################################################

# Parametric, Additive

## Supermarket

### Taxed

stmap_att = run_parametric_estimates(st, s_vars_att_ps, s_vars_att_or, get_season=T, manual=T, 
                                     ratio=F, get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")
stdap_att = run_parametric_estimates(st, s_vars_att_ps, s_vars_att_or, get_season=T, manual=F, ratio=F, 
                                     get_sens=F, trt_cols="Philadelphia", control_cols="Baltimore")

### Non-taxed

snmap_att = run_parametric_estimates(sn, s_vars_att_ps, s_vars_att_or, get_season=T, manual=T, ratio=F, 
                                     get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")
sndap_att = run_parametric_estimates(sn, s_vars_att_ps, s_vars_att_or, get_season=T, manual=F, 
                                     ratio=F, get_sens=F, trt_cols="Philadelphia", control_cols="Baltimore")

## Pharmacy

### Taxed

ptmap_att = run_parametric_estimates(pt, p_vars_att_ps, p_vars_att_or, get_season=T, manual=T, 
                                     ratio=F, get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")
ptdap_att = run_parametric_estimates(pt, p_vars_att_ps, p_vars_att_or, get_season=T, manual=F, 
                                     ratio=F, get_sens=F, trt_cols="Philadelphia", control_cols="Baltimore")

### Non-taxed

pnmap_att = run_parametric_estimates(pn, p_vars_att_ps, p_vars_att_or, get_season=T, manual=T, 
                                     ratio=F, get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")
pndap_att = run_parametric_estimates(pn, p_vars_att_ps, p_vars_att_or, get_season=T, manual=F, 
                                     ratio=F, get_sens=F, trt_cols="Philadelphia", control_cols="Baltimore")

# Bootstrap, Additive

## Supermarket

### Taxed

stmab_att = run_freq_bootstrap(st, s_vars_att_ps, s_vars_att_or, nboots=1000, get_season=T, manual=T, ratio=F, 
                               get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")

### Non-taxed

snmab_att = run_freq_bootstrap(sn, s_vars_att_ps, s_vars_att_or, nboots=1000, get_season=T, manual=T, ratio=F, 
                               get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")

# Bootstrap, Additive

## Pharmacy

### Taxed

ptmab_att = run_freq_bootstrap(pt, p_vars_att_ps, p_vars_att_or, nboots=1000, get_season=T, manual=T, ratio=F, 
                               get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")

### Non-taxed

pnmab_att = run_freq_bootstrap(pn, p_vars_att_ps, p_vars_att_or, nboots=1000, get_season=T, manual=T, ratio=F, 
                               get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")

# Bootstrap, Multiplicative

## Supermarket

### Taxed

stmmb_att = run_freq_bootstrap(st, s_vars_att_ps, s_vars_att_or, nboots=1000, get_season=T, manual=T, ratio=T, 
                               get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")

### Non-taxed

snmmb_att = run_freq_bootstrap(sn, s_vars_att_ps, s_vars_att_or, nboots=1000, get_season=T, 
                               manual=T, ratio=T, get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")

# Bootstrap, Multiplicative

## Pharmacy

### Taxed

ptmmb_att = run_freq_bootstrap(pt, p_vars_att_ps, p_vars_att_or, nboots=1000, get_season=T, manual=T, ratio=T, 
                               get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")

### Non-taxed

pnmmb_att = run_freq_bootstrap(pn, p_vars_att_ps, p_vars_att_or, nboots=1000, get_season=T, manual=T, ratio=T, 
                               get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")

#################### ATN #######################################################

# Parametric, Additive

## Supermarket

### Taxed

stmap_atn = run_parametric_estimates(st, s_vars_atn_ps, s_vars_atn_or, get_season=T, manual=T, 
                                     ratio=F, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")
stdap_atn = run_parametric_estimates(st, s_vars_atn_ps, s_vars_atn_or, get_season=T, manual=F, 
                                     ratio=F, get_sens=F, trt_cols="Border Counties", control_cols="Non-Border Counties")

### Non-taxed

snmap_atn = run_parametric_estimates(sn, s_vars_atn_ps, s_vars_atn_or, get_season=T, manual=T, 
                                     ratio=F, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")
sndap_atn = run_parametric_estimates(sn, s_vars_atn_ps, s_vars_atn_or, get_season=T, manual=F, 
                                     ratio=F, get_sens=F, trt_cols="Border Counties", control_cols="Non-Border Counties")

## Pharmacy

### Taxed

ptmap_atn = run_parametric_estimates(pt, p_vars_atn_ps, p_vars_atn_or, get_season=T, manual=T, 
                                     ratio=F, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")
ptdap_atn = run_parametric_estimates(pt, p_vars_atn_ps, p_vars_atn_or, get_season=T, manual=F, 
                                     ratio=F, get_sens=F, trt_cols="Border Counties", control_cols="Non-Border Counties")

### Non-taxed

pnmap_atn = run_parametric_estimates(pn, p_vars_atn_ps, p_vars_atn_or, get_season=T, manual=T, 
                                     ratio=F, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")
pndap_atn = run_parametric_estimates(pn, p_vars_atn_ps, p_vars_atn_or, get_season=T, manual=F, 
                                     ratio=F, get_sens=F, trt_cols="Border Counties", control_cols="Non-Border Counties")

# Bootstrap, Additive

## Supermarket

### Taxed

stmab_atn = run_freq_bootstrap(st, s_vars_atn_ps, s_vars_atn_or, nboots=1000, get_season=T, manual=T, 
                               ratio=F, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")

### Non-taxed

snmab_atn = run_freq_bootstrap(sn, s_vars_atn_ps, s_vars_atn_or, nboots=1000, get_season=T, manual=T, 
                               ratio=F, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")

# Bootstrap, Additive

## Pharmacy

### Taxed

ptmab_atn = run_freq_bootstrap(pt, p_vars_atn_ps, p_vars_atn_or, nboots=1000, get_season=T, manual=T, 
                               ratio=F, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")

### Non-taxed

pnmab_atn = run_freq_bootstrap(pn, p_vars_atn_ps, p_vars_atn_or, nboots=1000, get_season=T, manual=T,
                               ratio=F, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")

# Bootstrap, Multiplicative

## Supermarket

### Taxed

stmmb_atn = run_freq_bootstrap(st, s_vars_atn_ps, s_vars_atn_or, nboots=1000, get_season=T, manual=T, 
                               ratio=T, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")

### Non-taxed

snmmb_atn = run_freq_bootstrap(sn, s_vars_atn_ps, s_vars_atn_or, nboots=1000, get_season=T, manual=T, 
                               ratio=T, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")

# Bootstrap, Multiplicative

## Pharmacy

### Taxed

ptmmb_atn = run_freq_bootstrap(pt, p_vars_atn_ps, p_vars_atn_or, nboots=1000, get_season=T, manual=T, ratio=T, 
                               get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")

### Non-taxed

pnmmb_atn = run_freq_bootstrap(pn, p_vars_atn_ps, p_vars_atn_or, nboots=1000, get_season=T, manual=T, ratio=T, 
                               get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")


main_results = list(stmap_att=stmap_att, stdap_att=stdap_att, snmap_att=snmap_att, sndap_att=sndap_att,
                    ptmap_att=ptmap_att, ptdap_att=ptdap_att, pnmap_att=pnmap_att, pndap_att=pndap_att,
                    stmab_att=stmab_att, snmab_att=snmab_att, ptmab_att=ptmab_att, pnmab_att=pnmab_att,
                    stmmb_att=stmmb_att, snmmb_att=snmmb_att, ptmmb_att=ptmmb_att, pnmmb_att=pnmmb_att,
                    stmap_atn=stmap_atn, stdap_atn=stdap_atn, snmap_atn=snmap_atn, sndap_atn=sndap_atn,
                    ptmap_atn=ptmap_atn, ptdap_atn=ptdap_atn, pnmap_atn=pnmap_atn, pndap_atn=pndap_atn,
                    stmab_atn=stmab_atn, snmab_atn=snmab_atn, ptmab_atn=ptmab_atn, pnmab_atn=pnmab_atn,
                    stmmb_atn=stmmb_atn, snmmb_atn=snmmb_atn, ptmmb_atn=ptmmb_atn, pnmmb_atn=pnmmb_atn)
list.save(main_results, "beverage_tax_main_results_11_07_22.RData")

main_results$stdap_att$year$effects
main_results$stdap_att$year$variances
main_results$stmab_att$year$effects
main_results$stmab_att$year$lowers
main_results$stmab_att$year$uppers

main_results$ptdap_att$year$effects
main_results$ptdap_att$year$variances
main_results$ptmab_att$year$effects
main_results$ptmab_att$year$lowers
main_results$ptmab_att$year$uppers

main_results$stdap_atn$year$effects
main_results$stdap_atn$year$variances
main_results$stmab_atn$year$effects
main_results$stmab_atn$year$lowers
main_results$stmab_atn$year$uppers

main_results$ptdap_atn$year$effects
main_results$ptdap_atn$year$variances
main_results$ptmab_atn$year$effects
main_results$ptmab_atn$year$lowers
main_results$ptmab_atn$year$uppers

###### Pre-trends ##############################################################

# Parametric

stdap_att_pre = run_parametric_pretrends(st, s_vars_att_ps, s_vars_att_or, trt_cols="Philadelphia", control_cols="Baltimore")
stdap_atn_pre = run_parametric_pretrends(st, s_vars_atn_ps, s_vars_atn_or, trt_cols="Border Counties", control_cols="Non-Border Counties")

ptdap_att_pre = run_parametric_pretrends(pt, p_vars_att_ps, p_vars_att_or, trt_cols="Philadelphia", control_cols="Baltimore")
ptdap_atn_pre = run_parametric_pretrends(pt, p_vars_atn_ps, p_vars_atn_or, trt_cols="Border Counties", control_cols="Non-Border Counties")

# Bootstrap

stmab_att_pre = run_freq_bootstrap_pretrends(st, s_vars_att_ps, s_vars_att_or, nboots=1000, 
                                             trt_cols="Philadelphia", control_cols="Baltimore")
stmab_atn_pre = run_freq_bootstrap_pretrends(st, s_vars_atn_ps, s_vars_atn_or, nboots=1000, 
                                             trt_cols="Border Counties", control_cols="Non-Border Counties")

ptmab_att_pre = run_freq_bootstrap_pretrends(pt, p_vars_att_ps, p_vars_att_or, nboots=1000, 
                                             trt_cols="Philadelphia", control_cols="Baltimore")
ptmab_atn_pre = run_freq_bootstrap_pretrends(pt, p_vars_atn_ps, p_vars_atn_or, nboots=1000, 
                                             trt_cols="Border Counties", control_cols="Non-Border Counties")

pretrend_results = list(stdap_att_pre=stdap_att_pre, stdap_atn_pre=stdap_atn_pre, ptdap_att_pre=ptdap_att_pre, ptdap_atn_pre=ptdap_atn_pre,
                        stmab_att_pre=stmab_att_pre, stmab_atn_pre=stmab_atn_pre, ptmab_att_pre=ptmab_att_pre, ptmab_atn_pre=ptmab_atn_pre)
list.save(pretrend_results, "beverage_tax_pretrend_results_11_07_22.RData")


###### Bayes Boot ##############################################################

# Additive

stmab_att_bayes = run_bayes_bootstrap_wt(st, s_vars_att_ps, s_vars_att_or, nboots=10, get_season=T, manual=T, 
                                         ratio=F, get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")
stmab_atn_bayes = run_bayes_bootstrap_wt(st, s_vars_atn_ps, s_vars_atn_or, nboots=1000, get_season=T, manual=T, 
                                         ratio=F, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")

ptmab_att_bayes = run_bayes_bootstrap_wt(pt, p_vars_att_ps, p_vars_att_or, nboots=1000, get_season=T, manual=T, 
                                         ratio=F, get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")
ptmab_att_bayes_sample = run_bayes_bootstrap(pt, p_vars_att_ps, p_vars_att_or, nboots=500, get_season=T, manual=T, 
                                         ratio=F, get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")
ptmab_atn_bayes = run_bayes_bootstrap_wt(pt, p_vars_atn_ps, p_vars_atn_or, nboots=1000, get_season=T, manual=T, 
                                         ratio=F, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")

# Multiplicative

stmmb_att_bayes = run_bayes_bootstrap_wt(st, s_vars_att_ps, s_vars_att_or, nboots=1000, get_season=T, manual=T, 
                                         ratio=T, get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")
stmmb_atn_bayes = run_bayes_bootstrap_wt(st, s_vars_atn_ps, s_vars_atn_or, nboots=1000, get_season=T, manual=T, 
                                         ratio=T, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")

ptmmb_att_bayes = run_bayes_bootstrap_wt(pt, p_vars_att_ps, p_vars_att_or, nboots=1000, get_season=T, manual=T, 
                                         ratio=T, get_sens=T, trt_cols="Philadelphia", control_cols="Baltimore")
ptmmb_atn_bayes = run_bayes_bootstrap_wt(pt, p_vars_atn_ps, p_vars_atn_or, nboots=1000, get_season=T, manual=T, 
                                         ratio=T, get_sens=T, trt_cols="Border Counties", control_cols="Non-Border Counties")

bayes_results = list(stmab_att_bayes=stmab_att_bayes, stmab_atn_bayes=stmab_atn_bayes, 
                     ptmab_att_bayes=ptmab_att_bayes, ptmab_atn_bayes=ptmab_atn_bayes,
                     stmmb_att_bayes=stmmb_att_bayes, stmmb_atn_bayes=stmmb_atn_bayes, 
                     ptmmb_att_bayes=ptmmb_att_bayes, ptmmb_atn_bayes=ptmmb_atn_bayes,
                     ptmab_att_bayes_sample=ptmab_att_bayes_sample)
list.save(bayes_results, "beverage_tax_bayes_main_results_11_07_22.RData")


# Pre trends
stmab_att_pre_bayes = run_bayes_bootstrap_wt_pretrends(st, s_vars_att_ps, s_vars_att_or, nboots=1000, 
                                                       trt_cols="Philadelphia", control_cols="Baltimore")
stmab_atn_pre_bayes = run_bayes_bootstrap_wt_pretrends(st, s_vars_atn_ps, s_vars_atn_or, nboots=1000, 
                                                       trt_cols="Border Counties", control_cols="Non-Border Counties")

ptmab_att_pre_bayes = run_bayes_bootstrap_wt_pretrends(pt, p_vars_att_ps, p_vars_att_or, nboots=1000, 
                                                       trt_cols="Philadelphia", control_cols="Baltimore")
ptmab_atn_pre_bayes = run_bayes_bootstrap_wt_pretrends(pt, p_vars_atn_ps, p_vars_atn_or, nboots=1000, 
                                                       trt_cols="Border Counties", control_cols="Non-Border Counties")

bayes_results_pre = list(stmab_att_pre_bayes=stmab_att_pre_bayes, stmab_atn_pre_bayes=stmab_atn_pre_bayes, 
                         ptmab_att_pre_bayes=ptmab_att_pre_bayes, ptmab_atn_pre_bayes=ptmab_atn_pre_bayes)#,
                         # ptmab_att_bayes_sample=ptmab_att_bayes_sample)
list.save(bayes_results_pre, "beverage_tax_bayes_pretrend_results_11_07_22.RData")



