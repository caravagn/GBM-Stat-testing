
# install.packages('VGAM', dependencies = T)
# library(VGAM)
# 
# install.packages('vcfR', dependencies = T)
# library(vcfR)
# 
# patient = 'gbm_Giulio/42mergedwes.vcf'
# ID = '42'
# 
# data = read.table(patient, sep = '\t', stringsAsFactors = FALSE, header = TRUE)
# head(data)
# 
# ## VCF codes
# WES_cols = paste(ID, c('M', 'S', 'T1', 'T2', 'T3', 'T4'), 'WES', sep = '')    ## Whole Exome
# TS1_cols = paste(ID, c('M', 'S'), '*', sep = '')                              ## Target Sequencing High Depth - Panel One
# TS2_cols = paste(ID, c('M', 'S', 'T1', 'T2', 'T3', 'T4', 'B'), '*', sep = '') ## Target Sequencing High Depth - Panel Two

## Example -- we fit data from the WES, and test on Deep sequencing
N_wes  = 300                          # Number of mutations detected in the WES sample
N_deep = 50                           # Number of mutations to test via Deep sequencing
wes_coverage_mean = 100               # Mean coverage
deep_coverage_mean = 3000             # Mean coverage
BetaBin = c(prob = 0.3, rho = 0.1)   # Beta-Binomial parameters for the true model
deep_purity = .7

repeat{
  # Coverage values for WES -- BetaBinomial samples around wes_coverage_mean
  wes_coverage_values = VGAM::rbetabinom(N_wes, 2 * wes_coverage_mean, prob = .5, rho = 0.1)
    
  # Mutant reads detected by WES -- Binomial samples with parameters BetaBin
  wes_rmutants_values = VGAM::rbetabinom(N_wes, wes_coverage_mean, prob = BetaBin['prob'], rho = BetaBin['rho'])

  if(all(wes_coverage_values - wes_rmutants_values > 0)) break
}

hist(wes_coverage_values)
hist(wes_rmutants_values)
hist(wes_rmutants_values/wes_coverage_values)


# ... just to test Tarone's Z statistics
# wes_rmutants_values = rbinom(N_wes, n, prob = .4)

# Deep sequencing -- Binomial samples around wes_coverage_mean, with slightly less variation than WES
deep_coverage_values = rbinom(N_deep, deep_coverage_mean, runif(1, min = 0.3, max = 0.9))


# 
# fitBinom = fitdist(data=wes_rmutants_values, dist="binom", fix.arg=list(size=100), start=list(prob=0.35))
# 
# require (MASS)
# require(ggpubr)
# qqp(wes_rmutants_values, "binom", size=100, fitBinom$estimate)
# 
# fit = VGAM::vglm(cbind(wes_rmutants_values, wes_coverage_values - wes_rmutants_values) ~ 1, betabinomial, trace = FALSE)
# plot(fit)
# 
# qqp(wes_rmutants_values, "binom", size=100, fitBinom$estimate)


# dev.new(heigth=22)
pplot_fit_power(s_n = wes_rmutants_values, 
                t_n = wes_coverage_values, 
                tests_t_n = deep_coverage_values,
                main = bquote('Simulated Beta-Binomial data with ' ~ N == .(N_wes) ~ ~ mu == .(BetaBin['prob']) ~ ',' ~ rho == .(BetaBin['rho'])),
                psign = 0.05,
                purity = .05,
                range_coverage = 10:4000
                )

