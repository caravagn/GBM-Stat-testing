pplot_fit_power = function(s_n, t_n, tests_t_n,
                      main, 
                      psign = 0.05,
                      purity = 1,
                      range_coverage = 10:5000)
{
  # Multi-layout
  layout(matrix(c(1, 6, 1, 6, 2, 3, 4, 5, 4, 5), 5, 2, byrow = TRUE))
  
  require(VGAM)
  
  # Fitting a Beta-Binomial to data
  cat('*** Beta-Binomial model fit\n')  
  fit = Coef(VGAM::vglm(cbind(s_n, t_n - s_n) ~ 1, betabinomial, trace = FALSE))
  fit_prob = round(fit[1], 3)
  fit_disp = round(fit[2], 3)
  print(fit)  

  # empirical density  
  N_emp = round(mean(t_n))
  x = round(seq(0, 1, 0.01) * N_emp)
  bbin_den = VGAM::dbetabinom(x, N_emp, prob = fit_prob, rho = fit_disp, log = FALSE)
  cat("*** Empirical coverage WES", N_emp, '\n')
  
  # Fitting a Binomial to data
  cat('*** Binomial model fit with Empirical coverage WES \n')  
  
  require(fitdistrplus)
  fitBinom = fitdist(data = s_n, dist="binom", fix.arg=list(size=N_emp), start=list(prob=.3))
  print(fitBinom$estimate)  
  bin_den = dbinom(x, N_emp, fitBinom$estimate, log = FALSE)
  
  # density and Histogram of success probability
  h = hist(s_n / t_n, breaks = seq(0, 1.1, 0.01), plot = FALSE)
  h$counts=h$counts/sum(h$counts)
  
  plot(h,
    col = 'lightblue',
    main = main,
    xlab = 'VAF',
    border = NA
  )
  
  lines(x/N_emp, bbin_den, lwd = 2, col = 'darkred')
  lines(x/N_emp, bin_den, lwd = 2, lty = 2, col = 'darkblue')
  
  legend(
    "topright",
    c('data', "Beta-Binomial (MLE fit)", "Binomial (MLE fit)"),
    col = c('lightblue', 'darkred', "darkblue"),
    lwd = 2,
    bty = 'n'
  )
  
  legend(
    "right",
    title = 'MLE fit',
    legend = 
      sapply(c(
      bquote(mu == .(fit_prob) ~ ',' ~ rho == .(fit_disp) ~ C[symbol("\052")] == .(N_emp)),
      bquote(p == .(fitBinom$estimate) ~ C[symbol("\052")] == .(N_emp))
      ), as.expression),
    
    col = c('darkred', 'darkblue'),
    bty = 'n',
    pch = 19,
    cex = 1
  )
  
  # Histogram of trials
  hist(t_n,
     col = 'lightgray',
     freq = F,
     breaks = 22,
     border = NA,
     xlab = '',
     main = 'Coverage distribution')
  
  abline(v = mean(t_n), col = "royalblue", lty = 2)
  abline(v = median(t_n), col = "red", lty = 2)
  legend(x = "topright", c("Mean", "Median"), col = c("royalblue", "red"), lwd = c(2, 2), cex = .6, bty = 'n')
  
  
  # Histogram of successes
  hist(s_n,
     col = 'orange',
     freq = F,
     border = NA,
     breaks = 22,
     xlab = '',
     main = 'Mutants distribution')
  
  abline(v = mean(s_n), col = "royalblue", lty = 2)
  abline(v = median(s_n), col = "red", lty = 2)
  legend(x = "topright", c("Mean", "Median"), col = c("royalblue", "red"), lwd = c(2, 2), cex = .6, bty = 'n')
 
  # P-values for the null model  
  pvalues = VGAM::dbetabinom(0, range_coverage, prob = fit_prob, rho = fit_disp, log = FALSE)
  names(range_coverage) = range_coverage
  min_coverage_sign = range_coverage[min(which(pvalues < psign))]
  
  # Bonferroni correction
  Num_tests = length(tests_t_n)
  min_coverage_sign_FWER = range_coverage[min(which(pvalues < psign/Num_tests))]
  
  plot(pvalues,
       xlab = 'log(Coverage c)',
       ylab = bquote(log( H[0] ~ ':' ~ 'BetaBin[ 0 | c,' ~ mu ~','~ rho ~']' )),
       main = bquote(bold('P-values with FWER:') ~ alpha ~ '=' ~  .(psign) ~ ','~ .(Num_tests) ~ 'tests,'~ pi == .(purity)),
       type = 'l',
       log = 'xy')
  
  abline(h = psign, col = "red", lty = 2)
  abline(h = psign/Num_tests, col = "orange", lty = 2)
  
  abline(v = min_coverage_sign, col = "blue", lty = 2)
  abline(v = min_coverage_sign_FWER, col = "black", lty = 2) 
  
  legend('topright', 
         legend = 
           sapply(c(
             bquote(p[symbol("\052")] ~ '<' ~  .(psign)), 
             bquote(c[symbol("\052")] ~ '=' ~  .(min_coverage_sign)),
             bquote(p[symbol("\052")]^{FWER} ~ '<' ~  .(psign/Num_tests)), 
             bquote(c[symbol("\052")]^{FWER} ~ '=' ~  .(min_coverage_sign_FWER))
             ), as.expression),
         col = c('red', 'blue', 'orange', 'black'), 
         lty = 2, 
         box.lwd = 0,
         box.col = "white",
         bg = "white"
         ) 
  
    cat('*** Mulitple Hypothesis Testing\n\tAlpha:', psign, 'Tests:', Num_tests, 'is p =', psign/Num_tests, '\n',
        '\tMinimum coverage for H_0 < p:', min_coverage_sign_FWER, '\n',
        '\tPurity correction: ', purity, '\n'
        )  

    mean_precorrections = mean(tests_t_n)
    tests_t_n = tests_t_n - (1-purity) * tests_t_n
    mean_postcorrections =  mean(tests_t_n)
    
    cat('*** Purity correction (mean testing coverage) :',
        mean_postcorrections, '-- it was', mean_precorrections, '\n')  
    
    
    rejected = tests_t_n[tests_t_n <= min_coverage_sign_FWER]
    points(rejected, pvalues[rejected], col = 'darkred', pch = 19)
  
    # cat('*** Rejected FWER :', rejected, '\n') 
    cat('*** Rejected FWER : ', length(rejected), ' mutations\n') 
    
    accepted = tests_t_n[tests_t_n > min_coverage_sign_FWER]
    points(accepted, pvalues[accepted], col = 'darkgreen', pch = 19)
    
    # cat('*** Accepted FWER :', accepted, '\n')
    cat('*** Accepted FWER :', length(accepted), 'mutations\n')
    
    # Histogram of trials for test
    h = hist(tests_t_n, breaks = 22, plot = FALSE)
    cuts = cut(h$breaks, c(-Inf, min_coverage_sign_FWER, Inf))
    
    # print(h$breaks)
    # print(c('darkred', 'darkgreen')[cuts])
    
    plot(h,
         col = c('darkred', 'darkgreen')[cuts], #c('darkred', 'darkgreen'),
         # freq = F,
         # breaks = 22,
         border = NA,
         xlab = '',
         main = 
           bquote(bold('Coverage distribution corrected for') ~ pi),
         sub = paste('Mean: pre-correction', round(mean_precorrections), ' -- post-correction', round(mean_postcorrections))
    )
   
     legend('topright', 
           legend = c(length(rejected), length(accepted)),
           col = c('darkred', 'darkgreen'),
           bty = 'n',
           pch = 19
    )
     
    
    ret = list(accepted = pvalues[accepted], rejected = pvalues[rejected])
    
    ## Tarone's Z statistics for overdispersion Beta Binomial
    cat('*** Bootstrapping Tarone\'s Z statistics for overdispersion.')
    M <- 1000
    alt_hyp = null_hyp <- vector("numeric", length = M)
    for (i in 1:M){
      boostrap_samples = sample(length(t_n), replace = TRUE)
      
      t = t_n[boostrap_samples]
      s = s_n[boostrap_samples]
      
      p_hat = sum(s) / sum(t)
      
      S = sum( (s - t * p_hat)^2 / (p_hat * (1 - p_hat)) )
      Z_score = (S - sum(t)) / sqrt(2 * sum(t * (t - 1)))

      null_hyp[i] <- Z_score
    }
    
    H0 = dnorm(seq(-5,5,0.01), 0, 1)
    
    hist(null_hyp, 
         breaks = 100, 
         col = 'khaki',
         border = NA,
         freq = F,
         xlim = c(-5, max(null_hyp)),
         ylim = c(0, max(H0)),
         xlab = 'Z score',
         main = paste('Tarone\'s Z statistics\n for overdispersion (M=', M, ')', sep ='') 
         )
    
   lines(seq(-5,5,0.01), H0, lty = 2, col = 'red')
   legend('topright', 
          sapply(c(
            bquote(H[0] ~ ': N(0,1)'), 
            'Beta-Binomial'
          ), as.expression),
          col = c('red', 'khaki'),
          bty = 'n',
          pch = 19
          )
    
    
    
    return(ret)
}



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
BetaBin = c(prob = 0.35, rho = 0.15)   # Beta-Binomial parameters for the true model
deep_purity = .7

# Coverage values for WES -- Binomial samples around wes_coverage_mean
wes_coverage_values = rbinom(N_wes, wes_coverage_mean, runif(1, min = 0.3, max = 0.8))
# Mutant reads detected by WES -- Binomial samples with parameters BetaBin
n = 1
wes_rmutants_values = VGAM::rbetabinom(N_wes, wes_coverage_mean, prob = BetaBin['prob'], rho = BetaBin['rho'])
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

