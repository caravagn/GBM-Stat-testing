jamPDF = function(in.files, out.file = 'jamPDF.pdf', layout = '3x3', delete.original = TRUE, hide.output = TRUE)
{
  in.files = in.files[sapply(in.files, file.exists)]
  if(length(in.files) == 0) return()
  
  cmd = paste('pdfjam ', 
              paste(in.files, collapse = ' '), 
              ' --nup ', layout, ' --landscape --outfile ', out.file, sep = ' ')
  
  aa = system(cmd, intern = hide.output, ignore.stderr = TRUE)
  
  # print(cmd)
  
  if(delete.original) file.remove(setdiff(in.files, out.file))
}
correctReadCounts = function(data, column, purity)
{
  tumourSpecificReads = function(coverage, Ct = 2, purity = 1) {
    
    #Ct = tumour copy number
    #purity = tumour purity
    #coverage = total coverage at the locus
    
    #How reads are from the tumour?
    tumour.reads = (Ct*purity) / ((Ct*purity) + (2*(1 - purity))) * coverage
    
    #Return
    return(tumour.reads)
    
  }
  
  for(i in 1:nrow(data)){
    if(!is.na(data[i, column]))
      data[i, column] =  ceiling(tumourSpecificReads(data[i, column], data$CN[i], purity))
  }
  
  return(data)
}

BBMLE = function(data, column, main){
  require(VGAM)
  require(fitdistrplus)
  require(crayon)
  
  cat(bgRed(column), '\n')
  col.R = paste('TRNR', column, sep = '.') 
  col.V = paste('TRNV', column, sep = '.') 
  
  
  df = lapply(1:length(data), function(w, m) {
    cat(blue('CN = ', names(data)[w], '\n'))
    m = paste(m, '- CN =', names(data)[w], ' - npoints =', nrow( data[[w]]), 'from sample', column)
    
    w = data[[w]]
    s_n = w[, col.V]
    t_n = w[, col.R]
    xx  = cbind(s_n, t_n)
    
    # cat(col.V, '\n')
    # print(xx)
    # 
    if(any(t_n < s_n)){
      cat(red("\n----- Correction for entries\n"))
      print(xx[t_n < s_n])
      cat(red("\n-----\n"))
      t_n[t_n < s_n] = s_n[t_n < s_n]
    }
    
    # MLE BetaBin    
    cat(bgGreen('*** BetaBinomial MLE fit:') )
    fit = Coef(VGAM::vglm(cbind(s_n, t_n - s_n) ~ 1, betabinomial, trace = FALSE))
    fit_prob = round(fit[1], 3)
    fit_disp = round(fit[2], 3)
    cat(cyan('  mu ='), fit['mu'], cyan('rho ='), fit['rho'], '\n')  
    
    # empirical density  
    N_emp = round(mean(t_n))
    x = round(seq(0, 1, 0.01) * N_emp)
    bbin_den = VGAM::dbetabinom(x, N_emp, prob = fit_prob, rho = fit_disp, log = FALSE)
    cat("Density from empirical coverage: ", N_emp, '\n')
    
    # density and Histogram of success probability
    h = hist(s_n / t_n, breaks = seq(0, 1.1, 0.01), plot = FALSE)
    h$counts=h$counts/sum(h$counts)
    
    plot(h,
         col = 'lightblue',
         main = m,
         xlab = 'VAF',
         border = NA
    )
    
    lines(x/N_emp, bbin_den, lwd = 2, col = 'darkred')
    
    legend(
      "topright",
      title = 'Beta-Binomial (MLE fit)',
      legend = bquote(mu == .(fit_prob) ~ ',' ~ rho == .(fit_disp) ~ C[symbol("\052")] == .(N_emp)),
      col = c('darkred'),
      bty = 'n',
      pch = 19,
      cex = 1
    )
    
    c(fit_prob, fit_disp)
  },
  m = main
  )
  
  df = Reduce(rbind, df)
  if(is.null(ncol(df))) df = t(data.frame(df))
  
  df = cbind(df, Region = column)
  df = cbind(df, CN = names(data))
  df = data.frame(df, stringsAsFactors = FALSE)
  if(is.null(ncol(df))) df = t(data.frame(df))
  
  df$mu = as.numeric(df$mu)
  df$rho = as.numeric(df$rho)
  
  rownames(df) = NULL
  df
}

BBpval = function(training.params, toTest, significance = 0.05)
{
  panel = lapply(names(toTest), function(w)
    {
      tr = training.params[training.params$CN == w, ]
      tt = toTest[[w]]
  
      pvalues = matrix(1, ncol = nrow(tr), nrow = nrow(tt))
      colnames(pvalues) = tr$Region
      rownames(pvalues) = rownames(tt)

      for(i in 1:nrow(pvalues))
      {
        for(j in 1:nrow(tr))
          pvalues[i, j] = VGAM::dbetabinom(0, tt[i, 'TEST'], prob = tr[j, 'mu'], rho = tr[j, 'rho'], log = FALSE)
      }
      
      cat(cyan('Pre-correction pvalues\n'))
      print(pvalues)
      
      p = cbind(tt, pvalues)
      colnames(p)[1] = 'Num. Reads'
      p
  })
  
  panel = Reduce(rbind, panel)
  
  num_of_tests = Reduce(prod, lapply(toTest, nrow)) * length(unique(training.params$Region))
  cutoff = significance/num_of_tests
  cat(bgRed('\nMHT -- Num of tests:'), num_of_tests, ' --> alpha/ntests = ', cutoff, '\n') 
  
  sign = apply(panel, 1, function(w){
    all(w[2:length(w)] < cutoff)
  })
  panel = cbind(panel, sign = sign)
  
  cat(bgGreen('\nP-values corrected for MHT\n'))
  print(panel)
  
  return(panel)
}

# pplot_fit_power = function(s_n, t_n, tests_t_n,
#                       main, 
#                       psign = 0.05,
#                       purity = c(1, 1),
#                       range_coverage = 1:(max(tests_t_n) + 100))
# {
#   # Multi-layout
#   layout(matrix(c(1, 6, 1, 6, 2, 3, 4, 5, 4, 5), 5, 2, byrow = TRUE))
#   
#   require(VGAM)
#   require(fitdistrplus)
#   require(crayon)
#   
#   purity.training = purity[1]
#   purity.test = purity[2]
#   
#   cat('Input coverage (head) t_n: ')
#   cat(head(t_n))
#   
#   # Correct for purity of the training sample, and saturate if required
#   t_n = ceiling(t_n - (1-purity.training) * t_n)
#   if(any(s_n > t_n)) {
#     w = which(s_n > t_n)
#     s_n[w] = t_n[w]
#     warning('Correction for purity of the training sample leads to s_n > t_n -- setting s_n = t_n.')
#   }
#   
#   cat('\nCorrected for purity of', purity.training, ':')
#   cat(head(t_n), '\n\n')
#   
#   # Fitting a Beta-Binomial to data
#   cat(bgGreen('Beta-Binomial MLE fit: '))
#   fit = Coef(VGAM::vglm(cbind(s_n, t_n - s_n) ~ 1, betabinomial, trace = FALSE))
#   fit_prob = round(fit[1], 3)
#   fit_disp = round(fit[2], 3)
#   cat(cyan('  mu ='), fit['mu'], cyan('rho ='), fit['rho'], '\n')  
#   cat('\t Corrected for purity', purity.training, '\n')
# 
#   # empirical density  
#   N_emp = round(mean(t_n))
#   x = round(seq(0, 1, 0.01) * N_emp)
#   bbin_den = VGAM::dbetabinom(x, N_emp, prob = fit_prob, rho = fit_disp, log = FALSE)
#   cat("\t Generated a density at empirical coverage WES", N_emp, '\n')
#   
#   # Fitting a Binomial to data
#   cat(bgGreen('\nBinomial MLE fit: '))
#   # fitBinom = fitdist(data = s_n, dist = "binom", fix.arg = list(size = N_emp), start = list(prob = mean(s_n/t_n)))
#   fitBinom = fitdist(data = s_n, dist = "binom", fix.arg = list(size = max(t_n)), start = list(prob = mean(s_n/t_n)))
#   
#   # print(fitBinom$estimate['prob'])
#   cat(cyan('  p ='), fitBinom$estimate['prob'], '\n')  
#   
#   bin_den = dbinom(x, N_emp, fitBinom$estimate, log = FALSE)
#   
#   cat('\t Maximum coverage  C=', max(t_n), ' -- used for fit, try to use all values..\n')  
#   cat('\t Init empirical p=', mean(s_n/t_n), '\n')  
#   cat("\t Generated a density at empirical coverage WES", N_emp, '\n')
#   
#   # density and Histogram of success probability
#   h = hist(s_n / t_n, breaks = seq(0, 1.1, 0.01), plot = FALSE)
#   h$counts=h$counts/sum(h$counts)
#   
#   plot(h,
#     col = 'lightblue',
#     main = main,
#     xlab = 'VAF',
#     border = NA
#   )
#   
#   lines(x/N_emp, bbin_den, lwd = 2, col = 'darkred')
#   lines(x/N_emp, bin_den, lwd = 2, lty = 2, col = 'darkblue')
#   
#   legend(
#     "topright",
#     c('data', "Beta-Binomial (MLE fit)", "Binomial (MLE fit)"),
#     col = c('lightblue', 'darkred', "darkblue"),
#     lwd = 2,
#     bty = 'n'
#   )
#   
#   legend(
#     "right",
#     title = 'MLE fit',
#     legend = 
#       sapply(c(
#       bquote(mu == .(fit_prob) ~ ',' ~ rho == .(fit_disp) ~ C[symbol("\052")] == .(N_emp)),
#       bquote(p == .(fitBinom$estimate) ~ C[symbol("\052")] == .(N_emp))
#       ), as.expression),
#     
#     col = c('darkred', 'darkblue'),
#     bty = 'n',
#     pch = 19,
#     cex = 1
#   )
#   
#   # Histogram of trials
#   hist(t_n,
#      col = 'lightgray',
#      freq = F,
#      breaks = 22,
#      border = NA,
#      xlab = '',
#      main = 'Coverage distribution')
#   
#   abline(v = mean(t_n), col = "royalblue", lty = 2)
#   abline(v = median(t_n), col = "red", lty = 2)
#   legend(x = "topright", c("Mean", "Median"), col = c("royalblue", "red"), lwd = c(2, 2), cex = .6, bty = 'n')
#   
#   # Histogram of successes
#   hist(s_n,
#      col = 'orange',
#      freq = F,
#      border = NA,
#      breaks = 22,
#      xlab = '',
#      main = 'Mutants distribution')
#   
#   abline(v = mean(s_n), col = "royalblue", lty = 2)
#   abline(v = median(s_n), col = "red", lty = 2)
#   legend(x = "topright", c("Mean", "Median"), col = c("royalblue", "red"), lwd = c(2, 2), cex = .6, bty = 'n')
#  
#   # P-values for the null model H_0
#   pvalues = VGAM::dbetabinom(0, range_coverage, prob = fit_prob, rho = fit_disp, log = FALSE)
#   names(pvalues) = range_coverage
#   names(range_coverage) = range_coverage
#   min_coverage_sign = range_coverage[min(which(pvalues < psign))]
#   
#   # Bonferroni correction
#   Num_tests = length(tests_t_n)
#   min_coverage_sign_FWER = range_coverage[min(which(pvalues < psign/Num_tests))]
# 
#   if(is.na(min_coverage_sign_FWER)) {
#     cat(bgRed('\nYou do not have minimum coverage to pass this test, all SNVs are rejected\n'))
#     stop('Interrupted.')
#   }
# 
#     
#   plot(pvalues,
#        xlab = 'log(Coverage c)',
#        ylab = bquote(log( H[0] ~ ':' ~ 'BetaBin[ 0 | c,' ~ mu ~','~ rho ~']' )),
#        main = bquote(bold('P-values with FWER:') ~ alpha ~ '=' ~  .(psign) ~ ','~ .(Num_tests) ~ 'tests,'~ pi == .(purity.test)),
#        type = 'l',
#        log = 'xy')
#   
#   abline(h = psign, col = "red", lty = 2)
#   abline(h = psign/Num_tests, col = "orange", lty = 2)
#   
#   abline(v = min_coverage_sign, col = "blue", lty = 2)
#   abline(v = min_coverage_sign_FWER, col = "black", lty = 2) 
#   
#   legend('topright', 
#          legend = 
#            sapply(c(
#              bquote(p[symbol("\052")] ~ '<' ~  .(psign)), 
#              bquote(c[symbol("\052")] ~ '=' ~  .(min_coverage_sign)),
#              bquote(p[symbol("\052")]^{FWER} ~ '<' ~  .(psign/Num_tests)), 
#              bquote(c[symbol("\052")]^{FWER} ~ '=' ~  .(min_coverage_sign_FWER))
#              ), as.expression),
#          col = c('red', 'blue', 'orange', 'black'), 
#          lty = 2, 
#          box.lwd = 0,
#          box.col = "white",
#          bg = "white"
#          ) 
#   
#     cat(bgRed('\nMulitple Hypothesis Testing'), '\n\talpha:', psign, 'with Tests:', Num_tests, ' -->', cyan('FWER p ='), psign/Num_tests, '\n',
#         '\tMinimum coverage for H_0 < p:', min_coverage_sign_FWER, '\n',
#         '\tCoverage Range:', min(range_coverage), '--', max(range_coverage), '\n',
#         '\tPurity correction: ', purity.test
#         )  
# 
#     SNVs = data.frame(NR=tests_t_n)
#     
#     mean_precorrections = mean(tests_t_n)
#     tests_t_n = tests_t_n - (1-purity.test) * tests_t_n
#     mean_postcorrections =  mean(tests_t_n)
#     
#     SNVs = cbind(SNVs, NR.adj = tests_t_n)
#     
#     cat(' - mean testing coverage is now ',
#         mean_postcorrections, ', it was', mean_precorrections, '\n')  
#     
#     rejected = tests_t_n[tests_t_n <= min_coverage_sign_FWER]
#     points(rejected, pvalues[rejected], col = 'darkred', pch = 19)
#   
#     cat(bgMagenta('\nSummary FWER\n')) 
#     cat('\t Rejected SNVs:', red(length(rejected)), '\n') 
#     
#     
#     SNVs = cbind(SNVs, pvalue = pvalues[SNVs$NR.adj])
#     aster = function(p, pmin){
#       if(p>=pmin) return('')
#       if(p < 0.001) return('***')
#       if(p < 0.01) return('**')
#       if(p < 0.05) return('*')
#       return('')
#     }
#     
#     SNVs = cbind(SNVs, sapply(SNVs$pvalue, aster, pmin = psign/Num_tests))
#     colnames(SNVs)[4] = ''
#     
#     accepted = tests_t_n[tests_t_n > min_coverage_sign_FWER]
#     
#     points(accepted, pvalues[accepted], col = 'darkgreen', pch = 19)
#     
#     # cat('*** Accepted FWER :', accepted, '\n')
#     cat('\t Accepted SNVs:', green(length(accepted)), '\n')
#     
#     cat('\nTable\n')
#     print(SNVs)
#     
#     # Histogram of trials for test
#     h = hist(tests_t_n, breaks = 22, plot = FALSE)
#     cuts = cut(h$breaks, c(-Inf, min_coverage_sign_FWER, Inf))
#     
#     # print(h$breaks)
#     # print(c('darkred', 'darkgreen')[cuts])
#     
#     plot(h,
#          col = c('darkred', 'darkgreen')[cuts], #c('darkred', 'darkgreen'),
#          # freq = F,
#          # breaks = 22,
#          border = NA,
#          xlab = '',
#          main = 
#            bquote(bold('Coverage distribution corrected for') ~ pi),
#          sub = paste('Mean: pre-correction', round(mean_precorrections), ' -- post-correction', round(mean_postcorrections))
#     )
#    
#      legend('topright', 
#            legend = c(length(rejected), length(accepted)),
#            col = c('darkred', 'darkgreen'),
#            bty = 'n',
#            pch = 19
#     )
#      
#     
#     ret = list(accepted = pvalues[accepted], rejected = pvalues[rejected])
#     
#     ## Tarone's Z statistics for overdispersion Beta Binomial
#     cat('\nExtra: Bootstrapping Tarone\'s Z statistics for overdispersion.\n')
#     M <- 1000
#     alt_hyp = null_hyp <- vector("numeric", length = M)
#     for (i in 1:M){
#       boostrap_samples = sample(length(t_n), replace = TRUE)
#       
#       t = t_n[boostrap_samples]
#       s = s_n[boostrap_samples]
#       
#       p_hat = sum(s) / sum(t)
#       
#       S = sum( (s - t * p_hat)^2 / (p_hat * (1 - p_hat)) )
#       Z_score = (S - sum(t)) / sqrt(2 * sum(t * (t - 1)))
# 
#       null_hyp[i] <- Z_score
#     }
#     
#     H0 = dnorm(seq(-5,5,0.01), 0, 1)
#     
#     hist(null_hyp, 
#          breaks = 100, 
#          col = 'khaki',
#          border = NA,
#          freq = F,
#          xlim = c(-5, max(null_hyp)),
#          ylim = c(0, max(H0)),
#          xlab = 'Z score',
#          main = paste('Tarone\'s Z statistics\n for overdispersion (M=', M, ')', sep ='') 
#          )
#     
#    lines(seq(-5,5,0.01), H0, lty = 2, col = 'red')
#    legend('topright', 
#           sapply(c(
#             bquote(H[0] ~ ': N(0,1)'), 
#             'Beta-Binomial'
#           ), as.expression),
#           col = c('red', 'khaki'),
#           bty = 'n',
#           pch = 19
#           )
#     
#     
#     
#     return(SNVs)
# }
# 
# loadVCF = function(file)
# {
#   f = read.vcfR(file)
# 
#   f.data = NULL
#   f.data$NV <- extract.gt(f, element='NV', as.numeric=TRUE)
#   f.data$NR <- extract.gt(f, element='NR', as.numeric=TRUE)
#   
#   # check out rownames now -- add chr
#   N = rownames(f.data$NV)
#   for(n in 1:length(N))
#     if(!startsWith(N[n], 'chr')) 
#       N[n] = paste('chr', N[n], sep ='')
#     
#   rownames(f.data$NR) = rownames(f.data$NV) = N
#   
#   # colnames according to panels/ WES
#   if(all(startsWith(colnames(f.data$NV), 'NG')))
#   {
#     colnames(f.data$NV) = colnames(f.data$NR) = unlist(
#       lapply(
#         strsplit(colnames(f.data$NV), '_'), 
#         function(x) return(x[2])))
#     
#     cat('* Colnames edited:', colnames(f.data$NV), '\n')
#   }
#   
#   if(all(startsWith(colnames(f.data$NV), 'SP')))
#   {
#     colnames(f.data$NV) = colnames(f.data$NR) = substr(colnames(f.data$NV), 3, nchar(colnames(f.data$NV)))
#     cat('* Colnames edited:', colnames(f.data$NV), '\n')
#   }
#   
#   f.data$VAF = f.data$NV / f.data$NR
#   
#   if(any(is.na(f.data$VAF))) {
#     warning('There are mutations with NR = 0, which leads to NaN VAFs -- forcing them to be 0.')
#     f.data$VAF[is.na(f.data$VAF)] = 0
#   }
#   
#   
#   return(f.data)
# }
# 
# 
# subset_zeroesM_mutations = function(D, min.coverage = 100)
# {
#   margin = endsWith(colnames(D$NV), 'M')
#   
#   selection = apply(D$NV, 1, function(x) x[margin] == 0)
# 
#   D$NV = D$NV[selection, , drop = FALSE]
#   D$NR = D$NR[selection, , drop = FALSE]
#   
#   selection = apply(D$NR, 1, function(x) x[margin] > min.coverage)
#   
#   D$NV = D$NV[selection, , drop = FALSE]
#   D$NR = D$NR[selection, , drop = FALSE]
#   
#   return(D)
# }
# 
# 
# correctForPurity = function(D, purity, PROCESS)
# {
#   T_purity = purity$TUM
#   M_purity = purity$M
#   
#   # Prepare names  
#   names(T_purity) = paste(PROCESS, names(T_purity), sep = '')
#   names(M_purity) = paste(PROCESS, 'M', sep = '')
#   
#   cat('[correctForPurity] Columns: ', colnames(D), '\nPurity estimates are\n')
#   print(T_purity)
#   print(M_purity)
#   
#   for(c in colnames(D))
#   {
#     if(c %in% names(T_purity)) 
#     {
#       cat('\nCorrecting tumor VAF for ', c, 'with purity', T_purity[c])
#       D[, c] = D[, c] / T_purity[c]
#     }
#     
#     if(c %in% names(M_purity)) 
#     {
#       cat('\nCorrecting margin VAF for ', c, 'with purity', M_purity[c])
#       D[, c] = D[, c] / M_purity[c]
#     }
#   }
#   
#   if(any(D > 1)) 
#   {
#     cat('\n')
#     warning('\nCorrection for purity leads to VAF > 1')
#     D[D > 1] = 1
#   }
#   
#   return(D)
# }
# 
# asmuts  = function(D) {return(rownames(D$NV))}
# assampl = function(D) {return(colnames(D$NV))}
# 
# subset_clonal_mutations = function(D, clonal.cutoff = 0.2, correction = NULL)
# {
#   samples = endsWith(assampl(D), 'M') | endsWith(assampl(D), 'B') | endsWith(assampl(D), 'S')
# 
#   vaf = D$VAF.adj[, !samples]
#   selection = apply(vaf, 1, function(x) all(x > clonal.cutoff))
#      
#   if(!any(is.null(correction)))
#   {
#     # check out correction
#     removed = names(selection)[!selection]
#     corrected = intersect(removed, asmuts(correction))
#     
#     if (length(corrected) > 0) {
#       cat('Correction used for', corrected, '\n')
#       selection[corrected] = TRUE
#     }
#   }
#   
#   D$NV = D$NV[selection, , drop = FALSE]
#   D$NR = D$NR[selection, , drop = FALSE]
#   D$VAF = D$VAF[selection, , drop = FALSE]
#   D$VAF.adj = D$VAF.adj[selection, , drop = FALSE]
#   
#   cat('Found', nrow(D$NV), 'clonal mutations.\n')
#   print(rownames(D$VAF.adj))
#   
#   return(D)
# }
# 
# crossCheck = function(muts, reference, reference.Zeroes)
# {
#   if(length(muts) == 0) return(muts)
#   
#   cat('Cross-checking:', paste(muts, collapse = ', '), '\n\t', paste(asmuts(reference), collapse ='\n\t '), '\n')
#   
#   which.testable = muts %in% asmuts(reference)
#   muts.testable = muts[which.testable]
#   muts.non.testable = muts[!which.testable]
#   
#   if(length(muts.testable) == 0) {
#     cat('Nothing to test, returning.')
#     return(muts.non.testable)
#   }
# 
#   cat('Testing:', paste(muts.testable, collapse = ', '), '\n\t', paste(asmuts(reference.Zeroes), collapse ='\n\t '), '\n')
#   
#   which.found = muts.testable %in% asmuts(reference.Zeroes)
#   muts.found = muts[which.found]
#   
#   cat('\nRejected: ', muts.testable[!which.found], '\n')
#   cat('Passed: ', c(muts.non.testable, muts.found), '\n')
#   
#   return(c(muts.non.testable, muts.found))
# }
# 
# 
# plot_VAF = function(W, 
#                     T1, 
#                     T2,
#                     W_clonal_mutations,
#                     T1_test,
#                     T2_test,
#                     clonal.cutoff,
#                     purity,
#                     PROCESS,
#                     show.SUBCLONAL = FALSE
#                     )
# {
#   require(RColorBrewer)
#   
#   Wvaf = W$VAF.adj
#   T1vaf = T1$VAF.adj
#   T2vaf = T2$VAF.adj
#   
#   ##### 
#   WES_rows = matrix('Subclonal', nrow = nrow(Wvaf), ncol = 1)
#   rownames(WES_rows) = rownames(Wvaf)
#   colnames(WES_rows) = 'Mutation'
#   WES_rows[W_clonal_mutations, ] = 'Clonal'
#   
#   WES_cols = c('darkgoldenrod2', 'darkblue')
#   names(WES_cols) = c('Subclonal', 'Clonal')
# 
#   Wvaf = Wvaf[order(rownames(Wvaf)), ]
#   
#   seqs = c(0, 1e-10, seq(0.01, 1.01, 0.01))
#   sseqs = seqs[seqs <= clonal.cutoff]
#   gseqs = seqs[seqs > clonal.cutoff]
#   colrMuts = c(
#     colorRampPalette(brewer.pal(9, 'Greens'))(length(sseqs)), 
#     colorRampPalette(brewer.pal(9, 'YlOrRd'))(length(gseqs))
#   )
#   
#   # WES_rows = WES_rows[order(WES_rows['Mutation']), , drop = FALSE]
#   # print(WES_rows)
#   # Wvaf = Wvaf[rownames(WES_rows), ]
# 
#   if(!show.SUBCLONAL){
#     WES_rows = WES_rows[WES_rows[, 'Mutation'] != 'Subclonal', , drop = FALSE ]
#     Wvaf = Wvaf[rownames(WES_rows), ]
#   }
#   print(Wvaf)
#     
#   tabWES = pheatmap(
#     Wvaf,
#     color = c('lightgray', colrMuts),
#     breaks = seqs,
#     cellwidth = 10,
#     cellheight = 5,
#     cluster_rows = FALSE,
#     cluster_cols = FALSE,
#     main = paste('WES data'),
#     fontsize_row = 5,
#     annotation_row = as.data.frame(WES_rows), 
#     annotation_colors = list(Mutation = WES_cols),
#     border = NA,
#     silent = TRUE
#   )
#   
#   #####   
#   
#   TES1_rows = matrix('No', nrow = nrow(T1vaf), ncol = 1)
#   rownames(TES1_rows) = rownames(T1vaf)
#   colnames(TES1_rows) = 'Test'
#   TES1_rows[T1_test, ] = 'Yes'
#   
#   T1vaf = T1vaf[order(rownames(T1vaf)), ]
#   
#   TES_cols = c('darkgreen', 'indianred1')
#   names(TES_cols) = c('Yes', 'No')
#   
#   tabT1 = pheatmap(
#     T1vaf,
#     color = c('lightgray', colrMuts),
#     breaks = seqs,
#     cellwidth = 10,
#     cellheight = 5,
#     cluster_rows = FALSE,
#     cluster_cols = FALSE,
#     main = paste('Panel TES1'),
#     fontsize_row = 5,
#     annotation_row = as.data.frame(TES1_rows),
#     annotation_colors = list(Test = TES_cols),
#     border = NA,
#     silent = TRUE
#   )
#   
#   #####   
#   
#   
#   TES2_rows = matrix('No', nrow = nrow(T2vaf), ncol = 1)
#   rownames(TES2_rows) = rownames(T2vaf)
#   colnames(TES2_rows) = 'Test'
#   TES2_rows[T2_test, ] = 'Yes'
#   
#   T2vaf = T2vaf[order(rownames(T2vaf)), ]
#   
#   tabT2 = pheatmap(
#     T2vaf,
#     color = c('lightgray', colrMuts),
#     breaks = seqs,
#     cellwidth = 10,
#     cellheight = 5,
#     cluster_rows = FALSE,
#     cluster_cols = FALSE,
#     main = paste('Panel TES2'),
#     fontsize_row = 5,
#     annotation_row = as.data.frame(TES2_rows),
#     annotation_colors = list(Test = TES_cols),
#     border = NA,
#     silent = TRUE
#   )
# 
#   grid.arrange(tabWES$gtable, tabT1$gtable, tabT2$gtable, ncol = 3)
# }
# 
# 
# # plot_training_set = function(D, file = 'plot_training_set.pdf')
# # {
# #   par(mfrow = c(2,1))
# #   
# #   D.cl = subset_clonal_mutations(D, type = 'clonal')
# #   D.scl = subset_clonal_mutations(D, type = 'subclonal')
# #   
# #   vaf.cl = D.cl$NV / D.cl$NR
# #   vaf.scl = D.scl$NV / D.scl$NR
# #   
# #   hist(vaf.cl, breaks = seq(0, 1.01, 0.01), main = 'VAF of Clonal mutations', border = NA, col = 'lightblue', xlab = 'VAF')
# #   hist(vaf.scl, breaks = seq(0, 1.01, 0.01), main = 'VAF of Subclonal mutations', border = NA, col = 'lightblue', xlab = 'VAF')
# #   dev.copy2pdf(file = file)
# # }
# 
# 
# batch_tests = function(WES.clonal, panel.zeroes, panel.testable, margin.sample, purity, psign, panel)
# {
#   require(crayon)
#   
#   if(length(panel.testable) == 0) stop('Nothing to test, aborting.')
#   
#   samples = endsWith(assampl(WES.clonal), 'M') | endsWith(assampl(WES.clonal), 'B') | endsWith(assampl(WES.clonal), 'S')
#   WES.clonal$NV = WES.clonal$NV[, !samples]
#   WES.clonal$NR = WES.clonal$NR[, !samples]
#   
# 
#     cat(bgGreen('Tests Matrix\n'))
#   # print(pvalues)
#   
#   cat('Purity', purity$TUM, ' and margin is ', purity$M, '\n')
#   # T_purity = purity$TUM
#   # names(T_purity) = assampl(WES.clonal)
#   cat('Tumour purity is unused here!\n')
#   
#   M_purity = purity$M
#   result = NULL
# 
#   for(t in 1:length(assampl(WES.clonal)))
#   {
#     tested.sample = assampl(WES.clonal)[t]
#     
#     cat('\n\n************** TRAINING ', tested.sample, '\n\n')
#     
#     quartz(height = 10)
#     
#     toTest = panel.zeroes$NR[panel.testable, margin.sample]
#     names(toTest) = panel.testable
#     
#     ex = pplot_fit_power(
#       s_n = WES.clonal$NV[, tested.sample],
#       t_n = WES.clonal$NR[, tested.sample],
#       tests_t_n = toTest,
#       main = tested.sample,
#       psign = psign,
#       purity = c(1, M_purity),
#       range_coverage = 1:20000)
#     
#     # remove last column
#     ex = ex[, 1:(ncol(ex)-1)]
#     colnames(ex)[ncol(ex)] = paste(tested.sample, '.pvalue', sep = '')
#     
#     if(t == 1) result = ex
#     else{
#       result = cbind(
#         result,
#         ex[, paste(tested.sample, '.pvalue', sep = ''), drop = FALSE]
#         )
#     }
#     
#     dev.copy2pdf(file = paste(panel, tested.sample, '.pdf', sep =''))
#     dev.off()
#     
#     readline('Hit for next test...')
#   }
#   
#   print(result)
# 
#   ## Very hard correction
#   numTests = nrow(result) * length(assampl(WES.clonal))
#   psignFWER = psign/numTests
#   
#   cat(bgRed('\nOverall Bonferroni FWER correction for', numTests,' number of tests\n'))
#   cat(cyan('Significance at level', psign, 'obtained with p <', psignFWER), '\n')
#   
#   means = apply(
#     result,
#     1,
#     function(x) mean(x[3:length(x)]))
# 
#   asts = apply(
#     result,
#     1,
#     function(x) {
#       ps = x[3:length(x)]
#       aps = mean(ps)
#       if(all(ps < psignFWER)) 
#       {
#         if(aps < 0.001) return(paste('***'))
#         if(aps < 0.01) return(paste('**'))
#         if(aps < 0.05) return(paste('*'))
#       }
#       else return(paste('reject'))
#     })
#   
#   result$means = means
#   result$sign = asts
#   
#   
#   result = result[order(result$means), ]
#   
#   return(result)
# }
