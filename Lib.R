pplot_fit_power = function(s_n, t_n, tests_t_n,
                      main, 
                      psign = 0.05,
                      purity = c(1, 1),
                      range_coverage = 1:(max(tests_t_n) + 100))
{
  # Multi-layout
  layout(matrix(c(1, 6, 1, 6, 2, 3, 4, 5, 4, 5), 5, 2, byrow = TRUE))
  
  require(VGAM)
  require(fitdistrplus)
  require(crayon)
  
  purity.training = purity[1]
  purity.test = purity[2]
  
  cat('Input coverage (head) t_n: ')
  cat(head(t_n))
  
  # Correct for purity of the training sample, and saturate if required
  t_n = ceiling(t_n - (1-purity.training) * t_n)
  if(any(s_n > t_n)) {
    w = which(s_n > t_n)
    s_n[w] = t_n[w]
    warning('Correction for purity of the training sample leads to s_n > t_n -- setting s_n = t_n.')
  }
  
  cat('\nCorrected for purity of', purity.training, ':')
  cat(head(t_n), '\n\n')
  
  # Fitting a Beta-Binomial to data
  cat(bgGreen('Beta-Binomial MLE fit: '))
  fit = Coef(VGAM::vglm(cbind(s_n, t_n - s_n) ~ 1, betabinomial, trace = FALSE))
  fit_prob = round(fit[1], 3)
  fit_disp = round(fit[2], 3)
  cat(cyan('  mu ='), fit['mu'], cyan('rho ='), fit['rho'], '\n')  
  cat('\t Corrected for purity', purity.training, '\n')

  # empirical density  
  N_emp = round(mean(t_n))
  x = round(seq(0, 1, 0.01) * N_emp)
  bbin_den = VGAM::dbetabinom(x, N_emp, prob = fit_prob, rho = fit_disp, log = FALSE)
  cat("\t Generated a density at empirical coverage WES", N_emp, '\n')
  
  # Fitting a Binomial to data
  cat(bgGreen('\nBinomial MLE fit: '))
  # fitBinom = fitdist(data = s_n, dist = "binom", fix.arg = list(size = N_emp), start = list(prob = mean(s_n/t_n)))
  fitBinom = fitdist(data = s_n, dist = "binom", fix.arg = list(size = max(t_n)), start = list(prob = mean(s_n/t_n)))
  
  # print(fitBinom$estimate['prob'])
  cat(cyan('  p ='), fitBinom$estimate['prob'], '\n')  
  
  bin_den = dbinom(x, N_emp, fitBinom$estimate, log = FALSE)
  
  cat('\t Maximum coverage  C=', max(t_n), ' -- used for fit, try to use all values..\n')  
  cat('\t Init empirical p=', mean(s_n/t_n), '\n')  
  cat("\t Generated a density at empirical coverage WES", N_emp, '\n')
  
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
 
  # P-values for the null model H_0
  pvalues = VGAM::dbetabinom(0, range_coverage, prob = fit_prob, rho = fit_disp, log = FALSE)
  names(pvalues) = range_coverage
  names(range_coverage) = range_coverage
  min_coverage_sign = range_coverage[min(which(pvalues < psign))]
  
  # Bonferroni correction
  Num_tests = length(tests_t_n)
  min_coverage_sign_FWER = range_coverage[min(which(pvalues < psign/Num_tests))]

  if(is.na(min_coverage_sign_FWER)) {
    cat(bgRed('\nYou do not have minimum coverage to pass this test, all SNVs are rejected\n'))
    stop('Interrupted.')
  }

    
  plot(pvalues,
       xlab = 'log(Coverage c)',
       ylab = bquote(log( H[0] ~ ':' ~ 'BetaBin[ 0 | c,' ~ mu ~','~ rho ~']' )),
       main = bquote(bold('P-values with FWER:') ~ alpha ~ '=' ~  .(psign) ~ ','~ .(Num_tests) ~ 'tests,'~ pi == .(purity.test)),
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
  
    cat(bgRed('\nMulitple Hypothesis Testing'), '\n\talpha:', psign, 'with Tests:', Num_tests, ' -->', cyan('FWER p ='), psign/Num_tests, '\n',
        '\tMinimum coverage for H_0 < p:', min_coverage_sign_FWER, '\n',
        '\tCoverage Range:', min(range_coverage), '--', max(range_coverage), '\n',
        '\tPurity correction: ', purity.test
        )  

    SNVs = data.frame(NR=tests_t_n)
    
    mean_precorrections = mean(tests_t_n)
    tests_t_n = tests_t_n - (1-purity.test) * tests_t_n
    mean_postcorrections =  mean(tests_t_n)
    
    SNVs = cbind(SNVs, NR.adj = tests_t_n)
    
    cat(' - mean testing coverage is now ',
        mean_postcorrections, ', it was', mean_precorrections, '\n')  
    
    rejected = tests_t_n[tests_t_n <= min_coverage_sign_FWER]
    points(rejected, pvalues[rejected], col = 'darkred', pch = 19)
  
    cat(bgMagenta('\nSummary FWER\n')) 
    cat('\t Rejected SNVs:', red(length(rejected)), '\n') 
    
    
    SNVs = cbind(SNVs, pvalue = pvalues[SNVs$NR.adj])
    aster = function(p, pmin){
      if(p>=pmin) return('')
      if(p < 0.001) return('***')
      if(p < 0.01) return('**')
      if(p < 0.05) return('*')
      return('')
    }
    
    SNVs = cbind(SNVs, sapply(SNVs$pvalue, aster, pmin = psign/Num_tests))
    colnames(SNVs)[4] = ''
    
    accepted = tests_t_n[tests_t_n > min_coverage_sign_FWER]
    
    points(accepted, pvalues[accepted], col = 'darkgreen', pch = 19)
    
    # cat('*** Accepted FWER :', accepted, '\n')
    cat('\t Accepted SNVs:', green(length(accepted)), '\n')
    
    cat('\nTable\n')
    print(SNVs)
    
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
    cat('\nExtra: Bootstrapping Tarone\'s Z statistics for overdispersion.\n')
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
    
    
    
    return(SNVs)
}

loadVCF = function(file)
{
  f = read.vcfR(file)

  f.data = NULL
  f.data$NV <- extract.gt(f, element='NV', as.numeric=TRUE)
  f.data$NR <- extract.gt(f, element='NR', as.numeric=TRUE)
  
  # check out rownames now -- add chr
  N = rownames(f.data$NV)
  for(n in 1:length(N))
    if(!startsWith(N[n], 'chr')) 
      N[n] = paste('chr', N[n], sep ='')
    
  rownames(f.data$NR) = rownames(f.data$NV) = N
  
  # colnames according to panels/ WES
  if(all(startsWith(colnames(f.data$NV), 'NG')))
  {
    colnames(f.data$NV) = colnames(f.data$NR) = unlist(
      lapply(
        strsplit(colnames(f.data$NV), '_'), 
        function(x) return(x[2])))
    
    cat('* Colnames edited:', colnames(f.data$NV), '\n')
  }
  
  if(all(startsWith(colnames(f.data$NV), 'SP')))
  {
    colnames(f.data$NV) = colnames(f.data$NR) = substr(colnames(f.data$NV), 3, nchar(colnames(f.data$NV)))
    cat('* Colnames edited:', colnames(f.data$NV), '\n')
  }
  
  return(f.data)
}

# table_clonal_in_T = function(W, T1, T2)
# {
#   mutations = unique(
#     c(rownames(W$NR),     
#     rownames(T1$NR),     
#     rownames(T2$NR))
#     )
#   
#   M = matrix(NA, nrow = length(mutations), ncol = 2)
#   rownames(M) = mutations
#   colnames(M) = c('WES', 'TES1')
#   
#   W = subset_clonal_mutations(W, type = 'clonal')
#   M[rownames(W$NR), 'WES'] = 1
# 
#   T1 = subset_clonal_mutations(T1, type = 'clonal')
#   M[rownames(T1$NR), 'TES1'] = 1
#   
#   w = which(is.na(M[, 'WES']) &  M[, 'TES1'] == 1)
#   M[w, 'WES'] = 1
#   
#   
#   M = M[order(rowSums(M, na.rm = T), decreasing = TRUE), ]
#   M = M[rowSums(M, na.rm = T) > 0, ]
#   
#   return(as.data.frame(M))
# }


subset_zeroesM_mutations = function(D, min.coverage = 100)
{
  margin = endsWith(colnames(D$NV), 'M')
  
  selection = apply(D$NV, 1, function(x) x[margin] == 0)

  D$NV = D$NV[selection, , drop = FALSE]
  D$NR = D$NR[selection, , drop = FALSE]
  
  selection = apply(D$NR, 1, function(x) x[margin] > min.coverage)
  
  D$NV = D$NV[selection, , drop = FALSE]
  D$NR = D$NR[selection, , drop = FALSE]
  
  return(D)
}


asmuts  = function(D) {return(rownames(D$NV))}
assampl = function(D) {return(colnames(D$NV))}

subset_clonal_mutations = function(D, type, correction = NULL)
{
  samples = endsWith(assampl(D), 'M') | endsWith(assampl(D), 'B') | endsWith(assampl(D), 'S')
  D$NV = D$NV[, !samples, drop = FALSE]
  D$NR = D$NR[, !samples, drop = FALSE]
  
  selection = NULL
  if(type == 'clonal') {
   
     selection = apply(D$NV, 1, function(x) all(x > 0))
     
     if(!any(is.null(correction)))
     {
       # check out correction
       removed = names(selection)[!selection]
       corrected = intersect(removed, asmuts(correction))
       
       if(length(corrected) > 0) {
         cat('Correction used for', corrected)
         selection[corrected] = TRUE
       }
      }
     
     # print(selection[asmuts(correction)])
  }
  else
    selection = apply(D$NV, 1, function(x) any(x == 0))
    
  D$NV = D$NV[selection, , drop = FALSE]
  D$NR = D$NR[selection, , drop = FALSE]
  
  
  return(D)
}


plot_VAF = function(W, T1, T2)
{
  require(RColorBrewer)
  
  Wvaf = W$NV / W$NR
  T1vaf = T1$NV / T1$NR
  T2vaf = T2$NV / T2$NR
  
  if(any(is.na(Wvaf)) || any(is.na(T1vaf)) || any(is.na(T2vaf))) {
    warning('There are mutations with NR = 0, which leads to NaN VAFs -- forcing them to be 0.')
    Wvaf[is.na(Wvaf)] = 0
    T1vaf[is.na(T1vaf)] = 0
    T2vaf[is.na(T2vaf)] = 0
  }
  
  # annotation_cols = matrix('Tumour', nrow = ncol(Wvaf$NV), ncol = 1)
  # rownames(annotation_cols) = colnames(Wvaf$NV)
  # colnames(annotation_cols) = 'Sample'
  # 
  # margin.sample = endsWith(colnames(Wvaf$NV), 'M') 
  # blood.sample = endsWith(colnames(Wvaf$NV), 'B') 
  # stem.sample = endsWith(colnames(Wvaf$NV), 'S') 
  # 
  # annotation_cols[margin.sample, ] = 'Margin'
  # annotation_cols[blood.sample, ] = 'Blood'
  # annotation_cols[stem.sample, ] = 'Stem'
  # 
  # annotation_rows = matrix('NO', nrow = nrow(Wvaf$NV), ncol = 2)
  # rownames(annotation_rows) = rownames(Wvaf$NV)
  # colnames(annotation_rows) = c('Mutation', mut.annot.label)
  # 
  # clonal = subset_clonal_mutations(D, type = 'clonal')
  # clonal = rownames(clonal$NV)
  # 
  # annotation_rows[clonal, 'Mutation'] = 'clonal in T'
  # if(length(mut.annot) > 0) annotation_rows[mut.annot, mut.annot.label] = 'YES'
  
  seqs = c(0, 1e-10, seq(0.01, 1.01, 0.01))
  
  Wvaf = vaf[order(rownames(vaf)), ]
  
  tab = pheatmap(
    vaf,
    color = c('lightgray', colorRampPalette(brewer.pal(9, 'Blues'))(length(seqs) - 2)),
    breaks = seqs,
    cellwidth = 10,
    cellheight = 5,
    cluster_rows = FALSE,
    main = main,
    fontsize_row = 5,
    # annotation_row = as.data.frame(annotation_rows),
    # annotation_col = as.data.frame(annotation_cols),
    border = NA,
    silent = TRUE
  )
  
  return(tab$gtable)
}

# plot_VAF = function(D, main = 'VAF plot', mut.annot, mut.annot.label)
# {
#   require(RColorBrewer)
#   
#   vaf = D$NV / D$NR
#   
#   if(any(is.na(vaf))) {
#     warning('There must be some mutations with NR = 0, which leads to NaN VAFs -- forcing them to be 0.')
#     vaf[is.na(vaf)] = 0
#   }
#   
#   annotation_cols = matrix('Tumour', nrow = ncol(D$NV), ncol = 1)
#   rownames(annotation_cols) = colnames(D$NV)
#   colnames(annotation_cols) = 'Sample'
#   
#   margin.sample = endsWith(colnames(D$NV), 'M') 
#   blood.sample = endsWith(colnames(D$NV), 'B') 
#   stem.sample = endsWith(colnames(D$NV), 'S') 
#   
#   annotation_cols[margin.sample, ] = 'Margin'
#   annotation_cols[blood.sample, ] = 'Blood'
#   annotation_cols[stem.sample, ] = 'Stem'
#   
#   annotation_rows = matrix('NO', nrow = nrow(D$NV), ncol = 2)
#   rownames(annotation_rows) = rownames(D$NV)
#   colnames(annotation_rows) = c('Mutation', mut.annot.label)
# 
#   clonal = subset_clonal_mutations(D, type = 'clonal')
#   clonal = rownames(clonal$NV)
#   
#   annotation_rows[clonal, 'Mutation'] = 'clonal in T'
#   if(length(mut.annot) > 0) annotation_rows[mut.annot, mut.annot.label] = 'YES'
#   
#   seqs = c(0, 1e-10, seq(0.01, 1.01, 0.01))
# 
#   vaf = vaf[order(rownames(vaf)), ]
#   tab = pheatmap(
#     vaf,
#     color = c('lightgray', colorRampPalette(brewer.pal(9, 'Blues'))(length(seqs) - 2)),
#     breaks = seqs,
#     cellwidth = 10,
#     cellheight = 5,
#     cluster_rows = FALSE,
#     main = main,
#     fontsize_row = 5,
#     annotation_row = as.data.frame(annotation_rows),
#     annotation_col = as.data.frame(annotation_cols),
#     border = NA,
#     silent = TRUE
#     )
#   
#   return(tab$gtable)
# }

plot_training_set = function(D, file = 'plot_training_set.pdf')
{
  par(mfrow = c(2,1))
  
  D.cl = subset_clonal_mutations(D, type = 'clonal')
  D.scl = subset_clonal_mutations(D, type = 'subclonal')
  
  vaf.cl = D.cl$NV / D.cl$NR
  vaf.scl = D.scl$NV / D.scl$NR
  
  hist(vaf.cl, breaks = seq(0, 1.01, 0.01), main = 'VAF of Clonal mutations', border = NA, col = 'lightblue', xlab = 'VAF')
  hist(vaf.scl, breaks = seq(0, 1.01, 0.01), main = 'VAF of Subclonal mutations', border = NA, col = 'lightblue', xlab = 'VAF')
  dev.copy2pdf(file = file)
}

# prepare_training_test_set = function(
#   margin.sample,
#   WES,
#   TES1,
#   TES2
# )
# {
#   just.clonal = TRUE
#   
#   primary.samples = setdiff(colnames(TES1$NV), margin.sample)
#   primary.samples = primary.samples[!endsWith(primary.samples, 'B')]
#   primary.samples = primary.samples[!endsWith(primary.samples, 'S')]
#   
#   cat(bgGreen('DATA\n'))
#   cat('* Margin :', margin.sample, '\n')
#   cat('* Samples:', primary.samples, '\n')
# 
#   if(just.clonal)
#   {
#     nM = nrow(WES$NR)
#     WES = subset_clonal_mutations(WES, type = 'clonal')
#     nMc = nrow(WES$NR)
#     cat('* Subsetting WES data to clonal mutations in T: ', nMc, 'out of', nM, '\n')
#   }
#   
#   # First, take SNVs in the targeted panels that have 0 variant reads
#   zeros.TES1 = rownames(TES1$NV)[which(TES1$NV[, margin.sample] == 0)]
#   zeros.TES2 = rownames(TES2$NV)[which(TES2$NV[, margin.sample] == 0)]
#   
#   cat(bgGreen('\nFILTER #1'), cyan('SNVs in targeted panel that have 0 variant reads\n'))
#   cat('\t TES1 : ', (zeros.TES1), '\n')
#   cat('\t TES2 : ', (zeros.TES2), '\n')
#   cat(bgRed('Counts'), 'TES1', length(zeros.TES1), 'TES2', length(zeros.TES2), '\n')
#   
#   # cat('VAF status in WES\n')
#   # print(WES$NV[c(zeros.TES1, zeros.TES2), ])
#   # 
#   
#   # Second, keep only targeted SNVs that are listed in the primary (if they are not, it means that 
#   # they had no mutations, and so we do not need them)
#   if(length(zeros.TES1) > 0) zeros.TES1 = zeros.TES1[zeros.TES1 %in% rownames(WES$NV)]
#   if(length(zeros.TES2) > 0) zeros.TES2 = zeros.TES2[zeros.TES2 %in% rownames(WES$NV)]
#   
#   cat(bgGreen('\nFILTER #2'), cyan('SNVs in targeted panel that are clonal in WES samples\n'))
#   cat('\t TES1 (head): ', head(zeros.TES1), '\n')
#   cat('\t TES2 (head): ', head(zeros.TES2), '\n')
#   cat(bgRed('Counts'), 'TES1', length(zeros.TES1), 'TES2', length(zeros.TES2), '\n')
#   
#   # Third, keep only targeted SNVs that have at least a mutant read in the primary. Most of them
#   # should have, but some might have mutatns in the margin from the primary. This contradicts the fact
#   # that they have no variants in the targeted panel, which has much higher coverage. We trust more the
#   # targeted panel, and thus remove such variants considering those reads false positives.
#   if(length(zeros.TES1) > 0) zeros.TES1 = zeros.TES1[rowSums(WES$NV[zeros.TES1, primary.samples, drop = FALSE]) > 0]
#   if(length(zeros.TES2) > 0) zeros.TES2 = zeros.TES2[rowSums(WES$NV[zeros.TES2, primary.samples, drop = FALSE]) > 0]
#   
#   cat(bgGreen('\nFILTER #3'), cyan('targeted SNVs that have at least a mutant read in the primary\n'))
#   cat('\t TES1 (head): ', head(zeros.TES1), '\n')
#   cat('\t TES2 (head): ', head(zeros.TES2), '\n')
#   cat(bgRed('Counts'), 'TES1', length(zeros.TES1), 'TES2', length(zeros.TES2), '\n')
#   
#   # There are some weird cases as well, in which the coverage at the targeted panel is below say 100x
#   # Since we correct for MHT and we are strict with FWER, we better remove these entries which will
#   # never be significant
#   if(length(zeros.TES1) > 0) zeros.TES1 = zeros.TES1[TES1$NR[zeros.TES1, margin.sample, drop = FALSE] > 100]
#   if(length(zeros.TES2) > 0) zeros.TES2 = zeros.TES2[TES2$NR[zeros.TES2, margin.sample, drop = FALSE] > 100]
#   
#   cat(bgGreen('\nFILTER #4'), cyan('targeted SNVs that have at minimum coverage of 100 reads\n'))
#   cat('\t TES1 (head): ', head(zeros.TES1), '\n')
#   cat('\t TES2 (head): ', head(zeros.TES2), '\n')
#   cat(bgRed('Counts'), 'TES1', length(zeros.TES1), 'TES2', length(zeros.TES2), '\n')
#   
#   # For the test, we need first a training model. We train on SNVs from the tested.sample that are not
#   # the ones that we are going to need for the pvalue. These variants are
#   training.variants = rownames(WES$NV)
#   training.variants = training.variants[!(training.variants %in% c(zeros.TES1, zeros.TES2))]
#   
#   training.WES.data = WES
#   training.WES.data$NV = training.WES.data$NV[training.variants, ]
#   training.WES.data$NR = training.WES.data$NR[training.variants, ]
#   
#   return(
#     list(
#       training = training.WES.data,
#       test.TES1 = zeros.TES1,
#       test.TES2 = zeros.TES2
#     )
#   )
#   
# }

batch_tests = function(WES.clonal, panel.zeroes, panel.testable, margin.sample, purity, psign = 0.05)
{
  require(crayon)
  
  if(length(panel.testable) == 0) stop('Nothing to test, aborting.')
  
  pvalues = matrix(NA, ncol = length(assampl(WES.clonal)), nrow = length(panel.testable))
  colnames(pvalues) = assampl(WES.clonal)
  rownames(pvalues) = panel.testable
  
  cat(bgGreen('Tests Matrix\n'))
  print(pvalues)
  
  T_purity = purity$TUM
  names(T_purity) = assampl(WES.clonal)
  
  M_purity = purity$M
  
  
  result = NULL

  for(tested.sample in assampl(WES.clonal))
  {
    cat('\n\n************** TRAINING ', tested.sample, '\n\n')
    
    quartz(height = 10)
    
    toTest = panel.zeroes$NR[panel.testable, margin.sample]
    names(toTest) = panel.testable

    ex = pplot_fit_power(
      s_n = WES.clonal$NV[, tested.sample],
      t_n = WES.clonal$NR[, tested.sample],
      tests_t_n = toTest,
      main = tested.sample,
      psign = psign,
      purity = c(1, M_purity),
      range_coverage = 1:20000)
    
    ex = cbind(ex, region = tested.sample, variant = rownames(ex))
    result = rbind(result, ex)
    
    dev.copy2pdf(file = paste(tested.sample, '.pdf', sep =''))
    dev.off()
    
    # readline('Hit for next test...')
  }
  
  ## Very hard correction
  numTests = nrow(result)
  psignFWER = psign/numTests
  
  cat(bgRed('\nOverall Bonferroni FWER correction for', numTests,' number of tests\n'))
  cat(cyan('Significance at level', psign, 'obtained with p <', psignFWER), '\n')
  
  result = cbind(
    result, 
    outcome = 'ACCEPT_H0')
  
  result$outcome = result$pvalue < psignFWER
  
  # 
  result = result[, c('variant', 'region', 'NR', 'NR.adj', 'pvalue', 'outcome', 'Var.4')]
  rownames(result) = NULL
  colnames(result)[c(6, 7)] = c(paste('p <', psignFWER), '')
  
  return(result)
}
