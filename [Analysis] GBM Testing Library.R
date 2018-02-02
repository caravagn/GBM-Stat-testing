########## Author: Giulio Caravagna, ICR. <giulio.caravagna@icr.ac.uk> or <gcaravagn@gmail.com>
##########


MIN_TRAINING_SIZE = 10
CLONALITY_CUTOFF = 0.8
WORST_CASE_SCENARIO_PURITY = 0.01


##########


primary = function(x) 
{
  x[, grepl('T', colnames(x)), drop = FALSE]
}

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

correctReadCounts = function(reads.table, purity)
{
  tumourSpecificReads = function(coverage, Ct = 2, purity = 1) {
    
    # Ct = tumour copy number
    # purity = tumour purity
    # coverage = total coverage at the locus
    
    # How many reads are from the tumour?
    tumour.reads = (Ct*purity) / ((Ct*purity) + (2*(1 - purity))) * coverage
    
    return(tumour.reads)
  }
  
  if(nrow(reads.table) == 0 | ncol(reads.table) == 1) return(reads.table)
  
  cna.status = reads.table[, ncol(reads.table)]
  reads = as.matrix(reads.table[, -ncol(reads.table), drop = FALSE])
  
  for(i in 1:nrow(reads))
  {
    for(j in 1:ncol(reads))
      reads[i, j] =  ceiling(tumourSpecificReads(reads[i, j], cna.status[i], purity[j, ]))
  }
  
  return(data.frame(reads, CNA = cna.status))
}

plot_density = function(s_n, t_n, cn, sample, fit_prob, fit_disp) 
{
  # empirical density  
  N_emp = round(mean(t_n))
  x = round(seq(0, 1, 0.01) * N_emp)
  
  bbin_den = VGAM::dbetabinom(
    x, 
    N_emp, 
    prob = fit_prob, 
    rho = fit_disp, 
    log = FALSE)
  
  # density and Histogram of success probability
  h = hist(s_n / t_n, breaks = seq(0, 1.1, 0.01), plot = FALSE)
  h$counts = h$counts / sum(h$counts)
  
  title = paste('BBMLE for ', sample, ' with CN ', cn, ' (n = ', length(s_n),')', sep = '')
  
  plot(h,
       col = 'lightblue',
       main = title,
       xlab = 'VAF',
       ylab = NA,
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
  
}


BBMLE = function(NR, NV, samples, patient)
{
  require(VGAM)
  require(fitdistrplus)
  require(crayon)
  
  cat(bgRed('BBMLE fit'),  '\n')
  cat('\t', green('CNA:'), names(NR), '\n')
  cat('\t', blue('samples:'), samples, '\n')
  
  
  df = lapply(1:length(NR), function(w) 
  {
    ret.fit = NULL
    CN = names(NR)[w]
    
    for(s in samples){
      
      # Data for fit
      s_n = NV[[w]][, s]
      t_n = NR[[w]][, s]
      
      print.correction = any(t_n < s_n)
      if(print.correction) t_n[t_n < s_n] = s_n[t_n < s_n]
      
      # MLE BetaBin    
      cat(bgGreen('\n BetaBinomial MLE fit ') , s, ' with CN =', CN, ' : ')
      fit = Coef(VGAM::vglm(cbind(s_n, t_n - s_n) ~ 1, betabinomial, trace = FALSE))
      fit_prob = round(fit[1], 3)
      fit_disp = round(fit[2], 3)
      cat(cyan('  mu ='), fit['mu'], cyan('rho ='), fit['rho'])  
      
      if(print.correction) cat(red(" [Correction k < n]"))
      
      plot_density(s_n, t_n, CN, s, fit_prob, fit_disp) 
      dev.copy2pdf(file = paste(patient, '-CN_', CN, '-sample_', s, '.pdf', sep = ''))
      
      ret.fit = rbind(ret.fit, c(fit_prob, fit_disp, s, CN)) 
    }
    
    colnames(ret.fit) = c('mu', 'rho', 'sample', 'CN')
    ret.fit
  })
  
  pdfs = list.files()
  pdfs = pdfs[endsWith(pdfs, '.pdf')]
  pdfs = pdfs[startsWith(pdfs, patient)]
  jamPDF(pdfs, out.file = paste('BBMLE-', patient, '.pdf', sep = ''))
  
  df = Reduce(rbind, df)
  df = data.frame(df, stringsAsFactors = FALSE)
  
  df$mu = as.numeric(df$mu)
  df$rho = as.numeric(df$rho)
  
  df
}




