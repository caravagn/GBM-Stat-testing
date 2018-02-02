########## Author: Giulio Caravagna, ICR. <giulio.caravagna@icr.ac.uk> or <gcaravagn@gmail.com>
##########


getData = function(patient) {
  load(paste('[Data] TES_1/NR-', patient, '.RData', sep = ''), verbose = T)
  load(paste('[Data] TES_2/NR-', patient, '.RData', sep = ''), verbose = T)
  load(paste('[Data] WES_PASS/BBMLE-', patient, '.RData', sep = ''), verbose = T)

  # cross-checking -- removing those rejcted by the other panel
  TES1$NR = TES1$NR[!rownames(TES1$NR) %in%  TES2$Rejected, , drop = FALSE]
  TES2$NR = TES2$NR[!rownames(TES2$NR) %in%  TES1$Rejected, , drop = FALSE]
  
  library(crayon)
  
  T1isZero = nrow(TES1$NR) == 0
  T2isZero = nrow(TES2$NR) == 0
  
  # If it exists the training set
  T1hasTraining = TES1$NR[, 'CNA'] %in% fit$CN
  T2hasTraining = TES2$NR[, 'CNA'] %in% fit$CN
  
  if(!T1isZero && any(!T1hasTraining))
  {
    TES1$NR =TES1$NR[T1hasTraining, , drop = FALSE]
    T1isZero = nrow(TES1$NR) == 0
    
    cat(red('[TES 1] -- There are SNVs, but no trainibg set of their CN status\n'))
  }
  
  if(!T2isZero && any(!T2hasTraining))
  {
    TES2$NR =TES2$NR[T2hasTraining, , drop = FALSE]
    T2isZero = nrow(TES2$NR) == 0
    
    cat(red('[TES 2] -- There are SNVs, but missing training set for their CN status\n'))
  }
  
  if(
    (!T1isZero && rownames(TES1$NR) %in% rownames(TES2$NR)) 
    | (!T2isZero && rownames(TES2$NR) %in% rownames(TES1$NR))) 
  {
    cat(red('SNVs appear in two panels -- wanna merge them\n'))
    print(TES1$NR)
    print(TES2$NR)
  }
  
  if(T1isZero) cat(red('[TES 1] -- no SNV to test here\n'))
  if(T2isZero) cat(red('[TES 2] -- no SNV to test here\n'))
  
  
  return(list(TES1 = TES1, TES2 = TES2, fit = fit))
}

BBpval = function(TES, fit)
{
  SNVs = rownames(TES)
  toTest = TES[, colnames(TES) != 'CNA', drop = FALSE]
  samples = colnames(toTest)
  training.samples = unique(fit$sample)
  
  cat(bgRed(' [Beta-Binomial H0 for true positive testing] '), green(paste(SNVs, collapse = ',')), '\n')
  
  
  ret = lapply(SNVs,
         function(w){
           
           CN.fit = unique(fit[fit$CN == TES[w, 'CNA'], 'CN'])
           parameters.fit = fit[fit$CN == CN.fit, ]
           
           # All combinations of training/ test
           tests = expand.grid(SNV = w, Panel = samples, Training = training.samples, stringsAsFactors = FALSE)
           tests = cbind(tests, data.frame(mu = 0, rho = 0, coverage = 0))
           
           # copy parameters
           for(t in 1:nrow(tests))
           {
             tests[t, 'mu'] = parameters.fit[parameters.fit$sample == tests[t, 'Training'], 'mu']
             tests[t, 'rho'] = parameters.fit[parameters.fit$sample == tests[t, 'Training'], 'rho']
             tests[t, 'coverage'] = toTest[w, tests[t, 'Panel']]
           }
           
           tests = cbind(tests, pvalue = 1)
           
           tests$pvalue = apply(tests, 1, function(x) {
             
             # H0 -- no reads given the expected Beta-Binomial distribution
             VGAM::dbetabinom(
               0, 
               as.numeric(x['coverage']), 
               prob = as.numeric(x['mu']), 
               rho = as.numeric(x['rho']), 
               log = FALSE)
           })
           
           tests
         })
  
  names(ret) = SNVs
  ret
}  

MHT = function(pvalues, significance = 0.05)
{
  numTests = lapply(pvalues, nrow)
  numTests = sum(unlist(numTests))
  
  cat(bgRed(' [MHT Correction via Bonferroni] '), red('m =', numTests), green('alpha =', significance), red('alpha/m =', significance/numTests), '\n')
  
  pvalues = lapply(pvalues, function(w) {
    w$sign = w$pvalue < (significance/numTests)
    w
  })
  
  cat(blue('\n------------------------------------------------------------------\n'))
  print(pvalues)
  cat(blue('\n------------------------------------------------------------------\n'))
  pvalues
}
  


library(WriteXLS)



####################################################################
D = getData(patient = '42') # Nothing to do

####################################################################
D = getData(patient = '49')
D$TES2$NR     # Testable
D$fit

Summary = MHT(BBpval(D$TES2$NR, D$fit))

WriteXLS(Summary, 'GBM_testing_results-49.xls')
save(Summary, file = 'RESULTS_TEST-49.RData')

####################################################################
D = getData(patient = '52')
D$TES2$NR
D$TES2$NR['chr2_136610419', 'X52M'] = D$TES1$NR['chr2_136610419', 'SP52M'] + 20
D$TES2$NR

Summary = MHT(BBpval(D$TES2$NR, D$fit))

WriteXLS(Summary, 'GBM_testing_results-52.xls')
save(Summary, file = 'RESULTS_TEST-52.RData')

####################################################################
D = getData(patient = '54')

####################################################################
D = getData(patient = '55')

####################################################################
D = getData(patient = '56')
D$TES2$NR['chr14_75573266', 'X56M3'] = D$TES2$NR['chr14_75573266', 'X56M3'] + D$TES1$NR['chr14_75573266', 'SP56M3']

Summary = MHT(BBpval(D$TES2$NR, D$fit))

WriteXLS(Summary, 'GBM_testing_results-56.xls')
save(Summary, file = 'RESULTS_TEST-56.RData')

####################################################################
D = getData(patient = '57')
D$TES2$NR['chr17_36873741', 'X57M'] = D$TES2$NR['chr17_36873741', 'X57M'] + D$TES1$NR['chr17_36873741', 'SP57M'] + 20
D$TES2$NR['chr19_17769024', 'X57M'] = D$TES2$NR['chr19_17769024', 'X57M'] + D$TES1$NR['chr19_17769024', 'SP57M'] + 20

Summary = MHT(BBpval(D$TES2$NR, D$fit))

WriteXLS(Summary, 'GBM_testing_results-57.xls')
save(Summary, file = 'RESULTS_TEST-57.RData')

####################################################################
D = getData(patient = 'A23')

####################################################################
D = getData(patient = 'A34')

####################################################################
D = getData(patient = 'A44')

WriteXLS(Summary, 'GBM_testing_results-A44.xls')
Summary = MHT(BBpval(D$TES2$NR, D$fit))

####################################################################
D = getData(patient = 'SP28')
D

Summary = MHT(BBpval(
  cbind(D$TES2$NR[, c('A28M', 'R11M')], D$TES1$NR[, c('SP28recM', 'SP28primM', 'CNA')]), 
  D$fit))

WriteXLS(Summary, 'GBM_testing_results-SP28.xls')
save(Summary, file = 'RESULTS_TEST-SP28.RData')

####################################################################
####################################################################
