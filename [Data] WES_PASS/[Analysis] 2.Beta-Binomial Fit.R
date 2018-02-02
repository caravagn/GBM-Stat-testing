
patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A44', 'SP28')

for(patient in patients)
{
  # load(paste('CLONAL-', patient, '.RData', sep = ''))
  load(paste('READCOUNTS-', patient, '.RData', sep = ''), verbose = TRUE)
  
  NR = data.frame(WES$NR)
  NV = data.frame(WES$NV)
  
  # samples: last column is CNA status
  samples = colnames(NV)[-ncol(NV)]
  
  # purity
  load(paste('../[Data] CCFs/purity-', patient, '.RData', sep = ''), verbose = TRUE)
  purity = purity[samples , , drop = FALSE]
  
  # correct read counts for purity
  NR = correctReadCounts(NR, purity)
  
  NR = split(NR, f = NR$CNA)
  NV = split(NV, f = NV$CNA)
  
  NR = NR[sapply(NR, function(w) nrow(w) >= MIN_TRAINING_SIZE)]
  NV = NV[sapply(NV, function(w) nrow(w) >= MIN_TRAINING_SIZE)]
  
  fit = BBMLE(NR, NV, samples, patient)
  save(fit, file = paste('BBMLE-', patient, '.RData', sep = ''))
}


