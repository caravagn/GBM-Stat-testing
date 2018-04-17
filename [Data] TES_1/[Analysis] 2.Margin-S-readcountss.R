library(vcfR)

margin.s = function(x) 
{
  x[, !grepl('T', colnames(x)), drop = FALSE]
}

for(patient in patients)
{  
  file = paste('patient_', patient, '_platypus.vcf', sep = '')
  TES1 = read.vcfR(file)
  
  TES1.VCF = NULL
  TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
  TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)
  
  # Get clonal ones
  # NV = margin.s(TES1.VCF$NV)
  # NR = margin.s(TES1.VCF$NR)
  NV = TES1.VCF$NV
  NR = TES1.VCF$NR
  
  colnames(NV) = colnames(NR) # same orderings
  
  TES1 = list(NV, NR)
  names(TES1) = c('NV', 'NR')
  
  save(TES1, file = paste('MARGIN_S_READCOUNTS-', patient, '.RData', sep = ''))
}
