library(vcfR)

margin.s = function(x) 
{
  x[, !grepl('T', colnames(x)), drop = FALSE]
}

## NOTE: file A23_SP19.platypus was renamed to A23.platypus because that's thes same patient
## NOTE: file SP28_R11.platypus was renamed to SP28.platypus because that's thes same patient
patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A34', 'A44', 'SP28')
library(vcfR)

for(patient in patients)
{  
  file = paste(patient, '.platypus.vcf', sep = '')
  TES2 = read.vcfR(file)
  
  TES2.VCF = NULL
  TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
  TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
  
  rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '') # this panel does not have this..
  rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '') # this panel does not have this..
  
  # TES2 contains SNVs from other patients, we do not want to use them
  load(paste('../[Data] CCFs/CCF-', patient, '.RData', sep = ''), verbose = TRUE)
  TES2.VCF$NV = TES2.VCF$NV[rownames(TES2.VCF$NV) %in% rownames(CCF), , drop = FALSE]
  TES2.VCF$NR = TES2.VCF$NR[rownames(TES2.VCF$NR) %in% rownames(CCF), , drop = FALSE]
  
  # Get margin/S ones
  # NV = margin.s(TES2.VCF$NV)
  # NR = margin.s(TES2.VCF$NR)
  NV = TES2.VCF$NV
  NR = TES2.VCF$NR
  
    
  colnames(NV) = colnames(NR) # same orderings
  
  TES2 = list(NV, NR)
  names(TES2) = c('NV', 'NR')
  
  save(TES2, file = paste('MARGIN_S_READCOUNTS-', patient, '.RData', sep = ''))
}
