library(vcfR)

margin.s = function(x) 
{
  x[, !grepl('T', colnames(x)), drop = FALSE]
}

patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A44', 'SP28')

for(patient in patients)
{  
  file = paste('NG-8132_', patient, '.mutect2.platypus_PASS.vcf', sep = '')
  WES = read.vcfR(file)
  
  WES.VCF = NULL
  WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
  WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)
  
  rownames(WES.VCF$NV) = paste('chr', rownames(WES.VCF$NV), sep = '')
  rownames(WES.VCF$NR) = paste('chr', rownames(WES.VCF$NR), sep = '')
  
  # Get clonal ones
  NV = margin.s(WES.VCF$NV)
  NR = margin.s(WES.VCF$NR)
  
  colnames(NV) = colnames(NR) # same orderings
  
  WES = list(NV, NR)
  names(WES) = c('NV', 'NR')
  
  save(WES, file = paste('MARGIN_S_READCOUNTS-', patient, '.RData', sep = ''))
}
