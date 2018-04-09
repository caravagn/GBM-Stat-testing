library(vcfR)

patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A44', 'SP28')

patient = '42'

for(patient in patients)
{  
  file = paste('NG-8132_', patient, '.mutect2.platypus_PASS.vcf', sep = '')
  WES = read.vcfR(file)
  
  WES
  
  WES.VCF = NULL
  WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
  WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)

  rownames(WES.VCF$NV) = paste('chr', rownames(WES.VCF$NV), sep = '')
  rownames(WES.VCF$NR) = paste('chr', rownames(WES.VCF$NR), sep = '')
  
  # Get clonal ones
  load(paste('../[Data] CCFs/CLONAL-', patient, '.RData', sep = ''))
  clonal
  
  NV = cbind(primary(WES.VCF$NV)[rownames(clonal), ], clonal$CNA)
  NR = cbind(primary(WES.VCF$NR)[rownames(clonal), ], clonal$CNA)

  colnames(NV) = colnames(NR) = colnames(clonal) # same orderings
  
  WES = list(NV, NR)
  names(WES) = c('NV', 'NR')
  
  save(WES, file = paste('READCOUNTS-', patient, '.RData', sep = ''))
}
