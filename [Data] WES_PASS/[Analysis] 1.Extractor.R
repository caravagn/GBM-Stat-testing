library(vcfR)

primary = function(x) 
{
  x[, grepl('T', colnames(x)), drop = FALSE]
}


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
  load(paste('../[Data] CCFs/CLONAL-', patient, '.RData', sep = ''))
  clonal
  
  dfNV = as.data.frame(primary(WES.VCF$NV)[rownames(clonal), ])
  dfNR = as.data.frame(primary(WES.VCF$NR)[rownames(clonal), ])
  
  if(nrow(clonal) == 1) {
    dfNV = t(dfNV)
    dfNR = t(dfNR)
  }
  
  NV = cbind(dfNV, data.frame(CNA = clonal$CNA))
  NR = cbind(dfNR, data.frame(CNA = clonal$CNA))
  
  colnames(NV) = colnames(NR) = colnames(clonal) # same orderings
  rownames(NV) = rownames(NR) = rownames(clonal)
  
  WES = list(NV, NR)
  names(WES) = c('NV', 'NR')
  
  save(WES, file = paste('READCOUNTS-', patient, '.RData', sep = ''))
}
