library(vcfR)

primary = function(x) 
{
  x[, grepl('T', colnames(x)), drop = FALSE]
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
  
  
  VAF = matrixcalc::hadamard.prod(WES.VCF$NV, 1/WES.VCF$NR)
  VAF = primary(VAF)
  head(VAF)
                            
  
  save(VAF, file = paste('VAF-', patient, '.RData', sep = ''))
}
