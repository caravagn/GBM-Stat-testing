library(vcfR)


files = list.files(pattern = '*')
files = files[endsWith(files, '.vcf')]

for(file in files)
{  
  WES = read.vcfR(file)
  
  WES.VCF = NULL
  WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
  WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)

  rownames(WES.VCF$NV) = paste('chr', rownames(WES.VCF$NV), sep = '')
  rownames(WES.VCF$NR) = paste('chr', rownames(WES.VCF$NR), sep = '')
  
  save(WES.VCF, file = paste(file, '.RData', sep = ''))
}
