library(vcfR)


files = list.files(pattern = '*')
files = files[endsWith(files, '.vcf')]

for(file in files)
{  
  TES2 = read.vcfR(file)
  
  TES2.VCF = NULL
  TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
  TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
  
  rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
  rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')
  
  save(TES2.VCF, file = paste(file, '.RData', sep = ''))
}
