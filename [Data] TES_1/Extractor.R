library(vcfR)


files = list.files(pattern = '*')
files = files[endsWith(files, '.vcf')]

for(file in files)
{  
  TES1 = read.vcfR(file)
  
  TES1.VCF = NULL
  TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
  TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)
  
  save(TES1.VCF, file = paste(file, '.RData', sep = ''))
}
