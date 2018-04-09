
ismargin = function(x) 
{
  x[, grepl('M', colnames(x)), drop = FALSE]
}

## NOTE: file A23_SP19.platypus was renamed to A23.platypus because that's thes same patient
## NOTE: file SP28_R11.platypus was renamed to SP28.platypus because that's thes same patient
patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A34', 'A44', 'SP28')
library(vcfR)

MIN.READS.CUTOFF.TESTABLE = 10

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
  TES2.VCF$NR = TES2.VCF$NR[rownames(TES2.VCF$NR) %in% rownames(CCF), , drop = FALSE]
  TES2.VCF$NV = TES2.VCF$NV[rownames(TES2.VCF$NV) %in% rownames(CCF), , drop = FALSE]
  
  # in the margin -- WithMutantAllele/ 0 reads in the margin
  NV = ismargin(TES2.VCF$NV)
  WithMutantAllele = rownames(NV[rowSums(NV) >= MIN.READS.CUTOFF.TESTABLE, , drop = F])
  NV = NV[rowSums(NV) < MIN.READS.CUTOFF.TESTABLE, , drop = F]
  
  TES2.NR = ismargin(TES2.VCF$NR)
  TES2.NR = TES2.NR[rownames(NV), , drop = F]
  
  # Purity -- we correct for purity the coverage (worst case scenario)
  purity = data.frame(
    pi = rep(WORST_CASE_SCENARIO_PURITY, ncol(TES2.NR)),
    row.names = colnames(TES2.NR))
  
  # Clonal ones are testable
  load(paste('../[Data] CCFs/CLONAL-', patient, '.RData', sep = ''), verbose = TRUE)
  TES2.NR = TES2.NR[rownames(TES2.NR) %in% rownames(clonal), , drop = FALSE]
  
  # Add CNA
  TES2.NR = cbind(TES2.NR, CNA = clonal[rownames(TES2.NR), 'CNA'])
  
  TES2 = NULL
  TES2$NR = correctReadCounts(TES2.NR, purity)
  TES2$Rejected = WithMutantAllele
  
  save(TES2, file = paste('NR-', patient, '.RData', sep = ''))
}

