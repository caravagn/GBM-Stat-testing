library(vcfR)

ismargin = function(x) 
{
  x[, grepl('M', colnames(x)), drop = FALSE]
}

## NOTE: file patient_SP19_platypus was renamed to patient_A23_platypus because that's thes same patient
patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A34', 'A44', 'SP28')


for(patient in patients)
{  
  file = paste('patient_', patient, '_platypus.vcf', sep = '')
  TES1 = read.vcfR(file)
  
  TES1.VCF = NULL
  TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
  TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

  # in the margin -- WithMutantAllele/ 0 reads in the margin
  NV = ismargin(TES1.VCF$NV)
  WithMutantAllele = rownames(NV[rowSums(NV) >= MIN.READS.CUTOFF.TESTABLE, , drop = F])
  NV = NV[rowSums(NV) < MIN.READS.CUTOFF.TESTABLE, , drop = F]
  
  TES1.NR = ismargin(TES1.VCF$NR)
  TES1.NR = TES1.NR[rownames(NV), , drop = F]
  
  # Clonal ones are testable
  load(paste('../[Data] CCFs/CLONAL-', patient, '.RData', sep = ''))
  TES1.NR = TES1.NR[rownames(TES1.NR) %in% rownames(clonal), , drop = FALSE]
  
  # Purity -- we correct for purity the coverage (worst case scenario)
  purity = data.frame(
    pi = rep(WORST_CASE_SCENARIO_PURITY, ncol(TES1.NR)),
    row.names = colnames(TES1.NR))
  
  # Add CNA
  TES1.NR = cbind(TES1.NR, CNA = clonal[rownames(TES1.NR), 'CNA'])
  
  TES1 = NULL
  TES1$NR = correctReadCounts(TES1.NR, purity)
  TES1$Rejected = WithMutantAllele
  
  save(TES1, file = paste('NR-', patient, '.RData', sep = ''))
}


