sv = function(CCF, patient) {
  
  rownames(CCF) = paste('chr', rownames(CCF), sep = '')
  CCF[apply(CCF, 2, is.infinite)] = 0
  CCF[CCF > 1] = 1
  
  save(CCF, file = paste('CCF-', patient, '.RData', sep = ''))
}

CCF = read.csv("Patient_42_cancer_cell_fractions.csv", header = TRUE)
colnames(CCF) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
sv(CCF, '42')

CCF = read.csv("Patient_49_cancer_cell_fractions.csv", header = TRUE)
colnames(CCF) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
sv(CCF, '49')

CCF = read.csv("Patient_52_cancer_cell_fractions.csv", header = TRUE)
colnames(CCF) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
sv(CCF, '52')

CCF = read.csv("Patient_54_cancer_cell_fractions.csv", header = TRUE)
colnames(CCF) = c('M1', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6')
sv(CCF, '54')

CCF = read.csv("Patient_55_cancer_cell_fractions.csv", header = TRUE)
colnames(CCF) = c('S', 'T1', 'T2', 'T3', 'T4')
sv(CCF, '55')

CCF = read.csv("Patient_56_cancer_cell_fractions.csv", header = TRUE)
colnames(CCF) = c('M3', 'S', 'T1', 'T2', 'T3', 'T4')
sv(CCF, '56')

CCF = read.csv("Patient_57_cancer_cell_fractions.csv", header = TRUE)
colnames(CCF) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
sv(CCF, '57')

CCF = read.csv("Patient_A23_cancer_cell_fractions.csv", header = TRUE)
colnames(CCF) = c('M1', 'S1', 'T1', 'M2', 'S2', 'T2')
sv(CCF, 'A23')

CCF = read.csv("Patient_A34_cancer_cell_fractions.csv", header = TRUE)
colnames(CCF) = c('S1', 'S2', 'S3', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6')
sv(CCF, 'A34')

CCF = read.csv("Patient_A44_cancer_cell_fractions.csv", header = TRUE)
colnames(CCF) = c('M', 'S', 'T1', 'T2', 'T3', 'T5')
sv(CCF, 'A44')

CCF = read.csv("Patient_SP28_cancer_cell_fractions.csv", header = TRUE)
colnames(CCF) = c('M1', 'S1', 'T1', 'M2', 'S2', 'T2')
sv(CCF, 'SP28')


############ CNA

sv = function(CNA, patient) {
  rownames(CNA) = paste('chr', rownames(CNA), sep = '')
  CNA[apply(CNA, 2, is.infinite)] = 0
  save(CNA, file = paste('CNA-', patient, '.RData', sep = ''))
}

CNA = read.csv("Patient_42_copy_number_per_mutation.csv", header = TRUE)
colnames(CNA) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
sv(CNA, '42')

CNA = read.csv("Patient_49_copy_number_per_mutation.csv", header = TRUE)
colnames(CNA) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
sv(CNA, '49')

CNA = read.csv("Patient_52_copy_number_per_mutation.csv", header = TRUE)
colnames(CNA) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
sv(CNA, '52')

CNA = read.csv("Patient_54_copy_number_per_mutation.csv", header = TRUE)
colnames(CNA) = c('M1', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6')
sv(CNA, '54')

CNA = read.csv("Patient_55_copy_number_per_mutation.csv", header = TRUE)
colnames(CNA) = c('S', 'T1', 'T2', 'T3', 'T4')
sv(CNA, '55')

CNA = read.csv("Patient_56_copy_number_per_mutation.csv", header = TRUE)
colnames(CNA) = c('M3', 'S', 'T1', 'T2', 'T3', 'T4')
sv(CNA, '56')

CNA = read.csv("Patient_57_copy_number_per_mutation.csv", header = TRUE)
colnames(CNA) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
sv(CNA, '57')

CNA = read.csv("Patient_A23_copy_number_per_mutation.csv", header = TRUE)
colnames(CNA) = c('M1', 'S1', 'T1', 'M2', 'S2', 'T2')
sv(CNA, 'A23')

CNA = read.csv("Patient_A34_copy_number_per_mutation.csv", header = TRUE)
colnames(CNA) = c('S1', 'S2', 'S3', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6')
sv(CNA, 'A34')

CCF = read.csv("Patient_A44_copy_number_per_mutation.csv", header = TRUE)
colnames(CCF) = c('M', 'S', 'T1', 'T2', 'T3', 'T5')
sv(CCF, 'A44')

CNA = read.csv("Patient_SP28_copy_number_per_mutation.csv", header = TRUE)
colnames(CNA) = c('M1', 'S1', 'T1', 'M2', 'S2', 'T2')
sv(CNA, 'SP28')


############ Purity
sv = function(purity, patient) {
  rownames(purity) = purity[, 1]
  purity$V1 = NULL
  save(purity, file = paste('purity-', patient, '.RData', sep = ''))
}

purity = read.csv("Patient_42_purities_cna_analysis.txt", header = FALSE, sep = '\t')
purity[, 1] = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
sv(purity, '42')

purity = read.csv("Patient_49_purities_cna_analysis.txt", header = FALSE, sep = '\t')
purity[, 1] = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
sv(purity, '49')

purity = read.csv("Patient_52_purities_cna_analysis.txt", header = FALSE, sep = '\t')
purity[, 1] = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
sv(purity, '52')

purity = read.csv("Patient_54_purities_cna_analysis.txt", header = FALSE, sep = '\t')
purity[, 1] = c('M1', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6')
sv(purity, '54')

purity = read.csv("Patient_55_purities_cna_analysis.txt", header = FALSE, sep = '\t')
purity[, 1] = c('S', 'T1', 'T2', 'T3', 'T4')
sv(purity, '55')

purity = read.csv("Patient_56_purities_cna_analysis.txt", header = FALSE, sep = '\t')
purity[, 1] = c('M3', 'S', 'T1', 'T2', 'T3', 'T4')
sv(purity, '56')

purity = read.csv("Patient_57_purities_cna_analysis.txt", header = FALSE, sep = '\t')
purity[, 1] = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
sv(purity, '57')

purity = read.csv("Patient_A23_purities_cna_analysis.txt", header = FALSE, sep = '\t')
purity[, 1] = c('M1', 'S1', 'T1', 'M2', 'S2', 'T2')
sv(purity, 'A23')

purity = read.csv("Patient_A34_purities_cna_analysis.txt", header = FALSE, sep = '\t')
purity[, 1] = c('S1', 'S2', 'S3', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6')
sv(purity, 'A34')

purity = read.csv("Patient_A44_purities_cna_analysis.txt", header = FALSE, sep = '\t')
purity[, 1] = c('M', 'S', 'T1', 'T2', 'T3', 'T5')
sv(purity, 'A44')

purity = read.csv("Patient_SP28_purities_cna_analysis.txt", header = FALSE, sep = '\t')
purity[, 1] = c('M1', 'S1', 'T1', 'M2', 'S2', 'T2')
sv(purity, 'SP28')




