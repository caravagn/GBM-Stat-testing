library(vcfR)
GIT = '~/Documents/GitHub/GBM-Stat-testing'
CCF.FOLDER = paste(GIT, '/CCFs/', sep = '')
TES1.FOLDER = paste(GIT, '/TES_1/', sep = '')
TES2.FOLDER = paste(GIT, '/TES_2/', sep = '')

######################################################################
###################################################################### Patient 42
######################################################################
setwd(CCF.FOLDER)

CCF = read.csv('Patient_42_cancer_cell_fractions.csv', header = TRUE)

CCF = CCF[, c(3:ncol(CCF))]
colnames(CCF) = c('T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

CCF = CCF[apply(CCF, 1, function(w) !any(is.na(w))), ]
CCF = CCF[apply(CCF, 1, function(w) all(w > 0.8)), ]  

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_42_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(3:ncol(TES1.VCF$NV))]
TES1.VCF$NR = TES1.VCF$NR[, c(3:ncol(TES1.VCF$NR))]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('T1', 'T2', 'T3', 'T4', 'M')

TES1.VCF.filter = TES1.VCF$NR[, 'M'] == 0
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

setwd(TES2.FOLDER)

TES2 = read.vcfR('42.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

TES2.VCF$NV = TES2.VCF$NV[, 1, drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[, 1, drop = FALSE]
colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M')

TES2.VCF.filter = TES2.VCF$NR[, 'M'] == 0 
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

TES2.VCF.filter = rownames(TES2.VCF$NR) %in% rownames(CCF)
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

######################################################################
###################################################################### Patient 49
######################################################################

setwd(CCF.FOLDER)

CCF = read.csv('Patient_49_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

CCF = CCF[, c(3:ncol(CCF))]
colnames(CCF) = c('T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

CCF = CCF[apply(CCF, 1, function(w) !any(is.na(w))), ]
CCF = CCF[apply(CCF, 1, function(w) all(w > 0.8)), ]  

CNA = read.csv('Patient_49_copy_number_per_mutation.csv', header = TRUE)
head(CNA)

CNA = CNA[, c(3,4,5,6)]
colnames(CNA) = c('T1', 'T2', 'T3', 'T4')
rownames(CNA) = paste('chr', rownames(CNA), sep = '')
CNA = CNA[rownames(CCF), ]
CNA = CNA[apply(CNA, 1, function(w) length(unique(w)) == 1),]
CCF = CCF[rownames(CNA), ]

purity = read.csv('Patient_49_purities_cna_analysis.txt', sep = '\t', header = FALSE)
purity = purity[c(1,3,4,5,6), ]
rownames(purity) = c('M', 'T1', 'T2', 'T3', 'T4')
purity$V1 = NULL
colnames(purity) = c('purity')
setwd(TES1.FOLDER)

TES1 = read.vcfR('patient_49_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2:6)]
TES1.VCF$NR = TES1.VCF$NR[, c(2:6)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M', 'T1', 'T2', 'T3', 'T4')

TES1.VCF.filter = TES1.VCF$NV[, 'M'] == 0
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

TES1.VCF.filter = rownames(TES1.VCF$NR) %in% rownames(CCF)
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

setwd(TES2.FOLDER)

TES2 = read.vcfR('49.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

TES2.VCF$NV = TES2.VCF$NV[, 1, drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[, 1, drop = FALSE]
colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M')

TES2.VCF.filter = TES2.VCF$NV[, 'M'] == 0 
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

TES2.VCF.filter = rownames(TES2.VCF$NR) %in% rownames(CCF)
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

hand.del = 'chr7_7106634'
TES2.VCF.filter = !(rownames(TES2.VCF$NR) %in% hand.del)
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

setwd(GIT)

patient = list(
  id = '49',
  CCF = CCF,
  CNA = CNA,
  purity = purity,
  TES1 = TES1.VCF,
  TES2 = TES2.VCF
)
save(patient, file = 'P49.RData')

######################################################################
###################################################################### Patient 52
######################################################################

setwd(CCF.FOLDER)

CCF = read.csv('Patient_52_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

CCF = CCF[, c(3:ncol(CCF))]
colnames(CCF) = c('T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

CCF = CCF[apply(CCF, 1, function(w) !any(is.na(w))), ]
CCF = CCF[apply(CCF, 1, function(w) all(w > 0.8)), ]  

CNA = read.csv('Patient_52_copy_number_per_mutation.csv', header = TRUE)
head(CNA)

CNA = CNA[, c(3,4,5,6)]
colnames(CNA) = c('T1', 'T2', 'T3', 'T4')
rownames(CNA) = paste('chr', rownames(CNA), sep = '')
CNA = CNA[rownames(CCF), ]
CNA = CNA[apply(CNA, 1, function(w) length(unique(w)) == 1),]
CCF = CCF[rownames(CNA), ]

purity = read.csv('Patient_52_purities_cna_analysis.txt', sep = '\t', header = FALSE)
purity = purity[c(1,3,4,5,6), ]
rownames(purity) = c('M', 'T1', 'T2', 'T3', 'T4')
purity$V1 = NULL
colnames(purity) = c('purity')

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_52_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2:6)]
TES1.VCF$NR = TES1.VCF$NR[, c(2:6)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M', 'T1', 'T2', 'T3', 'T4')

TES1.VCF.filter = TES1.VCF$NV[, 'M'] == 0
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

TES1.VCF.filter = rownames(TES1.VCF$NR) %in% rownames(CCF)
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

setwd(TES2.FOLDER)

TES2 = read.vcfR('52.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

TES2.VCF$NV = TES2.VCF$NV[, 1, drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[, 1, drop = FALSE]
colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M')

TES2.VCF.filter = TES2.VCF$NV[, 'M'] == 0 
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

TES2.VCF.filter = rownames(TES2.VCF$NR) %in% rownames(CCF)
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

hand.del = 'chr5_140167486'
TES2.VCF.filter = !(rownames(TES2.VCF$NR) %in% hand.del)
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

setwd(GIT)

patient = list(
  id = '52',
  CCF = CCF,
  CNA = CNA,
  purity = purity,
  TES1 = TES1.VCF,
  TES2 = TES2.VCF
)
save(patient, file = 'P52.RData')

######################################################################
###################################################################### Patient 54
######################################################################

setwd(CCF.FOLDER)

CCF = read.csv('Patient_54_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

CCF = CCF[, c(2:ncol(CCF))]
colnames(CCF) = c('T1', 'T2', 'T3', 'T4', 'T5', 'T6')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

CCF = CCF[apply(CCF, 1, function(w) !any(is.na(w))), ]
CCF = CCF[apply(CCF, 1, function(w) all(w > 0.8)), ]  

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_54_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2:8)]
TES1.VCF$NR = TES1.VCF$NR[, c(2:8)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6')

TES1.VCF.filter = TES1.VCF$NV[, 'M'] == 0
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

TES1.VCF.filter = rownames(TES1.VCF$NR) %in% rownames(CCF)
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

setwd(TES2.FOLDER)

TES2 = read.vcfR('54.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

TES2.VCF$NV = TES2.VCF$NV[, 1, drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[, 1, drop = FALSE]
colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M')

TES2.VCF.filter = TES2.VCF$NV[, 'M'] == 0 
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

TES2.VCF.filter = rownames(TES2.VCF$NR) %in% rownames(CCF)
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

######################################################################
###################################################################### Patient 55
######################################################################

setwd(CCF.FOLDER)

CCF = read.csv('Patient_55_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

CCF = CCF[, c(2:ncol(CCF))]
colnames(CCF) = c('T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

CCF = CCF[apply(CCF, 1, function(w) !any(is.na(w))), ]
CCF = CCF[apply(CCF, 1, function(w) all(w > 0.8)), ]  

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_55_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2:6)]
TES1.VCF$NR = TES1.VCF$NR[, c(2:6)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M', 'T1', 'T2', 'T3', 'T4')

TES1.VCF.filter = TES1.VCF$NV[, 'M'] == 0
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

TES1.VCF.filter = rownames(TES1.VCF$NR) %in% rownames(CCF)
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

setwd(TES2.FOLDER)

TES2 = read.vcfR('55.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')


######################################################################
###################################################################### Patient 56
######################################################################

setwd(CCF.FOLDER)

CCF = read.csv('Patient_56_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

CCF = CCF[, c(3:ncol(CCF))]
colnames(CCF) = c('T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

CCF = CCF[apply(CCF, 1, function(w) !any(is.na(w))), ]
CCF = CCF[apply(CCF, 1, function(w) all(w > 0.8)), ]  
CCF = CCF[apply(CCF, 1, function(w) !any(is.na(w))), ]

CNA = read.csv('Patient_56_copy_number_per_mutation.csv', header = TRUE)
head(CNA)

CNA = CNA[, c(3,4,5,6)]
colnames(CNA) = c('T1', 'T2', 'T3', 'T4')
rownames(CNA) = paste('chr', rownames(CNA), sep = '')
CNA = CNA[rownames(CCF), ]
CNA = CNA[apply(CNA, 1, function(w) length(unique(w)) == 1),]
CCF = CCF[rownames(CNA), ]

purity = read.csv('Patient_56_purities_cna_analysis.txt', sep = '\t', header = FALSE)
purity = purity[c(1,3,4,5,6), ]
rownames(purity) = c('M', 'T1', 'T2', 'T3', 'T4')
purity$V1 = NULL
colnames(purity) = c('purity')

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_56_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2,4,5,6,7)]
TES1.VCF$NR = TES1.VCF$NR[, c(2,4,5,6,7)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M', 'T1', 'T2', 'T3', 'T4')

TES1.VCF.filter = TES1.VCF$NV[, 'M'] == 0
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

TES1.VCF.filter = rownames(TES1.VCF$NR) %in% rownames(CCF)
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

setwd(TES2.FOLDER)

TES2 = read.vcfR('56.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

TES2.VCF$NV = TES2.VCF$NV[, 1, drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[, 1, drop = FALSE]
colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M')

TES2.VCF.filter = TES2.VCF$NV[, 'M'] == 0 
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

TES2.VCF.filter = rownames(TES2.VCF$NR) %in% rownames(CCF)
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

setwd(GIT)

patient = list(
  id = '56',
  CCF = CCF,
  CNA = CNA,
  purity = purity,
  TES1 = TES1.VCF,
  TES2 = TES2.VCF
)
save(patient, file = 'P56.RData')

######################################################################
###################################################################### Patient 57
######################################################################

setwd(CCF.FOLDER)

CCF = read.csv('Patient_57_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

CCF = CCF[, c(3:ncol(CCF))]
colnames(CCF) = c('T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

CCF = CCF[apply(CCF, 1, function(w) !any(is.na(w))), ]
CCF = CCF[apply(CCF, 1, function(w) all(w > 0.8)), ]  

CNA = read.csv('Patient_57_copy_number_per_mutation.csv', header = TRUE)
head(CNA)

CNA = CNA[, c(3:ncol(CNA))]
colnames(CNA) = c('T1', 'T2', 'T3', 'T4')
rownames(CNA) = paste('chr', rownames(CNA), sep = '')
CNA = CNA[rownames(CCF), ]
CNA = CNA[apply(CNA, 1, function(w) length(unique(w)) == 1),]
CCF = CCF[rownames(CNA), ]

purity = read.csv('Patient_57_purities_cna_analysis.txt', sep = '\t', header = FALSE)
purity = purity[c(1,3,4,5,6), ]
rownames(purity) = c('M', 'T1', 'T2', 'T3', 'T4')
purity$V1 = NULL
colnames(purity) = c('purity')

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_57_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2,4,5,6,7)]
TES1.VCF$NR = TES1.VCF$NR[, c(2,4,5,6,7)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M', 'T1', 'T2', 'T3', 'T4')

TES1.VCF.filter = TES1.VCF$NV[, 'M'] == 0
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

TES1.VCF.filter = rownames(TES1.VCF$NR) %in% rownames(CCF)
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

setwd(TES2.FOLDER)

TES2 = read.vcfR('57.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

TES2.VCF$NV = TES2.VCF$NV[, 1, drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[, 1, drop = FALSE]
colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M')

TES2.VCF.filter = TES2.VCF$NV[, 'M'] == 0 
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

TES2.VCF.filter = rownames(TES2.VCF$NR) %in% rownames(CCF)
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

setwd(GIT)

patient = list(
  id = '57',
  CCF = CCF,
  CNA = CNA,
  purity = purity,
  TES1 = TES1.VCF,
  TES2 = TES2.VCF
)
save(patient, file = 'P57.RData')

######################################################################
###################################################################### Patient A23
######################################################################

setwd(CCF.FOLDER)

CCF = read.csv('Patient_A23_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

CCF = CCF[, c(3,6)]
colnames(CCF) = c('T1', 'T2')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

CCF = CCF[apply(CCF, 1, function(w) !any(is.na(w))), ]
CCF = CCF[apply(CCF, 1, function(w) all(w > 0.8)), ]  

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_SP19_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2,3,4,6)]
TES1.VCF$NR = TES1.VCF$NR[, c(2,3,4,6)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M1', 'T1', 'M2', 'T2')

TES1.VCF.filter = TES1.VCF$NV[, 'M1'] == 0 &  TES1.VCF$NV[, 'M2'] == 0
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

TES1.VCF.filter = rownames(TES1.VCF$NR) %in% rownames(CCF)
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

setwd(TES2.FOLDER)

TES2 = read.vcfR('A23_SP19.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

TES2.VCF$NV = TES2.VCF$NV[, c(1, 3), drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[, c(1, 3), drop = FALSE]
colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M1', 'M2')

TES2.VCF.filter = TES2.VCF$NV[, 'M1'] == 0 & TES2.VCF$NV[, 'M2'] == 0
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

TES2.VCF.filter = rownames(TES2.VCF$NR) %in% rownames(CCF)
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

######################################################################
###################################################################### Patient A34
######################################################################

setwd(CCF.FOLDER)

CCF = read.csv('Patient_A34_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

CCF = CCF[, c(4:ncol(CCF))]
colnames(CCF) = c('T1', 'T2', 'T3', 'T4', 'T5', 'T6')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

CCF = CCF[apply(CCF, 1, function(w) !any(is.na(w))), ]
CCF = CCF[apply(CCF, 1, function(w) all(w > 0.8)), ]  

######################################################################
###################################################################### Patient A44
######################################################################

setwd(CCF.FOLDER)

CCF = read.csv('Patient_A44_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

CCF = CCF[, c(3:ncol(CCF))]
colnames(CCF) = c('T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

CCF = CCF[apply(CCF, 1, function(w) !any(is.na(w))), ]
CCF = CCF[apply(CCF, 1, function(w) all(w > 0.8)), ]  

CNA = read.csv('Patient_A44_copy_number_per_mutation.csv', header = TRUE)
head(CNA)
 
CNA = CNA[, c(3,4,5,6)]
colnames(CNA) = c('T1', 'T2', 'T3', 'T4')
rownames(CNA) = paste('chr', rownames(CNA), sep = '')
CNA = CNA[rownames(CCF), ]
CNA = CNA[apply(CNA, 1, function(w) length(unique(w)) == 1),]
CCF = CCF[rownames(CNA), ]

purity = read.csv('Patient_A44_purities_cna_analysis.txt', sep = '\t', header = FALSE)
purity = purity[c(1,3,4,5,6), ]
rownames(purity) = c('M', 'T1', 'T2', 'T3', 'T4')
purity$V1 = NULL
colnames(purity) = c('purity')

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_A44_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2,4,5,6,7)]
TES1.VCF$NR = TES1.VCF$NR[, c(2,4,5,6,7)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M', 'T1', 'T2', 'T3', 'T4')

TES1.VCF.filter = TES1.VCF$NV[, 'M'] == 0
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

TES1.VCF.filter = rownames(TES1.VCF$NR) %in% rownames(CCF)
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

setwd(TES2.FOLDER)

TES2 = read.vcfR('A44.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

TES2.VCF$NV = TES2.VCF$NV[, 1, drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[, 1, drop = FALSE]
colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M')

TES2.VCF.filter = TES2.VCF$NV[, 'M'] == 0 
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

TES2.VCF.filter = rownames(TES2.VCF$NR) %in% rownames(CCF)
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

setwd(GIT)

patient = list(
  id = 'A44',
  CCF = CCF,
  CNA = CNA,
  purity = purity,
  TES1 = TES1.VCF,
  TES2 = TES2.VCF
)
save(patient, file = 'PA44.RData')


######################################################################
###################################################################### Patient SP28
######################################################################

setwd(CCF.FOLDER)

CCF = read.csv('Patient_SP28_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

CCF = CCF[, c(3,6)]
colnames(CCF) = c('T1', 'T2')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

CCF = CCF[apply(CCF, 1, function(w) !any(is.na(w))), ]
CCF = CCF[apply(CCF, 1, function(w) all(w > 0.8)), ]  

CNA = read.csv('Patient_SP28_copy_number_per_mutation.csv', header = TRUE)
head(CNA)

CNA = CNA[, c(3,6)]
colnames(CNA) = c('T1','T2')
rownames(CNA) = paste('chr', rownames(CNA), sep = '')
CNA = CNA[rownames(CCF), ]
CNA = CNA[apply(CNA, 1, function(w) length(unique(w)) == 1),]
CCF = CCF[rownames(CNA), ]

purity = read.csv('Patient_SP28_purities_cna_analysis.txt', sep = '\t', header = FALSE)
purity = purity[c(1,3,4,6), ]
rownames(purity) = c('M1', 'T1', 'M2', 'T2')
purity$V1 = NULL
colnames(purity) = c('purity')

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_SP28_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(1,3,5,7)]
TES1.VCF$NR = TES1.VCF$NR[, c(1,3,5,7)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M1', 'T1', 'M2', 'T2')

TES1.VCF.filter = TES1.VCF$NV[, 'M1'] == 0 & TES1.VCF$NV[, 'M2'] == 0
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

TES1.VCF.filter = rownames(TES1.VCF$NR) %in% rownames(CCF)
TES1.VCF$NV = TES1.VCF$NV[TES1.VCF.filter, , drop = FALSE]
TES1.VCF$NR = TES1.VCF$NR[TES1.VCF.filter, , drop = FALSE]

setwd(TES2.FOLDER)

TES2 = read.vcfR('SP28_R11.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

TES2.VCF$NV = TES2.VCF$NV[, c(2,3), drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[, c(2,3), drop = FALSE]
colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M1', 'M2')

TES2.VCF.filter = TES2.VCF$NV[, 'M1'] == 0 & TES2.VCF$NV[, 'M2'] == 0 
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

TES2.VCF.filter = rownames(TES2.VCF$NR) %in% rownames(CCF)
TES2.VCF$NV = TES2.VCF$NV[TES2.VCF.filter, , drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[TES2.VCF.filter, , drop = FALSE]

setwd(GIT)

patient = list(
  id = 'SP28',
  CCF = CCF,
  CNA = CNA,
  purity = purity,
  TES1 = TES1.VCF,
  TES2 = TES2.VCF
)
save(patient, file = 'PSP28.RData')
