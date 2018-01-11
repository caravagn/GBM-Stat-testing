########## Author: Giulio Caravagna, ICR. <giulio.caravagna@icr.ac.uk> or <gcaravagn@gmail.com>
##########
########## This script will analyze data for each patient. It will:
########## - load data prepared by the previous script, plus some extra info, and will reshape it in some data frame
########## - create training and test

MINTR.GROUP.SIZE = 10
PVALUE = 0.05
WORST.PURITY = 0.01
library(vcfR)
source('[Script] Auxiliary.lib.R')

######################################################################
###################################################################### Patient 49
######################################################################

load('P49.RData')

training.VCF = read.vcfR('[Data] WES_PASS/NG-8132_49.mutect2.platypus_PASS.vcf')

training = NULL
training$NV <- extract.gt(training.VCF, element='NV', as.numeric=TRUE)
training$NR <- extract.gt(training.VCF, element='NR', as.numeric=TRUE)
rownames(training$NV) = rownames(training$NR) = paste('chr', rownames(training$NV), sep = '')

training$NV = training$NV[, 4:7]
training$NR = training$NR[, 4:7]
colnames(training$NV) = colnames(training$NR) = c('T1', 'T2', 'T3', 'T4') 
training$NV = training$NV[rownames(patient$CCF), ]
training$NR = training$NR[rownames(patient$CCF), ]

rownames(patient$CNA) == rownames(training$NR)

colnames(training$NR) = paste('TRNR-', colnames(training$NR), sep = '')
colnames(training$NV)= paste('TRNV-', colnames(training$NV), sep = '')
training$NV = cbind(training$NV, CN = patient$CNA[, 1])

data = cbind(training$NR, training$NV)
data = cbind(data, TEST = NA)
data[rownames(patient$TES1$NR), 'TEST'] = patient$TES1$NR[, 'M']
data[rownames(patient$TES2$NR), 'TEST'] = patient$TES2$NR[, 'M']

data = data.frame(data)
data = correctReadCounts(data, 'TRNR.T1', patient$purity['T1', ])
data = correctReadCounts(data, 'TRNR.T2', patient$purity['T2', ])
data = correctReadCounts(data, 'TRNR.T3', patient$purity['T3', ])
data = correctReadCounts(data, 'TRNR.T4', patient$purity['T4', ])
#data = correctReadCounts(data, 'TEST', 0.055)
data = correctReadCounts(data, 'TEST', WORST.PURITY)
write.csv(data, 'a.txt')

data = split(data, f = data$CN)
filter = sapply(data, function(w) all(is.na(w$TEST)) | nrow(w) < MINTR.GROUP.SIZE)
data = data[!filter]


pdf('Patient49-BBfit.pdf')
training.params = lapply(c('T1', 'T2', 'T3', 'T4'),
       BBMLE,
       main = 'Patient 49',
       data = data)
dev.off()
jamPDF('Patient49-BBfit.pdf', out.file = 'Patient49-BBfit.pdf', layout = '2x2')
training.params = Reduce(rbind, training.params)

toTest = lapply(data, function(w) w[!is.na(w$TEST),  'TEST', drop = FALSE])
TestsTable = BBpval(training.params, toTest, PVALUE)

patient$TestsTable = TestsTable
patient$data = data
patient$training.params = training.params

save(patient, file = 'P49.RData')

######################################################################
###################################################################### Patient 52
######################################################################

load('P52.RData')

training.VCF = read.vcfR('[Data] WES_PASS/NG-8132_52.mutect2.platypus_PASS.vcf')

training = NULL
training$NV <- extract.gt(training.VCF, element='NV', as.numeric=TRUE)
training$NR <- extract.gt(training.VCF, element='NR', as.numeric=TRUE)
rownames(training$NV) = rownames(training$NR) = paste('chr', rownames(training$NV), sep = '')
head(training$NV)

training$NV = training$NV[, 4:7]
training$NR = training$NR[, 4:7]
colnames(training$NV) = colnames(training$NR) = c('T1', 'T2', 'T3', 'T4') 
training$NV = training$NV[rownames(patient$CCF), ]
training$NR = training$NR[rownames(patient$CCF), ]

rownames(patient$CNA) == rownames(training$NR)

colnames(training$NR) = paste('TRNR-', colnames(training$NR), sep = '')
colnames(training$NV)= paste('TRNV-', colnames(training$NV), sep = '')
training$NV = cbind(training$NV, CN = patient$CNA[, 1])

data = cbind(training$NR, training$NV)
data = cbind(data, TEST = NA)

which.test = c(patient$TES1$NR[, 'M'], patient$TES2$NR[, 'M'])
which.test = tapply(unlist(which.test), names(unlist(which.test)), sum)

data[names(which.test), 'TEST'] = which.test

data = data.frame(data)
data = correctReadCounts(data, 'TRNR.T1', patient$purity['T1', ])
data = correctReadCounts(data, 'TRNR.T2', patient$purity['T2', ])
data = correctReadCounts(data, 'TRNR.T3', patient$purity['T3', ])
data = correctReadCounts(data, 'TRNR.T4', patient$purity['T4', ])
#data = correctReadCounts(data, 'TEST', 0.055)
data = correctReadCounts(data, 'TEST', WORST.PURITY)

data = split(data, f = data$CN)
filter = sapply(data, function(w) all(is.na(w$TEST)) | nrow(w) < MINTR.GROUP.SIZE)
data = data[!filter]

pdf('Patient52-BBfit.pdf')
training.params = lapply(c('T1', 'T2', 'T3', 'T4'),
                         BBMLE,
                         main = 'Patient 52',
                         data = data)
dev.off()
jamPDF('Patient52-BBfit.pdf', out.file = 'Patient52-BBfit.pdf', layout = '2x2')
training.params = Reduce(rbind, training.params)

toTest = lapply(data, function(w) w[!is.na(w$TEST),  'TEST', drop = FALSE])
TestsTable = BBpval(training.params, toTest, PVALUE)

patient$TestsTable = TestsTable
patient$data = data
patient$training.params = training.params

save(patient, file = 'P52.RData')

######################################################################
###################################################################### Patient 56
######################################################################

load('P56.RData')

training.VCF = read.vcfR('[Data] WES_PASS/NG-8132_56.mutect2.platypus_PASS.vcf')

training = NULL
training$NV <- extract.gt(training.VCF, element='NV', as.numeric=TRUE)
training$NR <- extract.gt(training.VCF, element='NR', as.numeric=TRUE)
rownames(training$NV) = rownames(training$NR) = paste('chr', rownames(training$NV), sep = '')
head(training$NV)

training$NV = training$NV[, 4:7]
training$NR = training$NR[, 4:7]
colnames(training$NV) = colnames(training$NR) = c('T1', 'T2', 'T3', 'T4') 
training$NV = training$NV[rownames(patient$CCF), ]
training$NR = training$NR[rownames(patient$CCF), ]

rownames(patient$CNA) == rownames(training$NR)

colnames(training$NR) = paste('TRNR-', colnames(training$NR), sep = '')
colnames(training$NV)= paste('TRNV-', colnames(training$NV), sep = '')
training$NV = cbind(training$NV, CN = patient$CNA[, 1])

data = cbind(training$NR, training$NV)
data = cbind(data, TEST = NA)

which.test = c(patient$TES1$NR[, 'M'], patient$TES2$NR[, 'M'])
which.test = tapply(unlist(which.test), names(unlist(which.test)), sum)

data[names(which.test), 'TEST'] = which.test

data = data.frame(data)
data = correctReadCounts(data, 'TRNR.T1', patient$purity['T1', ])
data = correctReadCounts(data, 'TRNR.T2', patient$purity['T2', ])
data = correctReadCounts(data, 'TRNR.T3', patient$purity['T3', ])
data = correctReadCounts(data, 'TRNR.T4', patient$purity['T4', ])
#data = correctReadCounts(data, 'TEST', 0.055)
data = correctReadCounts(data, 'TEST', WORST.PURITY)

data = split(data, f = data$CN)
filter = sapply(data, function(w) all(is.na(w$TEST)) | nrow(w) < MINTR.GROUP.SIZE)
data = data[!filter]

pdf('Patient56-BBfit.pdf')
training.params = lapply(c('T1', 'T2', 'T3', 'T4'),
                         BBMLE,
                         main = 'Patient 56',
                         data = data)
dev.off()
jamPDF('Patient56-BBfit.pdf', out.file = 'Patient56-BBfit.pdf', layout = '2x2')
training.params = Reduce(rbind, training.params)

toTest = lapply(data, function(w) w[!is.na(w$TEST),  'TEST', drop = FALSE])
TestsTable = BBpval(training.params, toTest, PVALUE)

patient$TestsTable = TestsTable
patient$data = data
patient$training.params = training.params

save(patient, file = 'P56.RData')

######################################################################
###################################################################### Patient 57
######################################################################

load('P57.RData')

training.VCF = read.vcfR('[Data] WES_PASS/NG-8132_57.mutect2.platypus_PASS.vcf')

training = NULL
training$NV <- extract.gt(training.VCF, element='NV', as.numeric=TRUE)
training$NR <- extract.gt(training.VCF, element='NR', as.numeric=TRUE)
rownames(training$NV) = rownames(training$NR) = paste('chr', rownames(training$NV), sep = '')
head(training$NV)

training$NV = training$NV[, 4:7]
training$NR = training$NR[, 4:7]
colnames(training$NV) = colnames(training$NR) = c('T1', 'T2', 'T3', 'T4') 
training$NV = training$NV[rownames(patient$CCF), ]
training$NR = training$NR[rownames(patient$CCF), ]

rownames(patient$CNA) == rownames(training$NR)

colnames(training$NR) = paste('TRNR-', colnames(training$NR), sep = '')
colnames(training$NV)= paste('TRNV-', colnames(training$NV), sep = '')
training$NV = cbind(training$NV, CN = patient$CNA[, 1])

data = cbind(training$NR, training$NV)
data = cbind(data, TEST = NA)

which.test = c(patient$TES1$NR[, 'M'], patient$TES2$NR[, 'M'])
which.test = tapply(unlist(which.test), names(unlist(which.test)), sum)

data[names(which.test), 'TEST'] = which.test

data = data.frame(data)
data = correctReadCounts(data, 'TRNR.T1', patient$purity['T1', ])
data = correctReadCounts(data, 'TRNR.T2', patient$purity['T2', ])
data = correctReadCounts(data, 'TRNR.T3', patient$purity['T3', ])
data = correctReadCounts(data, 'TRNR.T4', patient$purity['T4', ])
#data = correctReadCounts(data, 'TEST', 0.055)
data = correctReadCounts(data, 'TEST', WORST.PURITY)

data = split(data, f = data$CN)
filter = sapply(data, function(w) all(is.na(w$TEST)) | nrow(w) < MINTR.GROUP.SIZE)
data = data[!filter]

pdf('Patient57-BBfit.pdf')
training.params = lapply(c('T1', 'T2', 'T3', 'T4'),
                         BBMLE,
                         main = 'Patient 57',
                         data = data)
dev.off()
jamPDF('Patient57-BBfit.pdf', out.file = 'Patient57-BBfit.pdf', layout = '2x2')
training.params = Reduce(rbind, training.params)

toTest = lapply(data, function(w) w[!is.na(w$TEST),  'TEST', drop = FALSE])
TestsTable = BBpval(training.params, toTest, PVALUE)

patient$TestsTable = TestsTable
patient$data = data
patient$training.params = training.params

save(patient, file = 'P57.RData')

######################################################################
###################################################################### Patient A44
######################################################################

load('PA44.RData')

training.VCF = read.vcfR('[Data] WES_PASS/NG-8132_A44.mutect2.platypus_PASS.vcf')

training = NULL
training$NV <- extract.gt(training.VCF, element='NV', as.numeric=TRUE)
training$NR <- extract.gt(training.VCF, element='NR', as.numeric=TRUE)
rownames(training$NV) = rownames(training$NR) = paste('chr', rownames(training$NV), sep = '')
head(training$NV)

training$NV = training$NV[, 4:7]
training$NR = training$NR[, 4:7]
colnames(training$NV) = colnames(training$NR) = c('T1', 'T2', 'T3', 'T4') 
training$NV = training$NV[rownames(patient$CCF), ]
training$NR = training$NR[rownames(patient$CCF), ]

rownames(patient$CNA) == rownames(training$NR)

colnames(training$NR) = paste('TRNR-', colnames(training$NR), sep = '')
colnames(training$NV)= paste('TRNV-', colnames(training$NV), sep = '')
training$NV = cbind(training$NV, CN = patient$CNA[, 1])

data = cbind(training$NR, training$NV)
data = cbind(data, TEST = NA)

data[rownames(patient$TES2$NR), 'TEST'] = patient$TES2$NR[, 'M']

data = data.frame(data)
data = correctReadCounts(data, 'TRNR.T1', patient$purity['T1', ])
data = correctReadCounts(data, 'TRNR.T2', patient$purity['T2', ])
data = correctReadCounts(data, 'TRNR.T3', patient$purity['T3', ])
data = correctReadCounts(data, 'TRNR.T4', patient$purity['T4', ])
#data = correctReadCounts(data, 'TEST', 0.055)
data = correctReadCounts(data, 'TEST', WORST.PURITY)

data = split(data, f = data$CN)
filter = sapply(data, function(w) all(is.na(w$TEST)) | nrow(w) < MINTR.GROUP.SIZE)
data = data[!filter]

pdf('PatientA44-BBfit.pdf')
training.params = lapply(c('T1', 'T2', 'T3', 'T4'),
                         BBMLE,
                         main = 'Patient A44',
                         data = data)
dev.off()
jamPDF('PatientA44-BBfit.pdf', out.file = 'PatientA44-BBfit.pdf', layout = '2x2')
training.params = Reduce(rbind, training.params)

toTest = lapply(data, function(w) w[!is.na(w$TEST),  'TEST', drop = FALSE])
TestsTable = BBpval(training.params, toTest, PVALUE)

patient$TestsTable = TestsTable
patient$data = data
patient$training.params = training.params

save(patient, file = 'PA44.RData')

######################################################################
###################################################################### Patient SP28
######################################################################

load('PSP28.RData')

training.VCF = read.vcfR('[Data] WES_PASS/NG-8132_SP28.mutect2.platypus_PASS.vcf')

training = NULL
training$NV <- extract.gt(training.VCF, element='NV', as.numeric=TRUE)
training$NR <- extract.gt(training.VCF, element='NR', as.numeric=TRUE)
rownames(training$NV) = rownames(training$NR) = paste('chr', rownames(training$NV), sep = '')
head(training$NV)

training$NV = training$NV[, c(3,7)]
training$NR = training$NR[, c(3,7)]
colnames(training$NV) = colnames(training$NR) = c('T1', 'T2') 
training$NV = training$NV[rownames(patient$CCF), ]
training$NR = training$NR[rownames(patient$CCF), ]

rownames(patient$CNA) == rownames(training$NR)

colnames(training$NR) = paste('TRNR-', colnames(training$NR), sep = '')
colnames(training$NV)= paste('TRNV-', colnames(training$NV), sep = '')
training$NV = cbind(training$NV, CN = patient$CNA[, 1])

data = cbind(training$NR, training$NV)
data = cbind(data, TEST = NA)

data[rownames(patient$TES2$NR), 'TEST'] = patient$TES2$NR[, 'M1']

data = data.frame(data)
data = correctReadCounts(data, 'TRNR.T1', patient$purity['T1', ])
data = correctReadCounts(data, 'TRNR.T2', patient$purity['T2', ])
#data = correctReadCounts(data, 'TEST', 0.055)
data = correctReadCounts(data, 'TEST', WORST.PURITY)

data = split(data, f = data$CN)
filter = sapply(data, function(w) all(is.na(w$TEST)) | nrow(w) < MINTR.GROUP.SIZE)
data = data[!filter]

pdf('PatientSP28-BBfit.pdf')
training.params = lapply(c('T1', 'T2'),
                         BBMLE,
                         main = 'Patient SP28',
                         data = data)
dev.off()
jamPDF('PatientSP28-BBfit.pdf', out.file = 'PatientSP28-BBfit.pdf', layout = '2x2')
training.params = Reduce(rbind, training.params)

toTest = lapply(data, function(w) w[!is.na(w$TEST),  'TEST', drop = FALSE])
TestsTable = BBpval(training.params, toTest, PVALUE)

patient$TestsTable = TestsTable
patient$data = data
patient$training.params = training.params

save(patient, file = 'PSP28.RData')



