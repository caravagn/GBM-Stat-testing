source('Lib.R')

library(vcfR)
library(pheatmap)
require(gridExtra)

PATIENTS = c('42', '49', '52', '54', '56', '57', 'A44', 'SP28')

purity = 
  list(
    list(TUM = c(S = 0.701, T1 = 0.890, T2 = 0.932, T3 = 0.757, T4 = 0.294), M = 0.081), # 42
    list(TUM = c(S = 0.999, T1 = 0.999, T2 = 0.376, T3 = 0.453, T4 = 0.478), M = 0.055), # 49
    list(TUM = c(S = 0.185, T1 = 0.525, T2 = 0.651, T3 = 0.690, T4 = 0.686), M = 0.032), # 52
    list(TUM = c(T1 = 0.592, T2 = 0.659, T3 = 0.415, T4 = 0.838, T5 = 0.996, T6 = 0.730), M = 0.083), # 54
    # list(TUM = c(0.446, 0.751, 0.304, 0.712), M = 0.498), # 55 mistake in the Excel?
    list(TUM = c(S = 0.257, T1 = 0.594, T2 = 0.108, T3 = 0.404, T4 = 0.712), M = 0.49),  # 56
    list(TUM = c(S = 0.05, T1 = 0.173, T2 = 0.334, T3 = 0.315, T4 = 0.294), M = 0.041), # 57
    # list(TUM = c(S = 0.02, T1 = 0.228), M = 0.010), # A23
    list(TUM = c(S = 0.03, T1 = 0.487, T2 = 0.220, T3 = 0.167, T5 = 0.083), M = 0.038) # A44
    # list(TUM = c(S = 0.04, T1 = 0.585), M = 0.046) # SP28
  )
names(purity) = PATIENTS

WES.FOLDER = 'WES_PASS'
TES1.FOLDER = 'TES_1'
TES2.FOLDER = 'TES_2'

PROCESS = PATIENTS[8]

WES.file = paste(WES.FOLDER, '/NG-8132_', PROCESS, '.mutect2.platypus_PASS.vcf', sep = '') 
TES1.file = paste(TES1.FOLDER, '/patient_', PROCESS, '_platypus.vcf', sep = '') 
TES2.file = paste(TES2.FOLDER, '/', PROCESS, '.platypus.vcf', sep = '') 

dir.create(PROCESS)
sink(file = paste(PROCESS, "/log.txt", sep = ''), split = TRUE)

######### Extract from VCF files the SNVs information
WES.data = loadVCF(WES.file)
TES1.data = loadVCF(TES1.file)
TES2.data = loadVCF(TES2.file)

assampl(WES.data)
assampl(TES1.data)
assampl(TES2.data)

######### Exceptions to deal with
if(PROCESS == '49') colnames(WES.data$NR)[3] = colnames(WES.data$NV)[3] = colnames(WES.data$VAF)[3] = '49S'
if(PROCESS == '52') colnames(WES.data$NR)[3] = colnames(WES.data$NV)[3] = colnames(WES.data$VAF)[3] = '52S'
if(PROCESS == '54') {
  colnames(WES.data$NR)[2] = colnames(WES.data$NV)[2] = colnames(WES.data$VAF)[2] = '54M'
  colnames(TES1.data$NR)[2] = colnames(TES1.data$NV)[2] = colnames(TES1.data$VAF)[2]  = '54M'
  colnames(TES2.data$NR)[1] = colnames(TES2.data$NV)[1] = colnames(TES2.data$VAF)[1]  = '54M'
}
if(PROCESS == '56') {
  colnames(WES.data$NR)[2] = colnames(WES.data$NV)[2] = colnames(WES.data$VAF)[2] = '56M'
  colnames(TES1.data$NR)[2] = colnames(TES1.data$NV)[2] = colnames(TES1.data$VAF)[2]  = '56M'
  colnames(TES2.data$NR)[1] = colnames(TES2.data$NV)[1] = colnames(TES2.data$VAF)[1]  = '56M'
}

######### Panel TES2 contains mutations from all cohort -- we just need the ones for this patient
w = which(rownames(TES2.data$NV) %in% rownames(WES.data$NV))
TES2.data$NV = TES2.data$NV[w, ]  
TES2.data$NR = TES2.data$NR[w, ]  
TES2.data$VAF = TES2.data$VAF[w, ]  

######### Compute adjusted VAF for purity
WES.data$VAF.adj = correctForPurity(WES.data$VAF, purity[[PROCESS]], PROCESS)
TES1.data$VAF.adj = correctForPurity(TES1.data$VAF, purity[[PROCESS]], PROCESS)
TES2.data$VAF.adj = correctForPurity(TES2.data$VAF, purity[[PROCESS]], PROCESS)

######################################################################################################

margin.sample = paste(PROCESS, 'M', sep ='')

# Clonal from the targeted panel 1 -- no correction
TES1.clonal = subset_clonal_mutations(TES1.data, clonal.cutoff = 0.2,  correction = NULL)
print(TES1.clonal)

# Clonal from WES, but corrected according to the targeted panel - higher resolution
WES.clonal = subset_clonal_mutations(WES.data,  clonal.cutoff = 0.2,  correction = TES1.clonal)
print(WES.clonal)

# Zeroes in the margins of the panels, with minimum coverage of 100x
TES1.zeroes = subset_zeroesM_mutations(TES1.data, min.coverage = 100)
TES2.zeroes = subset_zeroesM_mutations(TES2.data, min.coverage = 100)
print(asmuts(TES1.zeroes))
print(asmuts(TES2.zeroes))

TES1.testable = intersect(asmuts(TES1.zeroes), asmuts(WES.clonal))
TES2.testable = intersect(asmuts(TES2.zeroes), asmuts(WES.clonal))
print(TES1.testable)
print(TES2.testable)

# Final test, whatever is testable in one panel, is also testable in the other (if it is present)
TES1.testable = crossCheck(TES1.testable, TES2.data, TES2.zeroes)
TES2.testable = crossCheck(TES2.testable, TES1.data, TES1.zeroes)
print(TES1.testable)
print(TES2.testable)

quartz(width = 10, height = 20)
plot_VAF(
  WES.data,
  TES1.data,
  TES2.data,
  asmuts(WES.clonal),
  TES1.testable,
  TES2.testable,
  clonal.cutoff = 0.2,
  purity = purity[[PROCESS]], 
  PROCESS = PROCESS,
  show.SUBCLONAL = FALSE
  )
dev.copy2pdf(file = paste(PROCESS, '/Data.pdf', sep =''))
dev.off()  

setwd(PROCESS)
res1 = batch_tests(WES.clonal, TES1.zeroes, TES1.testable, margin.sample, purity[[PROCESS]],  psign = 0.05, panel = 'TES1')
res2 = batch_tests(WES.clonal, TES2.zeroes, TES2.testable, margin.sample, purity[[PROCESS]], psign = 0.05,  panel = 'TES2')

write.csv(res1, "Stats-TES1.csv")
write.csv(res2, "Stats-TES2.csv")

library(gridExtra)

quartz(width = 10)
grid.table(res1)
dev.copy2pdf(file = paste('Stats-TES1.pdf', sep =''))
dev.off()

quartz(width = 10)
grid.table(res2)
dev.copy2pdf(file = paste('Stats-TES2.pdf', sep =''))
dev.off()

setwd('..')  

