source('Lib.R')

library(vcfR)
library(pheatmap)
require(gridExtra)

WES.file = '49/NG-8132_49.mutect2.platypus_PASS.vcf'
TES1.file = '49/patient_49_platypus.vcf'
TES2.file = '49/49.platypus.vcf'

######### Extract from VCF files the SNVs information
WES.data = loadVCF(WES.file)
TES1.data = loadVCF(TES1.file)
TES2.data = loadVCF(TES2.file)

assampl(WES.data)
assampl(TES1.data)
assampl(TES2.data)

colnames(WES.data$NR)[3] = colnames(WES.data$NV)[3] = '49S'


# Panel TES2 contains mutations from all cohort -- we just need the ones for this patient
w = which(rownames(TES2.data$NV) %in% rownames(WES.data$NV))
TES2.data$NV = TES2.data$NV[w, ]  
TES2.data$NR = TES2.data$NR[w, ]  

######################################################################################################

margin.sample = '49M'

# Clonal from the targeted panel
TES1.clonal = subset_clonal_mutations(TES1.data, type = 'clonal')

# Clonal from WES, but corrected according to the targeted panel - higher resolution
WES.clonal = subset_clonal_mutations(WES.data, type = 'clonal', correction = TES1.clonal)

# Zeroes in the margins of the panels, with minimum coverage of 100x
TES1.zeroes = subset_zeroesM_mutations(TES1.data, min.coverage = 100)
TES2.zeroes = subset_zeroesM_mutations(TES2.data, min.coverage = 100)

TES1.testable = intersect(asmuts(TES1.zeroes), asmuts(WES.clonal))
TES2.testable = intersect(asmuts(TES2.zeroes), asmuts(WES.clonal))

# ######### Show some data 
# quartz(width = 20, height = 30)
# grid.arrange(
#   plot_VAF(WES.data, main = 'WES patient 42', mut.annot = training, 'Training set'),
#   plot_VAF(TES1.data, main = 'TES1 patient 42', mut.annot = data$test.TES1, 'Test set'),
#   plot_VAF(TES2.data, main = 'TES2 patient 42', mut.annot = data$test.TES2, 'Test set'),
#   ncol = 3
# )
# dev.copy2pdf(file = 'Data.pdf')
# dev.off()
# 
# plot_training_set(WES.data, 'WES-training-set.pdf')
# 

tested.sample = '42T1'

# For the test, we need first a training model. We train on SNVs from the tested.sample that are not
# the ones that we are going to need for the pvalue. These variants are
purity = list(TUM = c(.87, .82, .67, .31), M = .34)
purity = list(TUM = c(0.890, 0.932, 0.757, 0.294), M = 0.081)


# 49 Mr Bayes
purity = list(TUM = c(0.999, 0.376, 0.453, 0.478), M = 0.055)





batch_tests(WES.clonal, TES1.zeroes, TES1.testable, margin.sample, purity)
batch_tests(WES.clonal, TES2.zeroes, TES2.testable, margin.sample, purity)
  

