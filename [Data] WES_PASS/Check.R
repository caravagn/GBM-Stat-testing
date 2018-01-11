require(pheatmap)
require(RColorBrewer)
library(vcfR)

outputPlotsFolder = '.'

plotter = function(TES1.VCF, panel, clonal, cols, order.by, out.file)
{
  cwd = getwd()
  setwd(outputPlotsFolder)
  
  ann_col = c('darkred', 'darkgreen')
  names(ann_col) = c('NO', 'YES')

  order = order(TES1.VCF$NV[, order.by])
  
  TES1 = NULL
  for(c in cols) {
    attach = cbind(TES1.VCF$NV[, c], TES1.VCF$NR[, c])
    attach = attach[order, , drop = F ]
    
    colnames(attach) = c(paste('Var. Reads in ', c), paste('Num. Reads in ', c))
    TES1 = cbind(TES1, attach) 
  }
  
  maxV = max(TES1, na.rm = T)
  mybreaks = c(0, seq(1, maxV, by = maxV/9))
  mycolors = c('white', brewer.pal(9, 'Blues'))
  names(mycolors) = mybreaks
  
  df = data.frame(clonal = matrix('NO', nrow = nrow(TES1), ncol = 1), stringsAsFactors = F)
  rownames(df) = rownames(TES1)
  df[rownames(df) %in% w, ] = 'YES'
                  
  pheatmap(TES1, 
           main = panel,
           na_col = 'gainsboro', 
           color = mycolors,
           breaks = mybreaks,
           cluster_rows = FALSE, 
           annotation_row = df,
           annotation_colors = list(clonal = ann_col),
           cluster_cols = FALSE, 
           show_rownames = TRUE,
           fontsize_row = 4,
           gaps_col = ifelse(length(cols) > 1, seq(2, length(cols), by = 2), 0),
           display_numbers = T,
           border_color = NA, 
           number_format = '%d',
           number_color = 'orange',
           cellwidth = 20,
           fontsize = 6,
           cellheight = 5,
           file = paste('Fig3', panel, out.file, sep = '-'))
  
  setwd(cwd)
  # jamPDF(in.files = c('mut.pdf','cna.pdf'), out.file = out.file, layout = '1x1')
}

setwd('../CCFs/')

GIT = '~/Documents/GitHub/GBM-Stat-testing'
CCF.FOLDER = paste(GIT, '/CCFs/', sep = '')
TES1.FOLDER = paste(GIT, '/TES_1/', sep = '')
TES2.FOLDER = paste(GIT, '/TES_2/', sep = '')


######################################################################
###################################################################### Patient 42
######################################################################
setwd(CCF.FOLDER)
CCF = read.csv('Patient_42_cancer_cell_fractions.csv', header = TRUE)
rownames(CCF) = paste('chr', rownames(CCF), sep = '')
CCF = rownames(CCF)

setwd(outputPlotsFolder)
w = read.table('WES-42.pdf-samples.txt', stringsAsFactors = FALSE)[, 1]

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_42_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2, 7)]
TES1.VCF$NR = TES1.VCF$NR[, c(2, 7)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('S', 'M')

setwd(TES2.FOLDER)

TES2 = read.vcfR('42.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M', 'S')
TES2.VCF$NV = TES2.VCF$NV[rownames(TES2.VCF$NV) %in% CCF, , drop = F]
TES2.VCF$NR = TES2.VCF$NR[rownames(TES2.VCF$NR) %in% CCF, , drop = F]

plotter(TES1.VCF, 'TES1', w,  c('M', 'S'), 'M', 'P42.pdf')
plotter(TES2.VCF, 'TES2', w,  c('M', 'S'), 'M', 'P42.pdf')

######################################################################
###################################################################### Patient 49
######################################################################

WES = read.vcfR('NG-8132_49.mutect2.platypus_PASS.vcf')
WES.VCF = NULL
WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)

selection = apply(WES.VCF$NV, 1, function(w){
  w['NG-8132_56M3_lib74108_3847_3'] > 0 &&
    sum(
      w[c("NG-8132_56T1_lib74110_3847_3", "NG-8132_56T2_lib74111_3832_4", "NG-8132_56T3_lib74112_3782_3", "NG-8132_56T4_lib74113_3832_4")]
    ) == 0
})

WES.VCF$NV[selection, ]
write.vcf(WES[selection, ], file = '56-weird.csv.gz')

colors = colorRampPalette(brewer.pal(9, 'Blues'))(100)
pheatmap(WES.VCF$NV, breaks = c(0, 1:100, max(WES.VCF$NV)), color = c('gainsboro', colors), show_rownames = F, cluster_cols = F)


######################################################################
###################################################################### Patient 52
######################################################################
setwd(CCF.FOLDER)
CCF = read.csv('Patient_52_cancer_cell_fractions.csv', header = TRUE)
rownames(CCF) = paste('chr', rownames(CCF), sep = '')
CCF = rownames(CCF)

setwd(outputPlotsFolder)
w = read.table('WES-52.pdf-samples.txt', stringsAsFactors = FALSE)[, 1]

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_52_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2), drop = F]
TES1.VCF$NR = TES1.VCF$NR[, c(2), drop = F]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M')

setwd(TES2.FOLDER)

TES2 = read.vcfR('52.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M', 'S')

TES2.VCF$NV = TES2.VCF$NV[rownames(TES2.VCF$NV) %in% CCF, , drop = F]
TES2.VCF$NR = TES2.VCF$NR[rownames(TES2.VCF$NR) %in% CCF, , drop = F]

plotter(TES1.VCF, 'TES1', w,  c('M'), 'M', 'P52.pdf')
plotter(TES2.VCF, 'TES2', w,  c('M', 'S'), 'M', 'P52.pdf')

######################################################################
###################################################################### Patient 54
######################################################################
setwd(CCF.FOLDER)
CCF = read.csv('Patient_54_cancer_cell_fractions.csv', header = TRUE)
rownames(CCF) = paste('chr', rownames(CCF), sep = '')
CCF = rownames(CCF)

setwd(outputPlotsFolder)
w = read.table('WES-54.pdf-samples.txt', stringsAsFactors = FALSE)[, 1]

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_54_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2,9)]
TES1.VCF$NR = TES1.VCF$NR[, c(2,9)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M', 'S')

setwd(TES2.FOLDER)

TES2 = read.vcfR('54.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M1', 'M2', 'S')

TES2.VCF$NV = TES2.VCF$NV[rownames(TES2.VCF$NV) %in% CCF, , drop = F]
TES2.VCF$NR = TES2.VCF$NR[rownames(TES2.VCF$NR) %in% CCF, , drop = F]

plotter(TES1.VCF, 'TES1', w,  c('M', 'S'), 'M', 'P54.pdf')
plotter(TES2.VCF, 'TES2', w,  c('M1', 'M2', 'S'), 'M1', 'P54.pdf')

######################################################################
###################################################################### Patient 55
######################################################################
setwd(CCF.FOLDER)
CCF = read.csv('Patient_55_cancer_cell_fractions.csv', header = TRUE)
rownames(CCF) = paste('chr', rownames(CCF), sep = '')
CCF = rownames(CCF)

setwd(outputPlotsFolder)
w = read.table('WES-55.pdf-samples.txt', stringsAsFactors = FALSE)[, 1]

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_55_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2,3)]
TES1.VCF$NR = TES1.VCF$NR[, c(2,3)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M', 'S')

setwd(TES2.FOLDER)

TES2 = read.vcfR('55.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('S')

TES2.VCF$NV = TES2.VCF$NV[rownames(TES2.VCF$NV) %in% CCF, , drop = F]
TES2.VCF$NR = TES2.VCF$NR[rownames(TES2.VCF$NR) %in% CCF, , drop = F]

plotter(TES1.VCF, 'TES1', w,  c('M', 'S'), 'M', 'P55.pdf')
plotter(TES2.VCF, 'TES2', w,  c('S'), 'S', 'P55.pdf')


######################################################################
###################################################################### Patient 56
######################################################################

WES = read.vcfR('NG-8132_56.mutect2.platypus_PASS.vcf')
WES.VCF = NULL
WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)


heatmap.bp(WES.VCF$NV, rlabels = F, col.ramp =  colorRampPalette(c("white", "orange", "red"))(10))

selection = apply(WES.VCF$NV, 1, function(w){
  w['NG-8132_56M3_lib74108_3847_3'] > 0 &&
    sum(
      w[c("NG-8132_56T1_lib74110_3847_3", "NG-8132_56T2_lib74111_3832_4", "NG-8132_56T3_lib74112_3782_3", "NG-8132_56T4_lib74113_3832_4")]
    ) == 0
})

WES.VCF$NV[selection, ]
write.vcf(WES[selection, ], file = '56-weird.csv.gz')

colors = colorRampPalette(brewer.pal(9, 'Blues'))(100)
pheatmap(WES.VCF$NV, breaks = c(0, 1:100, max(WES.VCF$NV)), color = c('gainsboro', colors), show_rownames = F, cluster_cols = F)


######################################################################
###################################################################### Patient 57
######################################################################
setwd(CCF.FOLDER)
CCF = read.csv('Patient_57_cancer_cell_fractions.csv', header = TRUE)
rownames(CCF) = paste('chr', rownames(CCF), sep = '')
CCF = rownames(CCF)

setwd(outputPlotsFolder)
w = read.table('WES-57.pdf-samples.txt', stringsAsFactors = FALSE)[, 1]

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_57_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2,3)]
TES1.VCF$NR = TES1.VCF$NR[, c(2,3)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M', 'S')

setwd(TES2.FOLDER)

TES2 = read.vcfR('57.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

TES2.VCF$NV = TES2.VCF$NV[, c(1,2), drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[, 1:2, drop = FALSE]
colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M', 'S')

TES2.VCF$NV = TES2.VCF$NV[rownames(TES2.VCF$NV) %in% CCF, , drop = F]
TES2.VCF$NR = TES2.VCF$NR[rownames(TES2.VCF$NR) %in% CCF, , drop = F]

plotter(TES1.VCF, 'TES1', w,  c('M', 'S'), 'M', 'P57.pdf')
plotter(TES2.VCF, 'TES2', w,  c('M', 'S'), 'M', 'P57.pdf')
######################################################################
###################################################################### Patient A23
######################################################################
setwd(CCF.FOLDER)
CCF = read.csv('Patient_A23_cancer_cell_fractions.csv', header = TRUE)
rownames(CCF) = paste('chr', rownames(CCF), sep = '')
CCF = rownames(CCF)

setwd(outputPlotsFolder)
w = read.table('WES-A23.pdf-samples.txt', stringsAsFactors = FALSE)[, 1]

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_SP19_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2,4,5)]
TES1.VCF$NR = TES1.VCF$NR[, c(2,4,5)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M1',  'M2', 'S')

setwd(TES2.FOLDER)

TES2 = read.vcfR('A23_SP19.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

TES2.VCF$NV = TES2.VCF$NV[, 1:4, drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[, 1:4, drop = FALSE]
colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M1', 'S1', 'M2', 'S2')

TES2.VCF$NV = TES2.VCF$NV[rownames(TES2.VCF$NV) %in% CCF, , drop = F]
TES2.VCF$NR = TES2.VCF$NR[rownames(TES2.VCF$NR) %in% CCF, , drop = F]

plotter(TES1.VCF, 'TES1', w,  c('M1', 'M2', 'S'), 'M1', 'PA23.pdf')
plotter(TES2.VCF, 'TES2', w,  c('M1', 'S1', 'M2', 'S2'), 'M1', 'P57.pdf')

######################################################################
###################################################################### Patient A34
######################################################################

setwd(CCF.FOLDER)
CCF = read.csv('Patient_A34_cancer_cell_fractions.csv', header = TRUE)
rownames(CCF) = paste('chr', rownames(CCF), sep = '')
CCF = rownames(CCF)

setwd(outputPlotsFolder)
w = read.table('WES-A34.pdf-samples.txt', stringsAsFactors = FALSE)[, 1]

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_A34_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, 1:3]
TES1.VCF$NR = TES1.VCF$NR[, 1:3]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('S1',  'S2', 'S3')

plotter(TES2.VCF, 'TES2', w,  c('S1', 'S2', 'S3'), 'S1', 'PA34.pdf')

######################################################################
###################################################################### Patient A44
######################################################################
setwd(CCF.FOLDER)
CCF = read.csv('Patient_A44_cancer_cell_fractions.csv', header = TRUE)
rownames(CCF) = paste('chr', rownames(CCF), sep = '')
CCF = rownames(CCF)

setwd(outputPlotsFolder)
w = read.table('WES-A44.pdf-samples.txt', stringsAsFactors = FALSE)[, 1]

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_A44_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(2,3)]
TES1.VCF$NR = TES1.VCF$NR[, c(2,3)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M', 'S')

setwd(TES2.FOLDER)
TES2 = read.vcfR('A44.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

TES2.VCF$NV = TES2.VCF$NV[, c(1,2), drop = FALSE]
TES2.VCF$NR = TES2.VCF$NR[, 1:2, drop = FALSE]
colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('M', 'S')

TES2.VCF$NV = TES2.VCF$NV[rownames(TES2.VCF$NV) %in% CCF, , drop = F]
TES2.VCF$NR = TES2.VCF$NR[rownames(TES2.VCF$NR) %in% CCF, , drop = F]

plotter(TES1.VCF, 'TES1', w,  c('M', 'S'), 'M', 'PA44.pdf')
plotter(TES2.VCF, 'TES2', w,  c('M', 'S'), 'M', 'PA44.pdf')
######################################################################
###################################################################### Patient SP28
######################################################################
setwd(CCF.FOLDER)
CCF = read.csv('Patient_SP28_cancer_cell_fractions.csv', header = TRUE)
rownames(CCF) = paste('chr', rownames(CCF), sep = '')
CCF = rownames(CCF)

setwd(outputPlotsFolder)
w = read.table('WES-SP28.pdf-samples.txt', stringsAsFactors = FALSE)[, 1]

setwd(TES1.FOLDER)
TES1 = read.vcfR('patient_SP28_platypus.vcf')

TES1.VCF = NULL
TES1.VCF$NV <- extract.gt(TES1, element='NV', as.numeric=TRUE)
TES1.VCF$NR <- extract.gt(TES1, element='NR', as.numeric=TRUE)

TES1.VCF$NV = TES1.VCF$NV[, c(1,2,5,6)]
TES1.VCF$NR = TES1.VCF$NR[, c(1,2,5,6)]
colnames(TES1.VCF$NV) = colnames(TES1.VCF$NR) = c('M1', 'S1', 'M2', 'S2')

setwd(TES2.FOLDER)
TES2 = read.vcfR('SP28_R11.platypus.vcf')

TES2.VCF = NULL
TES2.VCF$NV <- extract.gt(TES2, element='NV', as.numeric=TRUE)
TES2.VCF$NR <- extract.gt(TES2, element='NR', as.numeric=TRUE)
rownames(TES2.VCF$NV) = paste('chr', rownames(TES2.VCF$NV), sep = '')
rownames(TES2.VCF$NR) = paste('chr', rownames(TES2.VCF$NR), sep = '')

colnames(TES2.VCF$NV) = colnames(TES2.VCF$NR) = c('S1', 'M1', 'M2', 'S2')

TES2.VCF$NV = TES2.VCF$NV[rownames(TES2.VCF$NV) %in% CCF, , drop = F]
TES2.VCF$NR = TES2.VCF$NR[rownames(TES2.VCF$NR) %in% CCF, , drop = F]

plotter(TES1.VCF, 'TES1', w,  c('M1', 'S1', 'M2', 'S2'), 'M1', 'PSP28.pdf')
plotter(TES2.VCF, 'TES2', w,  c('M1', 'S1', 'M2', 'S2'), 'M1', 'PSP28.pdf')
