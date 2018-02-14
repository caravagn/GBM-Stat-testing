require(pheatmap)
require(RColorBrewer)
library(vcfR)

CCF.plot = function(CCF, patient)
{
  is.tmr = function(CCF) 
  {
    t = grepl('T', colnames(CCF))
    t1 = grepl('T1', colnames(CCF))
    t2 = grepl('T2', colnames(CCF))
    t3 = grepl('T3', colnames(CCF))
    t4 = grepl('T4', colnames(CCF))
    t5 = grepl('T5', colnames(CCF))
    t6 = grepl('T6', colnames(CCF))
    
    df = rbind(t, t1, t2, t3, t4, t5, t6)
    CCF[, apply(df, 2, any), drop = FALSE]
  }
  
  CCF = is.tmr(CCF)
  CCF = CCF[myorder(CCF), , drop = F]
  
  cwd = getwd()
  setwd(outputPlotsFolder)
  
  pheatmap(CCF, 
           main = patient,
           na_col = 'darkgray', 
           color = c('gainsboro', 'steelblue'),
           # breaks = mybreaks,
           cluster_rows = FALSE, 
           # annotation_row = df,
           # annotation_colors = list(clonal = ann_col),
           cluster_cols = FALSE, 
           show_rownames = FALSE,
           fontsize_row = 4,
           legend = FALSE,
           gaps_col = 1:ncol(CCF),
           border_color = NA, 
           cellwidth = 20,
           fontsize = 6,
           cellheight = .5,
           file = paste('FinalPlot', patient, '.pdf', sep = '-'))
  
  setwd(cwd)
}

VAF.plot = function(CCF, patient)
{
  is.tmr = function(CCF) 
  {
    t = grepl('T', colnames(CCF))
    t1 = grepl('T1', colnames(CCF))
    t2 = grepl('T2', colnames(CCF))
    t3 = grepl('T3', colnames(CCF))
    t4 = grepl('T4', colnames(CCF))
    t5 = grepl('T5', colnames(CCF))
    t6 = grepl('T6', colnames(CCF))
    
    df = rbind(t, t1, t2, t3, t4, t5, t6)
    CCF[, apply(df, 2, any), drop = FALSE]
  }
  
  CCF = is.tmr(CCF)
  CCF = CCF[myorder(CCF), , drop = F]
  
  cwd = getwd()
  setwd(outputPlotsFolder)
  
  pheatmap(CCF, 
           main = patient,
           na_col = 'darkgray', 
           color = c('gainsboro', 'steelblue'),
           # breaks = mybreaks,
           cluster_rows = FALSE, 
           # annotation_row = df,
           # annotation_colors = list(clonal = ann_col),
           cluster_cols = FALSE, 
           show_rownames = FALSE,
           fontsize_row = 4,
           legend = FALSE,
           gaps_col = 1:ncol(CCF),
           border_color = NA, 
           cellwidth = 20,
           fontsize = 6,
           cellheight = .5,
           file = paste('FinalPlot-VAF', patient, '.pdf', sep = '-'))
  
  setwd(cwd)
}

myorder = function(data) {
  data[is.na(data)] = 0
  
  private = rowSums(data)
  
  datap = data[names(private[private == 1]), , drop = F]
  datanp = data[names(private[private != 1]), , drop = F]
  
  scoreCol = function(x) {
    score = 0
    for(i in 1:length(x)) {
      if(x[i]) {
        score = score + 2^(length(x)-i*1/x[i])
      }
    }
    return(score)
  }
  
  scoresnp = apply(datanp, 1, scoreCol)
  onp = order(c(scoresnp), decreasing=TRUE)
  
  scoresp = apply(datap, 1, scoreCol) 
  op = order(c(scoresp), decreasing=TRUE)
  
  
  fo = c(onp, op + max(onp))
  fo  
  
  order(apply(data, 1, scoreCol), decreasing=TRUE)
}


GIT = '~/Documents/GitHub/GBM-Stat-testing'
outputPlotsFolder = paste(GIT, '/Plots', sep = '')

CCF.FOLDER = paste(GIT, '/[Data] CCFs/', sep = '')
WES.FOLDER = paste(GIT, '/[Data] WES_PASS/', sep = '')
# TES1.FOLDER = paste(GIT, '/[Data] TES_1/', sep = '')
# TES2.FOLDER = paste(GIT, '/[Data] TES_2/', sep = '')
# 

#################### FIRST PLOT
setwd(CCF.FOLDER)

files = list.files()
files = files[endsWith(files, '.RData')]
files = files[startsWith(files, 'CCF')]

for(f in files) 
{
  load(f, verbose = T)
  
  pname = strsplit(f, split = '\\.')[[1]][1]
  
  CCF[CCF < 0.2] = 0
  CCF[CCF >= 0.2] = 1
  
  CCF.plot(CCF, pname)
}

#################### FIRST PLOT with VAF isntead of CCF
setwd(WES.FOLDER)

files = list.files()
files = files[endsWith(files, '.RData')]
files = files[startsWith(files, 'VAF')]

for(f in files) 
{
  load(f, verbose = T)
  
  pname = strsplit(f, split = '\\.')[[1]][1]
  
  VAF[VAF < 0.05] = 0
  VAF[VAF >= 0.05] = 1
  
  VAF.plot(VAF, pname)
}


#################### SECOND PLOT

ms = function(x) 
{
  m = grepl('M', colnames(x))
  s = grepl('S', colnames(x))

  df = rbind(m, s)
  x[, apply(df, 2, any), drop = FALSE]
}

######### P52
# 
# setwd(WES.FOLDER)
# load('NG-8132_52.mutect2.platypus_PASS.vcf.RData', verbose = T)
# 
# WES.VCF = ms(WES.VCF$NV)
# WES.VCF[WES.VCF > 0] = 1
# colnames(WES.VCF) = c('M', 'S')
# 
# setwd(TES1.FOLDER)
# load('patient_52_platypus.vcf.RData', verbose = T)
# 
# TES1.VCF = ms(TES1.VCF$NV)[, 2, drop = FALSE]
# TES1.VCF[TES1.VCF > 0] = 1
# colnames(TES1.VCF) = 'M'
# 
# setwd(TES2.FOLDER)
# load('52.platypus.vcf.RData', verbose = T)
# 
# TES2.VCF = ms(TES2.VCF$NV)
# TES2.VCF[TES2.VCF > 0] = 1
# colnames(TES2.VCF) = c('M', 'S')
# 
# load('../P52.RData')
# 
# block_A_M = function(WES.VCF, TES1.VCF, TES2.VCF){
#   x = WES.VCF[, 'M']
#   y = TES1.VCF[, 'M']
#   z = TES2.VCF[, 'M']
#   unique(names(c(x[x == 1], y[y == 1], z[z == 1])))
# }
# 
# bAM = block_A_M(WES.VCF, TES1.VCF, TES2.VCF)
# bCM = rownames(patient$TestsTable[patient$TestsTable$sign, ])
# bDM = rownames(patient$TestsTable[!patient$TestsTable$sign, ])
# bBM = setdiff(rownames(WES.VCF), c(bAM, bCM, bDM))
# 
# bCM %in% bAM
# 
# c(bAM, bCM, bDM)
# 
# rownames(WES.VCF)
# 
# m = c(bAM, bBM, bCM, bDM)
# duplicated(m)
# 
# mat = data.frame(row.names = c(bAM, bBM, bCM, bDM))
# df = data.frame()
# 
# TES2.VCF[bCM, ]
# 
# 
# block_A_M = names(WES.VCF[WES.VCF[, 'M'] == 1, 'M'])
# block_A_M = c(block_A_M, )
# 
# 
# 
