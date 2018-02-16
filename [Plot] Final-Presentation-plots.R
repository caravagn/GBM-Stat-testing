require(pheatmap)
require(RColorBrewer)
library(vcfR)

swantonOrder = function(data) {
  data[is.na(data)] = 0
  
  scoreCol = function(x) {
    score = 0
    for(i in 1:length(x)) {
      if(x[i]) {
        score = score + 2^(length(x)-i*1/x[i])
      }
    }
    return(score)
  }
  
  sharedData  = data[which(apply(data, 1, sum) > 1),]
  sharedIndex = which(apply(data, 1, sum) > 1)
  
  scores   = apply(sharedData, 1, scoreCol)
  topOrder = sharedIndex[order(scores, decreasing=TRUE)]
  
  privateData  = data[which(apply(data, 1, sum) <= 1),]
  privateIndex = which(apply(data, 1, sum) <= 1)
  
  scores   = apply(privateData, 1, scoreCol)
  bottomOrder = privateIndex[order(scores, decreasing=TRUE)]
  
  return(c(topOrder, bottomOrder))
}

is.tmr = function(CCF) {
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


CCF.plot = function(CCF, patient, annotation = NULL)
{

  cwd = getwd()
  setwd(outputPlotsFolder)
  
  col = c(`YES` = 'forestgreen', `NO` = 'brown4', `Not Testable` = 'lightgray')
  col = list(`Significant P-Value` = col)
  
  pheatmap(CCF, 
           main = patient,
           na_col = 'darkgray', 
           color = c('gainsboro', 'steelblue'),
           # breaks = mybreaks,
           cluster_rows = FALSE, 
           annotation_row = annotation,
           annotation_colors = col,
           cluster_cols = FALSE, 
           show_rownames = TRUE,
           fontsize_row = 4,
           legend = FALSE,
           gaps_col = 1:ncol(CCF),
           border_color = NA, 
           cellwidth = 20,
           fontsize = 6,
           cellheight = 4,
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
  CCF = CCF[swantonOrder(CCF), , drop = F]
  
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

ms = function(x) 
{
  m = grepl('M', colnames(x))
  s = grepl('S', colnames(x))
  
  df = rbind(m, s)
  x[, apply(df, 2, any), drop = FALSE]
}


subsetPanel = function(CCF, PANEL) {
  
  PANEL = PANEL[rownames(PANEL) %in% rownames(CCF), , drop = FALSE]
  
  df = data.frame(matrix(0, nrow = nrow(CCF), ncol = ncol(PANEL)))
  rownames(df) = rownames(CCF)
  colnames(df) = colnames(PANEL)
  
  for(i in 1:ncol(PANEL)) df[rownames(PANEL), i] = PANEL[, i]
  df
}


GIT = '~/Documents/GitHub/GBM-Stat-testing'
outputPlotsFolder = paste(GIT, '/Plots', sep = '')

CCF.FOLDER = paste(GIT, '/[Data] CCFs/', sep = '')
WES.FOLDER = paste(GIT, '/[Data] WES_PASS/', sep = '')
TES1.FOLDER = paste(GIT, '/[Data] TES_1/', sep = '')
TES2.FOLDER = paste(GIT, '/[Data] TES_2/', sep = '')
# 

#################### FIRST PLOT
setwd(CCF.FOLDER)

files = list.files()
files = files[endsWith(files, '.RData')]
files = files[startsWith(files, 'CCF')]

files = files[files !=  "CCF-A34.RData" ]

CCF.CUTOFF = 0.2
NV.CUTOFF = 1

for(f in files) 
{
  load(f, verbose = T)
  
  CCF[CCF < CCF.CUTOFF] = 0
  CCF[CCF >= CCF.CUTOFF] = 1
  
  CCF = is.tmr(CCF)
  
  pname = strsplit(f, split = '\\.')[[1]][1]
  patient = strsplit(pname, split = '-')[[1]][2]
  
  # Get the list of real SNVs (no indels) that are in exome regions, etc.
  load(paste(WES.FOLDER, '/Exone-SNVs-', patient, '.RData', sep = ''), verbose = TRUE)
  CCF = CCF[SNVs, , drop = FALSE]
  
  # Get read counts from Margin and S -- WES
  load(paste(WES.FOLDER, '/MARGIN_S_READCOUNTS-', patient, '.RData', sep = ''), verbose = TRUE)
  WES = WES$NV
  WES[WES < NV.CUTOFF] = 0
  WES[WES >= NV.CUTOFF] = 1
  
  # Get read counts from Margin and S -- TES1
  load(paste(TES1.FOLDER, '/MARGIN_S_READCOUNTS-', patient, '.RData', sep = ''), verbose = TRUE)
  TES1 = TES1$NV
  TES1[TES1 < NV.CUTOFF] = 0
  TES1[TES1 >= NV.CUTOFF] = 1
  TES1 = subsetPanel(CCF, TES1)
  
  # Get read counts from Margin and S -- TES1
  load(paste(TES2.FOLDER, '/MARGIN_S_READCOUNTS-', patient, '.RData', sep = ''), verbose = TRUE)
  TES2 = TES2$NV
  TES2[TES2 < NV.CUTOFF] = 0
  TES2[TES2 >= NV.CUTOFF] = 1
  TES2 = subsetPanel(CCF, TES2)
  
  CCF = CCF[swantonOrder(CCF), , drop = F]
  CCF = cbind(CCF, WES[rownames(CCF), ])
  CCF = cbind(CCF, TES1[rownames(CCF), ])
  CCF = cbind(CCF, TES2[rownames(CCF), ])
  
  # Test pass/ non pasas
  results.file = paste('../RESULTS_TEST-', patient, '.RData', sep = '')
  annotation = NULL
  
  if(file.exists(results.file)) {
    load(results.file, verbose = TRUE)
    T.Summary = Summary[sapply(Summary, function(w) all(w$sign))]
    F.Summary = Summary[sapply(Summary, function(w) any(!w$sign))]
    
    annotation = CCF[, 1, drop = FALSE]
    colnames(annotation) = 'Significant P-Value'
    annotation[TRUE] = 'Not Testable'
    annotation[names(T.Summary), 1] = 'YES'
    annotation[names(F.Summary), 1] = 'NO'
    
  }
  
  CCF.plot(CCF, pname, annotation)
}


#################### FIRST PLOT with VAF isntead of CCF
# setwd(WES.FOLDER)
# 
# files = list.files()
# files = files[endsWith(files, '.RData')]
# files = files[startsWith(files, 'VAF')]
# 
# files = files[files !=  "CCF-A34.RData" ]
# 
# for(f in files) 
# {
#   load(f, verbose = T)
#   
#   pname = strsplit(f, split = '\\.')[[1]][1]
#   patient = strsplit(pname, split = '-')[[1]][2]
#   
#   load(paste(WES.FOLDER, '/Exone-SNVs-', patient, '.RData', sep = ''), verbose = TRUE)
#   VAF = VAF[SNVs, , drop = FALSE]
#   
#   VAF[VAF < 0.05] = 0
#   VAF[VAF >= 0.05] = 1
#   
#   VAF.plot(VAF, pname)
# }
