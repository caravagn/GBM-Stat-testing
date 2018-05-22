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


CCF.plot = function(CCF, patient, annotation = NULL, CSQ)
{
  cwd = getwd()
  setwd(outputPlotsFolder)
  
  col = list(`Significant P-Value` = c(`YES` = 'gold', `NO` = 'slateblue4', `Not Testable` = 'whitesmoke'),
             `SNV Status` = c(`ubiquitous` = 'dodgerblue4', `shared` = 'goldenrod4', `private` = 'firebrick', `missing` = 'gainsboro', `NA` = 'darkgray')
             )
  
  CCF.values = CCF[, endsWith(colnames(CCF), 'CCF-WES')]
  ncols = ncol(CCF.values)
  CCF.values = rowSums(CCF.values)
  SNV.status = data.frame(row.names = rownames(CCF), stringsAsFactors = FALSE)
  SNV.status[names(which(CCF.values == ncols)), 'SNV Status'] = 'ubiquitous'
  SNV.status[names(which(CCF.values != ncols & CCF.values > 1)), 'SNV Status'] = 'shared'
  SNV.status[names(which(CCF.values == 1)), 'SNV Status'] = 'private'
  SNV.status[names(which(CCF.values == 0)), 'SNV Status'] = 'missing'
  SNV.status[is.na(SNV.status), 'SNV Status'] = 'NA'
  

  # Iavarone and Rabadan
  CSQ.specific.GBM = c('ATRX', 'TP53', 'MMR', 'LTBP4', 'PIK3CA', 'PIK3R1', 'PDGFRA', 'EGFR', 'NF1', 'PTPN11', 'PTEN', 'RB1', 'TERT')
  
  # TCGA
  CSQ.specific.GBM = c('PIK3R1', 'PIK3CA', 'PTEN', 'RB1', 'TP53', 'EGFR', 'IDH1', 'BRAF',
                            'NF1', 'SPTA1', 'GABRA6', 'KEL', 'CDH18', 'SEMA3C', 'PDGFRA', 'ATRX',
                            'COL1A2', 'LZTR1', 'ABCC9', 'NLRP5', 'DRD5', 'TCHH', 'SCN9A')
  
  
  CSQ.names = CSQ[rownames(CCF), 'CGC']
  CSQ.names = sapply(CSQ.names, function(s) 
    if(nchar(s) > 0 & (s %in% CSQ.specific.GBM)) {paste('-', s);} else "")
  names(CSQ.names) = NULL
  
  
  # CSQ.names = paste(rownames(CCF), CSQ.names)
  
  # !(rownames(annotation) %in%   rownames(SNV.status)) 

  if(is.null(annotation)) {
    annotation.null = data.frame(`Significant P-Value` = rep('Not Testable', nrow(SNV.status)), stringsAsFactors = FALSE)
    rownames(annotation.null) = rownames(SNV.status)
    colnames(annotation.null) = "Significant P-Value"
    
    annotation = cbind(SNV.status, annotation.null)
    # annotation = SNV.status
  } else  annotation = cbind(SNV.status, annotation[rownames(SNV.status), , drop = FALSE])
  
  # annotation = cbind(annotation, CGC = CSQ[rownames(annotation), 'CGC'])
  
  pheatmap::pheatmap(CCF, 
           main = patient,
           na_col = 'darkgray', 
           color = c('gainsboro', 'steelblue', 'darkcyan'),
           # breaks = mybreaks,
           cluster_rows = FALSE, 
           annotation_row = annotation,
           annotation_colors = col,
           cluster_cols = FALSE, 
           show_rownames = TRUE,
           fontsize_row = 10,
           legend = FALSE, 
           labels_row = CSQ.names,
           # legend_labels = c('Missing', 'CCF > .20', 'Binarized Read NV'),
           gaps_col = 1:ncol(CCF),
           border_color = NA, 
           cellwidth = 20,
           fontsize = 16,
           # cellheight = 10,
           cellheight = 4,
           file = paste('FinalPlot', patient, '.pdf', sep = '-'), 
           width = 20, height = 40 
           )

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

### Rename all S to SVZ
replace_S_SVZ = function(w) {
  S = substr(w, 1, 1)
  S2 = substr(w, 2, 2)
  
  wt = substr(w, 2, nchar(w))
  
  if(substr(w, 1, 3) != 'SVZ' && S == 'S' && S2 != 'V') 
    return(paste('SVZ', wt, sep = ''))
  else
    w
}

subsetPanel = function(CCF, PANEL) {
  
  PANEL = PANEL[rownames(PANEL) %in% rownames(CCF), , drop = FALSE]
  
  df = data.frame(matrix(0, nrow = nrow(CCF), ncol = ncol(PANEL)))
  rownames(df) = rownames(CCF)
  colnames(df) = colnames(PANEL)
  
  for(i in 1:ncol(PANEL)) df[rownames(PANEL), i] = PANEL[, i]

  df[!(rownames(CCF) %in% rownames(PANEL)), ] = NA
  
  df
}


namify.wesPanel = function(x, patient, code)
{
  for(cl in 1:ncol(x)) 
  {
    tke = strsplit(colnames(x)[cl], '_')[[1]]
    tke = tke[grepl(patient, tke)]
    tke = gsub(x = tke, patient, '')
    # colnames(x)[cl] = paste(patient, tke, code, sep ='-')
    colnames(x)[cl] = paste(tke, code, sep ='-')
  }
  x
}

namify.TES1Panel = function(x, patient, prefix, code)
{
  chpat = nchar(patient)
  chprf = nchar(prefix)
  
  for(cl in 1:ncol(x)) 
  {
    tke = substr(colnames(x)[cl], 1 + chpat + chprf, nchar(colnames(x)[cl]))
    # colnames(x)[cl] = paste(patient, tke, code, sep ='-')
    colnames(x)[cl] = paste(tke, code, sep ='-')
  }
  x
}

# Remove Blood, order as T -> S -> M
order.columns = function(x)
{
  cn = colnames(x)

  # leftmost
  CCF.values = cn[endsWith(colnames(x), 'CCF-WES')]
  B.values = cn[startsWith(colnames(x), 'B')]
  Margin.values = cn[startsWith(colnames(x), 'M')]
  S.values = cn[startsWith(colnames(x), 'S')]
  rightmost = setdiff(cn, c(CCF.values, B.values, Margin.values, S.values))
  
  Tpaneles.values = cn[startsWith(colnames(x), 'T') & (endsWith(colnames(x), '-TES1') | endsWith(colnames(x), '-TES2'))]
  
  
  # x[, c(CCF.values, B.values, Margin.values, S.values, rightmost)]
  x[, c(CCF.values, Tpaneles.values, S.values, Margin.values)]
}

dumpCoverage = function(panel, id){
  meanCovTES1 = apply(panel, 2, mean)
  medianCovTES1 = apply(panel, 2, median)
  
  file = paste('../Coverage-', id, '.txt', sep = '')
  
  cat(paste('\n', patient, 'mean coverage\n'), file = file, append = TRUE, sep = "\n")
  cat(
    paste(names(meanCovTES1), meanCovTES1), 
    file = file, append = TRUE, sep = "\n")
  
  cat(paste('\n', patient, 'median coverage\n'), file = file, append = TRUE, sep = "\n")
  cat(
    paste(names(medianCovTES1), medianCovTES1), 
    file = file, append = TRUE, sep = "\n")
  
  cat(paste('\n', patient, 'mean coverage across all samples\n'), file = file, append = TRUE, sep = "\n")
  cat(mean(meanCovTES1), 
      file = file, append = TRUE, sep = "\n")
  
  cat(paste('\n', patient, 'median coverage across all samples\n'), file = file, append = TRUE, sep = "\n")
  cat(mean(medianCovTES1), 
      file = file, append = TRUE, sep = "\n")
}


GIT = '~/Documents/GitHub/GBM-Stat-testing'
outputPlotsFolder = paste(GIT, '/Plots', sep = '')

CCF.FOLDER = paste(GIT, '/[Data] CCFs/', sep = '')
WES.FOLDER = paste(GIT, '/[Data] WES_PASS/', sep = '')
TES1.FOLDER = paste(GIT, '/[Data] TES_1/', sep = '')
TES2.FOLDER = paste(GIT, '/[Data] TES_2/', sep = '')
ANNOTATIONS.FOLDER = paste(GIT, '/[Data] CCFs Annotated', sep = '')

#################### FIRST PLOT
setwd(CCF.FOLDER)

files = list.files()
files = files[endsWith(files, '.RData')]
files = files[startsWith(files, 'CCF')]

# files = files[files !=  "CCF-A34.RData" ]

CCF.CUTOFF = 0.2
NV.CUTOFF = 2
TNV.CUTOFF = 10
TES.VAF.CUTOFF = 0.001


for(f in files) 
{
  ########################################## CCF values
  load(f, verbose = T)
  CCF.rownames = rownames(CCF)
  pname = strsplit(f, split = '\\.')[[1]][1]
  patient = strsplit(pname, split = '-')[[1]][2]
  
  CCF[CCF < CCF.CUTOFF] = 0
  CCF[CCF >= CCF.CUTOFF] = 1
  head(CCF)
  
  CCF = is.tmr(CCF)
  head(CCF)
  # colnames(CCF) = paste(patient, colnames(CCF), 'CCF from WES', sep ='-')
  colnames(CCF) = paste(colnames(CCF), 'CCF-WES', sep ='-')
  
  # if(patient == 'A23') colnames(CCF) = c('M(p)-TES2', 'S(p)-TES2', 'M(r)-TES2',  'S(r)-TES2', 'T(r)-TES2')
  if(patient == 'A23') {
    colnames(CCF) = c('T recurrent-CCF-WES', 'T primary-CCF-WES')
    CCF = CCF[, c(2, 1)]
  }
  if(patient == 'SP28') {
    colnames(CCF) = c('T recurrent-CCF-WES', 'T primary-CCF-WES')
    CCF = CCF[, c(2, 1)]
  }
  
  ########################################## Get the list of real SNVs (no indels) that are in exome regions, etc.
  load(paste(WES.FOLDER, '/Exone-SNVs-', patient, '.RData', sep = ''), verbose = TRUE)
  CCF = CCF[SNVs, , drop = FALSE]
  
  ########################################## Get read counts from Margin and S -- WES
  load(paste(WES.FOLDER, '/MARGIN_S_READCOUNTS-', patient, '.RData', sep = ''), verbose = TRUE)
  WES = WES$NV
  WES[WES < NV.CUTOFF] = 0
  WES[WES >= NV.CUTOFF] = 2
  head(WES)
  
  WES = namify.wesPanel(WES, patient, 'WES')
  head(WES)
  
  if(patient == '56') colnames(WES)[2] = 'M-WES'
  if(patient == 'A23') {
    colnames(WES)[2:5] = c('M recurrent-WES', 'S recurrent-WES', 'M primary-WES', 'S primary-WES')
    WES = WES[, c('B-WES', 'M primary-WES', 'S primary-WES', 'M recurrent-WES', 'S recurrent-WES')] 
  }
  if(patient == 'SP28') {
    colnames(WES) = c('M recurrent-WES', 'S recurrent-WES', 'B-WES', 'M primary-WES', 'S primary-WES')
    WES = WES[, c('B-WES', 'M primary-WES', 'S primary-WES', 'M recurrent-WES', 'S recurrent-WES')]
  }
  
  ########################################## Get read counts from Margin and S -- TES1
  load(paste(TES1.FOLDER, '/MARGIN_S_READCOUNTS-', patient, '.RData', sep = ''), verbose = TRUE)

  # cutoff based on VAF
  # TES1.VAF = TES1$NV/TES1$NR
  # TES1.VAF.0 = which(TES1.VAF < TES.VAF.CUTOFF)
  # TES1.VAF.1 = which(TES1.VAF >= TES.VAF.CUTOFF)
  # TES1 = TES1$NV
  # TES1[TES1.VAF.0] = 0
  # TES1[TES1.VAF.1] = 2
  
  dumpCoverage(TES1$NR, 'TES1')

  # cutoff based on read-counts
  TES1 = TES1$NV
  TES1[TES1 < TNV.CUTOFF] = 0
  TES1[TES1 >= TNV.CUTOFF] = 2
   
  TES1 = subsetPanel(CCF, TES1)
  head(TES1)
  
  TES1 = namify.TES1Panel(TES1, patient, 'SP', 'TES1')
  head(TES1)
  
  if(patient == 'A34') colnames(TES1) = c('S1-TES1', 'S2-TES1', 'S3-TES1', 'T1-TES1',  'T2-TES1', 'T3-TES1','T5-TES1', 'T6-TES1')
  if(patient == 'A23') {
    colnames(TES1) = c('B-TES1', 'M recurrent-TES1', 'T recurrent-TES1', 'M primary-TES1',  'S primary-TES1', 'T primary-TES1')
    TES1 = TES1[, c('B-TES1',  'M primary-TES1',  'S primary-TES1', 'T primary-TES1', 'M recurrent-TES1', 'T recurrent-TES1')]
  }
  if(patient == '56') colnames(TES1)[2] = 'M-TES1'
  if(patient == 'A44')colnames(TES1) = c('B-TES1', 'M-TES1', 'S-TES1', 'T1-TES1', 'T2-TES1','T3-TES1', 'T5-TES1')
  if(patient == 'SP28') {
    colnames(TES1) = c('M recurrent-TES1', 'S recurrent-TES1', 'T recurrent-TES1', 'B-TES1', 'M primary-TES1', 'S primary-TES1', 'T primary-TES1')
    TES1 = TES1[, c('M primary-TES1', 'S primary-TES1', 'T primary-TES1', 'M recurrent-TES1', 'S recurrent-TES1', 'T recurrent-TES1', 'B-TES1')]
  }
  
  ########################################## Get read counts from Margin and S -- TES1
  load(paste(TES2.FOLDER, '/MARGIN_S_READCOUNTS-', patient, '.RData', sep = ''), verbose = TRUE)
  
  # cutoff based on VAF
  # TES2.VAF = TES2$NV/TES2$NR
  # TES2.VAF.0 = which(TES2.VAF < TES.VAF.CUTOFF)
  # TES2.VAF.1 = which(TES2.VAF >= TES.VAF.CUTOFF)
  # TES2 = TES2$NV
  # TES2[TES2.VAF.0] = 0
  # TES2[TES2.VAF.1] = 2
  
  dumpCoverage(TES2$NR, 'TES2')
  
  # cutoff based on read-counts
  TES2 = TES2$NV
  TES2[TES2 < TNV.CUTOFF] = 0
  TES2[TES2 >= TNV.CUTOFF] = 2
  
  TES2 = subsetPanel(CCF, TES2)
  head(TES2)

  if(patient == '56') colnames(TES2)[1] = '56M'
  if(patient == '55') colnames(TES2) = '55S'
  
  TES2 = namify.TES1Panel(TES2, patient, '', 'TES2')
  head(TES2)
  
  if(patient == 'A23') colnames(TES2) = c('M primary-TES2', 'S primary-TES2', 'M recurrent-TES2',  'S recurrent-TES2', 'T recurrent-TES2')
  if(patient == 'SP28') colnames(TES2) = c('S primary-TES2', 'M primary-TES2', 'M recurrent-TES2', 'S recurrent-TES2')
  # if(patient == 'A34') colnames(TES2) = c('T2-TES2', 'T3-TES2')
  
  ########################################## Ordering
  CCF = CCF[swantonOrder(CCF), , drop = F]
  CCF = cbind(CCF, WES[rownames(CCF), , drop = FALSE])
  CCF = cbind(CCF, TES1[rownames(CCF), , drop = FALSE])
  CCF = cbind(CCF, TES2[rownames(CCF), , drop = FALSE])
  head(CCF)
  
  ########################################## Annotations
  CSQ = NULL
  if(file.exists(paste(ANNOTATIONS.FOLDER, '/ANNOTATED-CGC-', patient, '.RData', sep = ''))) {
    load(paste(ANNOTATIONS.FOLDER, '/ANNOTATED-CGC-', patient, '.RData', sep = ''), verbose = T)
    rownames(CSQ) = CCF.rownames
    CSG = CSQ[SNVs, ]
    head(CSQ)
  }
  
  # rownames(CCF) = paste(CSQ$CGC, CCF.rownames)
  
  ########################################## Test pass/ non pasas
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
  
  
  
  ########################################## Plot order, rename + plot
  CCF = order.columns(CCF)

  ### Rename all S to SVZ
  colnames(CCF) = sapply(colnames(CCF), replace_S_SVZ)
  
  CCF.plot(CCF, paste('PATIENT', patient), annotation, CSQ)
  
  
  ### Export for phylogenetic inference
  print(head(CCF))
  print(head(annotation))
  CCF[CCF > 1] = 1
  if(!is.null(annotation)) CCF = cbind(CCF, annotation[rownames(CCF), , drop = F])
  write.csv(CCF, file = paste('../Phylogeny-Table-PATIENT', patient,'.txt', sep = ''))
}





###### EXAMPLE POWER-PLOT
# setwd(GIT)
# mu = 0.5
# rho = .05
# 
# x = seq(10, 500, by = 10)
# y = NULL
# for(c in x)
#   y = c(y, sum(dbetabinom(0:10, size = c, prob = mu, rho = rho)))
# 
# # plot(x, y, log = 'xy', type = 'l')
# 
# plot(log(x), log(y), type = 'l', xaxt = 'n', yaxt = 'n', xlab = 'Coverage (adjusted for purity)', ylab = 'P-value')
# points(log(x), log(y), pch = 18, col = 'orange')
# 
# axis(1, x, at = log(x))
# 
# vals = c(1, 0.05, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12) 
# axis(2, vals, at = log(vals))
# abline(h = log(0.05), col = 'red', lty = 2)
# title(bquote(bold('Test power for') ~ mu ~'= 0.5 and'~ rho~ '=5 x'~10^{-2} ~ italic('at significance level')~ alpha ~' = 0.05'))
#       # sub = bquote(bold('Significance level:'~ alpha ~' = 0.05'))
#       # )
# 
# load('RESULTS_TEST-52.RData', verbose = T)
# Summary
# 
# real.points = unlist(lapply(Summary, function(w) w$coverage[1]))
# y = NULL
# for(c in real.points)
#   y = c(y, sum(dbetabinom(0:10, size = c, prob = mu, rho = rho)))
# 
# red.p = log(y) > log(0.05)
# points(log(real.points)[red.p], log(y)[red.p], pch = 19, col = 'red', cex = 2)
# points(log(real.points)[!red.p], log(y)[!red.p], pch = 19, col = 'darkgreen', cex = 2)
# 
# dev.copy2pdf(file = 'test.pdf')


