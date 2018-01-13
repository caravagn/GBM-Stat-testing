require(pheatmap)
require(RColorBrewer)


plotter = function(CCF, tum.cols, out.file, gaps_col = 2)
{
  cwd = getwd()
  setwd(outputPlotsFolder)
  
  CCF[CCF > 1] = 1
  # CCF = CCF[
  #   order(rowSums(CCF[, tum.cols, drop = F]), na.last = TRUE, decreasing = TRUE), , drop = FALSE]
  
  mybreaks = unique(c(0, 0.0001, seq(0.0001, 1.1, by = .1)))
  mycolor = colorRampPalette(brewer.pal(9, 'YlGnBu')) (length(mybreaks))
  mycolor = c('white', mycolor)
  
  pheatmap(CCF, 
           main = 'WES SNVs',
           na_col = 'gainsboro', 
           color = mycolor,
           breaks = mybreaks,
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           show_rownames = FALSE,
           border_color = NA, 
           gaps_col = gaps_col,
           cellwidth = 20,
           cellheight = .5,
           file = paste('CCF-ALL-', out.file, sep = ''))
  
  CCF = CCF[apply(CCF[, tum.cols, drop = F], 1, function(w) !any(is.na(w))), ]
  
  # pheatmap(CCF, 
  #          main = 'WES \nClonal CNA',
  #          na_col = 'gainsboro', 
  #          color = mycolor,
  #          breaks = mybreaks,
  #          cluster_rows = FALSE, 
  #          cluster_cols = FALSE, 
  #          show_rownames = TRUE,
  #          fontsize_row = 4,
  #          border_color = NA,
  #          cellwidth = 20,
  #          cellheight = 5,
  #          file = 'clonal-CNA.pdf')
  
  CCF = CCF[apply(CCF[, tum.cols, drop = F], 1, function(w) all(w > 0.8)), ]  
  
  pheatmap(CCF, 
           main = 'WES \n Clonal SNVs',
           na_col = 'gainsboro', 
           color = mycolor,
           breaks = mybreaks,
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           show_rownames = TRUE, 
           fontsize_row = 4,
           cellheight = 5,
           border_color = NA,
           gaps_col = gaps_col,
           cellwidth = 20,
           file = paste('CCF-Clonal', out.file, sep = ''))
  
  print(paste(out.file, '-samples.txt', sep = ''))
  write.table(rownames(CCF), paste(out.file, '-samples.txt', sep = ''))
  
  setwd(cwd)
}

GIT = '~/Documents/GitHub/GBM-Stat-testing'
outputPlotsFolder = paste(GIT, '/Plots', sep = '')


setwd('[Data] CCFs/')
######################################################################
###################################################################### Patient 42
######################################################################
CCF = read.csv('Patient_42_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

colnames(CCF) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

# CCF = CCF[, c('T1', 'T2', 'T3', 'T4')]

load('../Plots/WES-42.pdf.order.RData')
order.rows = paste('chr', unlist(order.rows), sep = '')
CCF = CCF[order.rows, ]

plotter(CCF, tum.cols = c('T1', 'T2', 'T3', 'T4'), out.file = 'WES-42.pdf')
######################################################################
###################################################################### Patient 49
######################################################################
CCF = read.csv('Patient_49_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

colnames(CCF) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

# CCF = CCF[, c('T1', 'T2', 'T3', 'T4')]

load('../Plots/WES-49.pdf.order.RData')
order.rows = paste('chr', unlist(order.rows), sep = '')
CCF = CCF[order.rows, ]

plotter(CCF, tum.cols = c('T1', 'T2', 'T3', 'T4'), out.file = 'WES-49.pdf')

######################################################################
###################################################################### Patient 52
######################################################################
CCF = read.csv('Patient_52_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

colnames(CCF) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

load('../Plots/WES-52.pdf.order.RData')
order.rows = paste('chr', unlist(order.rows), sep = '')
CCF = CCF[order.rows, ]

plotter(CCF, tum.cols = c('T1', 'T2', 'T3', 'T4'), 'WES-52.pdf')

######################################################################
###################################################################### Patient 54
######################################################################

CCF = read.csv('Patient_54_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

colnames(CCF) = c('M', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

load('../Plots/WES-54.pdf.order.RData')
order.rows = paste('chr', unlist(order.rows), sep = '')
CCF = CCF[order.rows, ]

plotter(CCF, tum.cols = c('T1', 'T2', 'T3', 'T4', 'T5', 'T6'), 'WES-54.pdf', gaps_col = 1)

######################################################################
###################################################################### Patient 55
######################################################################

CCF = read.csv('Patient_55_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

colnames(CCF) = c('S', 'T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

load('../Plots/WES-55.pdf.order.RData')
order.rows = paste('chr', unlist(order.rows), sep = '')
CCF = CCF[order.rows, ]

plotter(CCF, tum.cols = c('T1', 'T2', 'T3', 'T4'), 'WES-55.pdf', gaps_col = 1)

######################################################################
###################################################################### Patient 56
######################################################################

CCF = read.csv('Patient_56_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

colnames(CCF) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

load('../Plots/WES-56.pdf.order.RData')
order.rows = paste('chr', unlist(order.rows), sep = '')
CCF = CCF[order.rows, ]

plotter(CCF, tum.cols = c('T1', 'T2', 'T3', 'T4'), 'WES-56.pdf')

######################################################################
###################################################################### Patient 57
######################################################################

CCF = read.csv('Patient_57_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

colnames(CCF) = c('M', 'S', 'T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

load('../Plots/WES-57.pdf.order.RData')
order.rows = paste('chr', unlist(order.rows), sep = '')
CCF = CCF[order.rows, ]

plotter(CCF, tum.cols = c('T1', 'T2', 'T3', 'T4'), 'WES-57.pdf')

######################################################################
###################################################################### Patient A23
######################################################################

CCF = read.csv('Patient_A23_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

colnames(CCF) = c('M1', 'S1', 'T1', 'M2', 'S2', 'T2')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')
CCF = CCF[, c('M1', 'S1',  'M2', 'S2', 'T1', 'T2')]

load('../Plots/WES-A23.pdf.order.RData')
order.rows = paste('chr', unlist(order.rows), sep = '')
CCF = CCF[order.rows, ]

plotter(CCF, tum.cols = c('T1', 'T2'), 'WES-A23.pdf', gaps_col = c(2, 4))

######################################################################
###################################################################### Patient A34
######################################################################

CCF = read.csv('Patient_A34_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

colnames(CCF) = c('S1', 'S2', 'S3', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

plotter(CCF, tum.cols = c('T1', 'T2', 'T3', 'T4', 'T5', 'T6'), 'WES-A34.pdf', gaps_col = 3)

######################################################################
###################################################################### Patient A44
######################################################################

CCF = read.csv('Patient_A44_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

colnames(CCF) = c("M", 'S', 'T1', 'T2', 'T3', 'T4')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

load('../Plots/WES-A44.pdf.order.RData')
order.rows = paste('chr', unlist(order.rows), sep = '')
CCF = CCF[order.rows, ]

plotter(CCF, tum.cols = c('T1', 'T2', 'T3', 'T4'), 'WES-A44.pdf')

######################################################################
###################################################################### Patient SP28
######################################################################

CCF = read.csv('Patient_SP28_cancer_cell_fractions.csv', header = TRUE)
head(CCF)

colnames(CCF) = c('M1', 'S1', 'T1', 'M2', 'S2', 'T2')
rownames(CCF) = paste('chr', rownames(CCF), sep = '')

CCF = CCF[, c('M1', 'S1',  'M2', 'S2', 'T1', 'T2')]

load('../Plots/WES-SP28.pdf.order.RData')
order.rows = paste('chr', unlist(order.rows), sep = '')
CCF = CCF[order.rows, ]

plotter(CCF, tum.cols = c('T1', 'T2'), 'WES-SP28.pdf', gaps_col = c(2, 4))
