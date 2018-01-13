require(pheatmap)
require(RColorBrewer)
library(vcfR)

plotter = function(WES.VCF, selection, out.file)
{
  cwd = getwd()
  setwd(outputPlotsFolder)
  
  bin = WES.VCF$NV
  bin[bin > 0] = 1
  
  myorder = function() {
    scoreCol = function(x) {
      score = 0
      for(i in 1:length(x)) {
        if(x[i]) {
          score = score + 2^(length(x)-i*1/x[i])
        }
      }
      return(score)
    }
    scores = apply(bin, 1, scoreCol)
    order(scores, decreasing=TRUE)
  }
  
  bin = bin[myorder(), , drop = FALSE]
  order.rows = list(rownames(bin))
  save(order.rows, file = paste(out.file, 'order', 'Rdata', sep = '.'))
  
  annotations = NULL
  if(!all(is.null(selection))) annotations = as.data.frame(as.factor(selection))
  
  pheatmap(bin, 
           main = 'Binarized',
           color = c('white', 'darkblue'), 
           # show_rownames = F, 
           annotation_row = annotations,
           cluster_cols = F,
           cluster_rows = F,
           cellwidth = 10, 
           fontsize_row = 3,
           cellheight = 5,
           filename = 'a.pdf')
  
  pheatmap(bin, 
           main = 'Binarized',
           color = c('white', 'darkblue'), 
           show_rownames = F, 
           annotation_row = annotations,
           cluster_cols = F,
           cluster_rows = F,
           cellwidth = 20,
           cellheight = .5,
           filename = 'a2.pdf')
  
  colors = colorRampPalette(brewer.pal(9, 'Blues'))(100)
  mybreaks = c(0, 1:100)
  if(max(WES.VCF$NV) > 100) mybreaks = c(mybreaks, max(WES.VCF$NV))
  
  pheatmap(WES.VCF$NV[rownames(bin), , drop = F], 
           main = 'Reads with the variant allele',
           breaks = mybreaks, 
           color = c('gainsboro', colors), 
           annotation_row = annotations,
           cluster_cols = F,
           cluster_rows = F,
           cellwidth = 10, 
           fontsize_row = 3,
           cellheight = 5,
           filename = 'b.pdf')
  
  pheatmap(WES.VCF$NV[rownames(bin), , drop = F], 
           main = 'Reads with the variant allele',
           breaks = mybreaks, 
           color = c('gainsboro', colors), 
           annotation_row = annotations,
           show_rownames = F, 
           cluster_cols = F,
           cluster_rows = F,
           cellwidth = 20,
           cellheight = .5,
           filename = 'b2.pdf')
  
  # 
  jamPDF(in.files = c('a.pdf','b.pdf'), out.file = paste('zoom-', out.file, sep = ''), layout = '2x1')
  jamPDF(in.files = c('a2.pdf','b2.pdf'), out.file = out.file, layout = '2x1')
  
  setwd(cwd)
}


GIT = '~/Documents/GitHub/GBM-Stat-testing'
outputPlotsFolder = paste(GIT, '/Plots', sep = '')

source('[Analysis] Auxiliary.lib.R')
setwd('[Data] WES_PASS')

######################################################################
###################################################################### Patient 42
######################################################################

WES = read.vcfR('NG-8132_42.mutect2.platypus_PASS.vcf')
WES.VCF = NULL
WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)

selection = apply(WES.VCF$NV, 1, function(w){
  w['NG-8132_42M_lib74072_3832_1'] > 0 &&
    sum(
      w[c("NG-8132_42T1", "NG-8132_42T2_lib74075_3782_7", "NG-8132_42T3_lib74076_3832_1", "NG-8132_42T4_lib74077_3766_1")]
    ) == 0
})

WES.VCF$NV[selection, ]
plotter(WES.VCF, selection, 'WES-42.pdf')

######################################################################
###################################################################### Patient 49
######################################################################

WES = read.vcfR('NG-8132_49.mutect2.platypus_PASS.vcf')
WES.VCF = NULL
WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)

selection = apply(WES.VCF$NV, 1, function(w){
  w['NG-8132_49M_lib74079_3782_1'] > 0 &&
    sum(
      w[c("NG-8132_49T1_lib74081_3847_1", "NG-8132_49T2_lib74082_3832_1", "NG-8132_49T3_lib74083_3847_1", "NG-8132_49T4_lib74084_3782_1")]
    ) == 0
})

WES.VCF$NV[selection, ]
plotter(WES.VCF, selection, 'WES-49.pdf')

######################################################################
###################################################################### Patient 52
######################################################################
WES = read.vcfR('NG-8132_52.mutect2.platypus_PASS.vcf')
WES.VCF = NULL
WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)

selection = apply(WES.VCF$NV, 1, function(w){
  w['NG-8132_52M_lib74086_3847_1'] > 0 &&
    sum(
      w[c("NG-8132_52T1_lib74088_3782_7", "NG-8132_52T2_lib74089_3832_2", "NG-8132_52T3_lib74090_3832_2", "NG-8132_52T4_lib74091_3847_1")]
    ) == 0
})

WES.VCF$NV[selection, ]
plotter(WES.VCF, selection, 'WES-52.pdf')
######################################################################
###################################################################### Patient 54
######################################################################
WES = read.vcfR('NG-8132_54.mutect2.platypus_PASS.vcf')
WES.VCF = NULL
WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)

selection = apply(WES.VCF$NV, 1, function(w){
  w['NG-8132_54M1_lib74093_3847_2'] > 0 && w['NG-8132_54M1_lib74093_3847_2'] > 0
    sum(
      w[c("NG-8132_54T1_lib74095_3786_1", "NG-8132_54T2_lib74096_3847_2", "NG-8132_54T3_lib74097_3782_2", "NG-8132_54T4_lib74098_3847_2",
          "NG-8132_54T5_lib74099_3832_2", "NG-8132_54T6_lib74100_3782_2")]
    ) == 0
})

WES.VCF$NV[selection, ]
plotter(WES.VCF, selection, 'WES-54.pdf')
######################################################################
###################################################################### Patient 55
######################################################################
WES = read.vcfR('NG-8132_55.mutect2.platypus_PASS.vcf')
WES.VCF = NULL
WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)

selection = NULL
plotter(WES.VCF, selection, 'WES-55.pdf')
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
plotter(WES.VCF, selection, 'WES-56.pdf')
######################################################################
###################################################################### Patient 57
######################################################################
WES = read.vcfR('NG-8132_57.mutect2.platypus_PASS.vcf')
WES.VCF = NULL
WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)

selection = apply(WES.VCF$NV, 1, function(w){
  w['NG-8132_57M_lib74115_3782_4'] > 0 &&
    sum(
      w[c("NG-8132_57T1_lib74117_3782_4", "NG-8132_57T2_lib74118_3782_4", "NG-8132_57T3_lib74119_3847_4", "NG-8132_57T4_lib74120_3782_5")]
    ) == 0
})

WES.VCF$NV[selection, ]
plotter(WES.VCF, selection, 'WES-57.pdf')

######################################################################
###################################################################### Patient A23
######################################################################
WES = read.vcfR('NG-8132_A23.mutect2.platypus_PASS.vcf')
WES.VCF = NULL
WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)

selection = apply(WES.VCF$NV, 1, function(w){
  w['NG-8132_A23M_lib74157_3847_5'] > 0 &&  w['NG-8132_SP19M_lib74154_3847_4'] > 0 &&
    sum(
      w[c("NG-8132_A23T_lib74159_3847_5", "NG-8132_SP19T_lib74156_3847_5")]
    ) == 0
})

WES.VCF$NV[selection, ]
plotter(WES.VCF, selection, 'WES-A23.pdf')
######################################################################
###################################################################### Patient A44
######################################################################
WES = read.vcfR('NG-8132_A44.mutect2.platypus_PASS.vcf')
WES.VCF = NULL
WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)

selection = apply(WES.VCF$NV, 1, function(w){
  w['NG-8132_A44M_lib74138_3786_3'] > 0 &&  
    sum(
      w[c("NG-8132_A44T1_lib74132_3786_1", "NG-8132_A44T2_lib74133_3786_1", "NG-8132_A44T3_lib74134_3786_2", "NG-8132_A44T5_lib74136_3786_2")]
    ) == 0
})

WES.VCF$NV[selection, ]
plotter(WES.VCF, selection, 'WES-A44.pdf')
######################################################################
###################################################################### Patient SP28
######################################################################
WES = read.vcfR('NG-8132_SP28.mutect2.platypus_PASS.vcf')
WES.VCF = NULL
WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)

selection = apply(WES.VCF$NV, 1, function(w){
  w['NG-8132_R11M_lib74150_3786_5'] > 0 && w["NG-8132_SP28M_lib74147_3786_6"]
    sum(
      w[c("NG-8132_R11T_lib74152_3786_5", "NG-8132_SP28T_lib74149_3786_5")]
    ) == 0
})

WES.VCF$NV[selection, ]
plotter(WES.VCF, selection, 'WES-SP28.pdf')


setwd(GIT)
