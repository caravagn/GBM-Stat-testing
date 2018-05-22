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

# Takes the input from read.vcfR and gives the index for the SNVs only
getSNVsonly = function(f) {
  
  require(vcfR)
  
  snvlocs = which(!unlist(lapply(strsplit(getREF(f), ""), function(x) length(x)>1)) & !unlist(lapply(strsplit(getALT(f), ""), function(x) length(x)>1)))
  
  return(snvlocs)
}

getnonSNV = function(f){
 
  snvlocs = which(!unlist(lapply(strsplit(getREF(f), ""), function(x) length(x)>1)) & !unlist(lapply(strsplit(getALT(f), ""), function(x) length(x)>1)))
  setdiff(1:nrow(f), snvlocs)
  
  # snvlocs %in% all.rows
  # all.rows[!(all.rows %in% snvlocs)]
}

# Subsets regions that overlap with a bed file, takes positions as input as "1_123456", also location of bed file
subsetTargetRegions = function(pos, bed_file) {
  
  #Need GenomicRanges to do this
  require(GenomicRanges)
  
  #Read in the exome bed file
  exome.bait = read.table(bed_file, skip = 2)
  
  #Split the pos
  pos_split = strsplit(pos, "_")
  
  #Write function for extracting chr and pos
  return_item = function(x, item) {unlist(lapply(x, function(x) x[item]))}
  
  #What are the chrs and locs
  chrs = return_item(pos_split, item = 1)
  locs = as.numeric(return_item(pos_split, item = 2))
  
  #Create genomic ranges for data and target regions
  vafs.granges = GRanges(seqnames = paste0("chr",chrs), 
                         IRanges(start = locs, 
                                 end = locs))
  exom.granges = GRanges(seqnames = exome.bait$V1, IRanges(start = exome.bait$V2, end = exome.bait$V3))
  
  #What is the overlap?
  overlap_index = unique(queryHits(findOverlaps(vafs.granges, exom.granges)))
  
  #Return 
  return(overlap_index)
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

replace_S_SVZ = function(w) {
  S = substr(w, 1, 1)
  S2 = substr(w, 2, 2)
  
  wt = substr(w, 2, nchar(w))
  
  if(substr(w, 1, 3) != 'SVZ' && S == 'S' && S2 != 'V') 
    return(paste('SVZ', wt, sep = ''))
  else
    w
}

BEDFILE = '../S04380110_Covered.bed'

# Iavarone and Rabadan
CSQ.specific.GBM = c('ATRX', 'TP53', 'MMR', 'LTBP4', 'PIK3CA', 'PIK3R1', 'PDGFRA', 'EGFR', 'NF1', 'PTPN11', 'PTEN', 'RB1', 'TERT')

# TCGA
CSQ.specific.GBM = c('PIK3R1', 'PIK3CA', 'PTEN', 'RB1', 'TP53', 'EGFR', 'IDH1', 'BRAF',
                     'NF1', 'SPTA1', 'GABRA6', 'KEL', 'CDH18', 'SEMA3C', 'PDGFRA', 'ATRX',
                     'COL1A2', 'LZTR1', 'ABCC9', 'NLRP5', 'DRD5', 'TCHH', 'SCN9A')

patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A34', 'A44', 'SP28')

for(patient in patients)
{  
  file = paste('NG-8132_', patient, '.mutect2.platypus_PASS.vcf', sep = '')
  WES = read.vcfR(file)
  
  WES.VCF = NULL
  WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
  WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)
  
  ### RESTRICT TO SNVS and TARGET REGIONS
  loc = getnonSNV(WES)
  locSNV = getSNVsonly(WES)
  intersect(loc, locSNV)
  
  annotation = data.frame(Type = rep('SNV', nrow(WES.VCF$NV)), stringsAsFactors = FALSE)
  annotation[loc, ] = 'Indel'
  rownames(annotation) = rownames(WES.VCF$NV)
  
  names.to.fix = which(!startsWith(rownames(annotation), 'chr'))
  rownames(annotation)[names.to.fix] = paste('chr', rownames(annotation)[names.to.fix], sep = '')
  
  # Data
  data = WES.VCF$NV/WES.VCF$NR
  dataB = data
  dataB[dataB > 0.05] = 1
  dataB[dataB <= 0.05] = 0
  
  names.to.fix = which(!startsWith(rownames(dataB), 'chr'))
  rownames(dataB)[names.to.fix] = paste('chr', rownames(dataB)[names.to.fix], sep = '')
  
  ### ANNOTATED CCF
  load(paste('../[Data] CCFs Annotated/ANNOTATED-CGC-', patient, '.RData', sep = ''), verbose = T)
  
  CSQ[!(CSQ$CGC %in% CSQ.specific.GBM), 'CGC'] = ''
  head(CSQ)

  rownames(CSQ) %in% rownames(dataB)
  
  dataB = dataB[swantonOrder(dataB), , drop = F]
  head(dataB)
  
  CSQ = CSQ[rownames(dataB), , drop = FALSE]
  annotation = annotation[rownames(dataB), , drop = FALSE]
  rownames(CSQ) == rownames(dataB)
  rownames(annotation) == rownames(dataB)
  
  
  idx = which(CSQ$CGC != '')
  for(i in idx) CSQ[i, 'CGC'] = paste(CSQ[i, 'CGC'], annotation[i, ], CSQ[i, 'Mutation'])
  
  head(dataB)
  dataB = namify.wesPanel(dataB, patient, 'WES')
  head(dataB)
  
  colnames(dataB) = sapply(colnames(dataB), replace_S_SVZ)
  
  if(patient == '55') colnames(dataB)[3] = 'T1-WES'
  if(patient == 'A23') {
    colnames(dataB) = c('B-WES', 'M recurrent-WES', 
                                              'SVZ recurrent-WES', 'T1 recurrent-WES', 'M primary-WES',
                                              'SVZ primary-WES', 'T1 primary-WES')
    dataB = dataB[, c(1, 5,6,7, 2:4)]
  }
  if(patient == 'A34') colnames(dataB)[1] = 'B-WES'
  if(patient == 'SP28') {
    colnames(dataB) = c('M recurrent-WES', 'SVZ recurrent-WES', 'T1 recurrent-WES',
                        'B-WES', 'M primary-WES', 'SVZ primary-WES', 'T1 primary-WES')
    dataB = dataB[, c(4, 5,6,7, 1:3)]
  }
  
  CSQ$CGC = paste(rownames(CSQ), ' ', CSQ$CGC)
  
  
  pheatmap(dataB, 
           main = paste(patient, '-- VAF cutoff 0.05'),
           cluster_rows = F, 
           cluster_cols = F,
           annotation_row = annotation[rownames(dataB), , drop = F],
           color = RColorBrewer::brewer.pal(9, 'Blues'),
           legend = F,
           cellwidth = 20,
           cellheight = 8,
           fontsize_row = 6,
           labels_row = CSQ$CGC,
           file = paste(patient, '-indels_SNV-plot.pdf', sep = ''))
  
  WES.VCF$NV = WES.VCF$NV[loc, ]
  WES.VCF$NR = WES.VCF$NR[loc, ]
  
  CSQ = CSQ[loc, ]
  colnames(CSQ) = c('Annotated gene', 'GBM driver')
   
  rownames(WES.VCF$NV) = paste('chr', rownames(WES.VCF$NV), sep = '')
  rownames(WES.VCF$NR) = paste('chr', rownames(WES.VCF$NR), sep = '')
  
  indels = data.frame(WES.VCF$NV/WES.VCF$NR, stringsAsFactors = FALSE)
  indels = cbind(indels, CSQ)
  
  WriteXLS::WriteXLS(indels, 
                     ExcelFileName = paste('Patient-', patient, '-indels.xlsx', sep = ''), 
                     row.names = TRUE)
  
  save(indels, file = paste('Exone-non-SNVs-', patient, '.RData', sep = ''))
}



