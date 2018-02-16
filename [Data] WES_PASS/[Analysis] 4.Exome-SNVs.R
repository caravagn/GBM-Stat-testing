library(vcfR)

# Takes the input from read.vcfR and gives the index for the SNVs only
getSNVsonly = function(f) {
  
  require(vcfR)
  
  snvlocs = which(!unlist(lapply(strsplit(getREF(f), ""), function(x) length(x)>1)) & !unlist(lapply(strsplit(getALT(f), ""), function(x) length(x)>1)))
  
  return(snvlocs)
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


BEDFILE = '../S04380110_Covered.bed'
patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A44', 'SP28')

for(patient in patients)
{  
  file = paste('NG-8132_', patient, '.mutect2.platypus_PASS.vcf', sep = '')
  WES = read.vcfR(file)
  
  WES.VCF = NULL
  WES.VCF$NV <- extract.gt(WES, element='NV', as.numeric=TRUE)
  WES.VCF$NR <- extract.gt(WES, element='NR', as.numeric=TRUE)
  
  ### RESTRICT TO SNVS and TARGET REGIONS
  SNVsloc = getSNVsonly(WES)
  WES.VCF$NV = WES.VCF$NV[SNVsloc, ]
  WES.VCF$NR = WES.VCF$NR[SNVsloc, ]
  
  Targetsloc = subsetTargetRegions(rownames(WES.VCF$NV), BEDFILE)
  WES.VCF$NV = WES.VCF$NV[Targetsloc, ]
  WES.VCF$NR = WES.VCF$NR[Targetsloc, ]
  
  rownames(WES.VCF$NV) = paste('chr', rownames(WES.VCF$NV), sep = '')
  rownames(WES.VCF$NR) = paste('chr', rownames(WES.VCF$NR), sep = '')
  
  SNVs = rownames(WES.VCF$NR)                    
  
  save(SNVs, file = paste('Exone-SNVs-', patient, '.RData', sep = ''))
}



