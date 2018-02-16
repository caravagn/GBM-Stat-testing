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

patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A34', 'A44', 'SP28')


for(patient in patients)
{
  load(paste('CCF-', patient, '.RData', sep = ''), verbose = T)
  
  
  
  
  load(paste('CNA-', patient, '.RData', sep = ''))
  
  clonal = primary(CCF)
  clonal = clonal[complete.cases(clonal), ]
  
  subclonal = rownames(clonal[apply(clonal, 1, function(x) any(x < CLONALITY_CUTOFF)), ])
  clonal = clonal[apply(clonal, 1, function(x) all(x >= CLONALITY_CUTOFF)), ]
  
  CNA = primary(CNA)
  CNA = CNA[complete.cases(CNA), ]
  CNA = CNA[apply(CNA, 1, function(x) length(unique(x)) == 1), ]
  
  clonalSNVs = rownames(clonal)[rownames(clonal) %in% rownames(CNA)]
  clonal = cbind(clonal[clonalSNVs, ], `CNA` = CNA[clonalSNVs, 1])

  save(clonal, file = paste('CLONAL-', patient, '.RData', sep = ''))
  save(subclonal, file = paste('SUBCLONAL-', patient, '.RData', sep = ''))
}


