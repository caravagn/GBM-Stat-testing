

for(patient in patients)
{
  load(paste('CCF-', patient, '.RData', sep = ''))
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


