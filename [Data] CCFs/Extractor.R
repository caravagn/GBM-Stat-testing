files = list.files(pattern = '*')
files = files[endsWith(files, 'cancer_cell_fractions.csv')]

for(file in files)
{  
  CCF = read.csv(file, header = TRUE)
  rownames(CCF) = paste('chr', rownames(CCF), sep = '')
  CCF[apply(CCF, 2, is.infinite)] = 0
  
  save(CCF, file = paste(file, '.RData', sep = ''))
}



