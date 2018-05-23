library(vcfR)

# COSMIC = read.csv('CosmicCompleteTargetedScreensMutantExport.tsv', header = TRUE, sep = '\t', stringsAsFactors = FALSE)
# head(COSMIC)
# save(COSMIC, file = 'COSMIC.RData')

CGC = read.csv('Census_allWed Mar 28 10-44-54 2018.csv', header = TRUE, sep = ',')
head(CGC)

# 
patients = c('42', '49', '52', '54', '55', '56', '57', 'A23', 'A34', 'A44', 'SP28')

for(patient in patients)
{  
  file = paste('anot_NG-8132_', patient, '.mutect2.platypus_PASS.vcf', sep = '')
  WES = read.vcfR(file)
  
  WES.NV = extract.gt(WES, element='NV', as.numeric=TRUE)
  head(WES.NV)
  
  names.to.fix = which(!startsWith(rownames(WES.NV), 'chr'))
  rownames(WES.NV)[names.to.fix] = paste('chr', rownames(WES.NV)[names.to.fix], sep = '')
  head(WES.NV)
  
  # nrow(WES.NV)
  # coords = rownames(WES.NV)
  # 
  # extract_info_tidy(WES, info_fields = NULL, info_types = TRUE,
  #                   info_sep = ";")
  # 
  # 
  CSQ = extract.info(WES, "CSQ", as.numeric = FALSE, mask = FALSE)
  CSQ = lapply(CSQ, function(w) {
    w = strsplit(w, ',')[[1]]
    w = strsplit(w, '\\|')
    
    if(length(w) > 1) w = w[sapply(w, function(z) z[2] != "upstream_gene_variant")]
    if(length(w) > 1) w = w[sapply(w, function(z) z[2] != "downstream_gene_variant")]
    if(length(w) > 1) paste(unique(sapply(w, function(z) z[4])), collapse = ':')
    else return("unknown")
    
    paste(unique(sapply(w, function(z) paste(z[4], z[2], collapse =' ') )), collapse = ':')
  })
  names(CSQ) = NULL
  CSQ[CSQ == ""] = "unknown"
  head(CSQ)
  
  # 
  # 
  # df = vcfR2tidy(WES)
  # df.COSMIC = NULL
  # 
  # for(i in 1:nrow(df$fix)){
  #   entry = as.character(df$fix[i, 'CSQ'])
  #   entry = strsplit(entry, '\\|')[[1]]
  #   
  #   CHROM = as.character(df$fix[i, 'CHROM'])
  #   POS = as.character(df$fix[i, 'POS'])
  #   REF = as.character(df$fix[i, 'REF'])
  #   ALT = as.character(df$fix[i, 'ALT'])
  #   
  #   if("MODIFIER" %in% entry) {
  #     which.Modifier = min(which("MODIFIER" == entry))
  #     Hugo_Symbol = entry[which.Modifier + 1]
  #     ENSG = entry[which.Modifier + 2]
  #     
  #     df.COSMIC = rbind(df.COSMIC, 
  #                       data.frame(Hugo_Symbol = Hugo_Symbol, ENSG = ENSG,
  #                                  CHROM = CHROM, POS = POS, REF = REF, ALT = ALT,
  #                         stringsAsFactors = FALSE))
  #   } 
  # }
  # 
  # df.COSMIC[df.COSMIC == ""] = "Unknown"
  # 
  # w = which(df.COSMIC$Hugo_Symbol %in% CGC$Gene.Symbol)
  # nrow(df.COSMIC)
  
  CSQ = data.frame(Full = unlist(CSQ), stringsAsFactors = FALSE)
  CSQ$CGC = NA
  CSQ$Gene = NA
  CSQ$Mutation = NA
  head(CSQ)
  
  for(i in 1:nrow(CSQ)) {
    genes = CSQ$Full[i] 
    genes = strsplit(genes, ':')[[1]]
    if(genes == 'unknown') next;
    
    genes.IDs = sapply(strsplit(genes, ' '), function(w) w[1])
    mut.IDs = sapply(strsplit(genes, ' '), function(w) w[2])
    
    
    where = NULL
    for(g in genes.IDs){
      CGC.sbs = CGC[CGC$Gene.Symbol == g, ]
      where = c(where, nrow(CGC.sbs) > 0)
      # COSMIC.sbs = COSMIC[COSMIC$Gene.name == g, ]
      # COSMIC.sbs = COSMIC.sbs[COSMIC.sbs$Mutation.genome.position != "", ]
      # COSMIC.sbs = COSMIC.sbs[duplicated(COSMIC.sbs$Mutation.genome.position) , ]
      
      # sapply(coords, function(w){
      #   pos = unique(COSMIC.sbs$Mutation.genome.position)
      #   hit = NULL
      #   
      #   tomatch = as.integer(strsplit(w, '_')[[1]][2])
      # 
      #   for(p in pos) {
      #     left = strsplit(p, '-')[[1]][1]
      #     left = strsplit(left, ':')[[1]][2]
      #     right = strsplit(p, '-')[[1]][2]
      #     left = as.integer(left)
      #     right = as.integer(right)
      #     hit = c(hit, left <= tomatch & tomatch <= right)
      #   }
      #   hit
      # })
      
      # COSMIC.sbs COSMIC
    }
    
    
    
    if(!is.null(where) && any(where))
    {
      CSQ$CGC[i] = unique(genes.IDs[where])
      CSQ$Gene[i] = paste(genes.IDs[where], collapse = ':')
      CSQ$Mutation[i] = paste(mut.IDs[where], collapse = ':')
    }
    # COSMIC.sbs = COSMIC[COSMIC$Gene.name == df.COSMIC[i, 'Hugo_Symbol'], ]
    # head(COSMIC.sbs)
  }
  head(CSQ)
  nrow(CSQ) == nrow(WES.NV)
  rownames(CSQ) = rownames(WES.NV)
  
  CSQ[is.na(CSQ)] = ''
  
  save(CSQ, file = paste('ANNOTATED-CGC-', patient, '.RData', sep = ''))
}
