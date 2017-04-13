library(dplyr)

# genes - vector of genes to be annotated
gene_info = function(genes){
  # Expression information
  load('~/Documents/git/david_gitlab/Human_eyeIntegration_App/www/mean_rank_decile.Rdata')
  
  eye.indices = grep("(Cornea)|(Retina)|(RPE)", mean_rank_decile$Sub_Tissue)
  
  expr.filt = mean_rank_decile[eye.indices, ] %>%
    filter(Gene.Name %in% genes) %>% 
    mutate(Tissue = `Sub_Tissue`) %>%
    ungroup() %>%
    dplyr::select(Gene.Name, Tissue, Rank, Decile) %>%
    arrange(Gene.Name, Tissue)
  
  expr.spread = expr.filt %>% dplyr::select(-Rank) %>% spread(Tissue, Decile)
  
  # Adding network information
  load("~/Documents/git/david_gitlab/Human_eyeIntegration_paper/data/retina_connectivity.Rdata")
  retina.connect = data.frame(id = character(0), kWithin = double(0))
  for(i in 1:length(retina.mod.connect.list)){
    retina.connect = bind_rows(retina.connect, retina.mod.connect.list[[i]] %>% dplyr::select(id, kWithin))
  }
  retina.connect = retina.connect %>% mutate(Gene.Name = id) %>% dplyr::select(Gene.Name, kWithin)
  
  load("~/Documents/git/david_gitlab/Human_eyeIntegration_paper/data/rpe_connectivity.Rdata")
  rpe.connect = data.frame(id = character(0), kWithin = double(0))
  for(i in 1:length(rpe.mod.connect.list)){
    rpe.connect = bind_rows(rpe.connect, rpe.mod.connect.list[[i]] %>% dplyr::select(id, kWithin))
  }
  rpe.connect = rpe.connect %>% mutate(Gene.Name = id) %>% dplyr::select(Gene.Name, kWithin)
  
  expr.net = left_join(expr.spread, retina.connect, by = "Gene.Name") %>% mutate(Retina.kWithin = kWithin) %>% dplyr::select(-kWithin) %>% 
    left_join(rpe.connect, by = "Gene.Name") %>% mutate(RPE.kWithin = kWithin) %>% dplyr::select(-kWithin)
  
  # Adding closest genes
  # Retina
  load("~/Documents/git/david_gitlab/Human_eyeIntegration_paper/data/retina_lengthScaledTPM_network_final/network_structure/retina_network_structure.Rdata")
  edges.retina = edges.retina %>% dplyr::select(from, to, length) %>% arrange(desc(length))
  
  k.edges = data.frame(Gene.Name = character(0), Retina.Gene.1st = character(0), Retina.Gene.2nd = character(0), 
                       Retina.Gene.3rd = character(0), Retina.Gene.4th = character(0), Retina.Gene.5th = character(0))
  for(gene in expr.net$Gene.Name){
    edges.retina.filt = filter(edges.retina, from == gene | to == gene)[1:k, ]
    k.genes = paste(edges.retina.filt$from, edges.retina.filt$to, sep = "=")
    k.genes = gsub(paste0("(^", gene, "=)|(=", gene, "$)"), "", k.genes)
    
    k.genes = gsub("NA=NA", NA, k.genes)
    
    k.edges = bind_rows(k.edges, data.frame(Gene.Name = gene, Retina.Gene.1st = k.genes[1], Retina.Gene.2nd = k.genes[2], 
                                            Retina.Gene.3rd = k.genes[3], Retina.Gene.4th = k.genes[4], Retina.Gene.5th = k.genes[5]))
  }
  
  expr.net.ret = left_join(expr.net, k.edges, by = "Gene.Name")
  
  # RPE
  load("~/Documents/git/david_gitlab/Human_eyeIntegration_paper/data/rpe_lengthScaledTPM_network_final/network_structure/rpe_network_structure.Rdata")
  edges.rpe = edges.rpe %>% dplyr::select(from, to, length) %>% arrange(desc(length))
  
  k.edges = data.frame(Gene.Name = character(0), RPE.Gene.1st = character(0), RPE.Gene.2nd = character(0), 
                       RPE.Gene.3rd = character(0), RPE.Gene.4th = character(0), RPE.Gene.5th = character(0))
  for(gene in expr.net$Gene.Name){
    edges.rpe.filt = filter(edges.rpe, from == gene | to == gene)[1:k, ]
    k.genes = paste(edges.rpe.filt$from, edges.rpe.filt$to, sep = "=")
    k.genes = gsub(paste0("(^", gene, "=)|(=", gene, "$)"), "", k.genes)
    
    k.genes = gsub("NA=NA", NA, k.genes)
    
    k.edges = bind_rows(k.edges, data.frame(Gene.Name = gene, RPE.Gene.1st = k.genes[1], RPE.Gene.2nd = k.genes[2], 
                                            RPE.Gene.3rd = k.genes[3], RPE.Gene.4th = k.genes[4], RPE.Gene.5th = k.genes[5]))
  }
  
  expr.net.ret.rpe = left_join(expr.net.ret, k.edges, by = "Gene.Name")
  
  expr.net.ret.rpe$Retina.kWithin[is.na(expr.net.ret.rpe$Retina.kWithin)] = 0
  expr.net.ret.rpe$RPE.kWithin[is.na(expr.net.ret.rpe$RPE.kWithin)] = 0
  
  expr.net.ret.rpe = expr.net.ret.rpe %>% mutate(Sum = `Cornea - Adult Tissue` + `Cornea - Cell Line` + `Cornea - Fetal Tissue` + `Retina - Adult Tissue` + 
                                                   `Retina - Stem Cell Line` + `RPE - Adult Tissue` + `RPE - Cell Line` + `RPE - Fetal Tissue` + 
                                                   `RPE - Stem Cell Line` + Retina.kWithin + RPE.kWithin) %>% 
    select(c(Gene.Name:RPE.kWithin, Sum, Retina.Gene.1st:RPE.Gene.5th))
  
  return(expr.net.ret.rpe)
}
