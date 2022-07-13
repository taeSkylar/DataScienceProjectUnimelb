

prs_calculation <- function(bed,debug,n,weights){
  
  bed_data <- snp_attach(paste0(bed,".rds"))
  
  bed2 <- bed_data$genotypes
  
  bed2 <- big_apply(bed2,a.FUN = function(X,ind){
    X[,ind] <- X[,ind]*weights[ind]
  },a.combine = 'plus')
  
  bed2 <- as_FBM(bed2)
  
  prs_scores <- big_apply(bed2,a.FUN=function(X,ind){
    rowSums(X[ind,])
  })

  return(prs_scores)
  
}




cal_prs <- function(bfile,assoc,p,pheno){
  
  assoc <- fread(assoc,h=T)
  bim <- fread(paste0(bfile,".bim"),h=F) %>% as.data.frame %>% select(V2)
  
  
  fam <- as.data.frame(fread(paste0(bfile,".fam"),h=F))
  
  
  beta_name <- "OR"
  p_name <- "P"
  snp_name <- "SNP"
  
  
  
  assoc <- assoc %>% select("SNP","OR","P")
  
  colnames(assoc) <- c("SNP","BETA","P")
  colnames(bim) <- c("SNP")
  ordered <- merge(bim,assoc,all.x=T,all.y=F,sort=F)
  
  weights <- as.matrix(ordered$BETA)
    
  bed <- paste0(bfile,".bed")
  
  nsnp <- as.integer(nrow(weights))
  
  n=as.integer(nrow(fam))
  
  s <- prs_calculation(bed,F,n,weights)
  
  
}

snp_readBed('h6.rds')
