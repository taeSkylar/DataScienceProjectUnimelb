library(bigsnpr)
setwd('C:/Users/olladapu/Downloads/plink_win64_20220402')

snp_readBed('filtered_H2_dataset_final_1000G_filtered.bed')

read_snp <- snp_attach('filtered_H2_dataset_final_1000G_filtered.rds')


G   <- read_snp$genotypes
CHR <- read_snp$map$chromosome
POS <- read_snp$map$physical.pos
y   <- read_snp$fam$affection - 1
NCORES <- nb_cores()


geno_stats <- fread('filtered_H2_dataset_final_1000G_filtered.assoc')

set.seed(1)
ind.train <- sample(nrow(G),3000)
ind.test <- setdiff(rows_along(G),ind.train)


names(geno_stats)<- c("chr","rsid","pos","a0","F_A","F_U","a1","CHISQ","p","beta")

map <- read_snp$map[,-(2:3)]

names(map) <- c("chr","pos","a0","a1")


info_snp <- snp_match(geno_stats,map)

beta <- rep(NA,ncol(G))

beta[info_snp$`_NUM_ID_`] <- log10(info_snp$beta)

lpval <- rep(NA,ncol(G))

lpval[info_snp$`_NUM_ID_`] <- -log10(info_snp$p)


big_apply(G,a.FUN=function(X,ind){
  X[,ind]<- replace(X[,ind],is.na(X[,ind]),0)},a.combine = 'plus')


all_keep <- snp_grid_clumping(G,CHR,POS,ind.row=ind.train,
                              lpS = lpval,exclude = which(is.na(lpval)),
                              ncores=NCORES)



multi_PRS <- snp_grid_PRS(G,all_keep,beta,lpval,ind.row = ind.train,
                          n_thr_lpS = 50,ncores = NCORES)

dim(multi_PRS)



