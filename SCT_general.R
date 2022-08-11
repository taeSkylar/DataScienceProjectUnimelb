library(tidyverse)
library(data.table)
library(qqman)
library(bigsnpr)

setwd('C:/Users/olladapu/Downloads/plink_win64_20220402')

obj.bigSNP<- snp_attach("filtered_H2_dataset_final_1000G_filtered.rds")

head(obj.bigSNP$genotypes)
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection - 1
NCORES <- nb_cores()

sumstats <-fread('filtered_H2_dataset_final_1000G_filtered.assoc')

set.seed(1)
ind.train <- sample(nrow(G),3000)
ind.test <- setdiff(rows_along(G),ind.train)

names(sumstats)<- c("chr","rsid","pos","a0","F_A","F_U","a1","CHISQ","p","beta")

map <- obj.bigSNP$map[,-(2:3)]

names(map) <- c("chr","pos","a0","a1")


info_snp <- snp_match(sumstats,map)

beta <- rep(NA,ncol(G))


beta[info_snp$`_NUM_ID_`] <- log(info_snp$beta)

lpval <- rep(NA,ncol(G))

lpval[info_snp$`_NUM_ID_`] <- -log10(info_snp$p)

# The clumping step might take some time to complete
all_keep <- snp_grid_clumping(G, CHR, POS, ind.row = ind.train,
                              lpS = lpval, exclude = which(is.na(lpval)),
                              ncores = NCORES)
attr(all_keep, "grid")



multi_PRS <- snp_grid_PRS(G,all_keep,beta,lpval,ind.row = ind.train,
                          n_thr_lpS = 50,ncores = NCORES)

dim(multi_PRS)

final_mod <- snp_grid_stacking(multi_PRS, y[ind.train], ncores = NCORES, K = 4)
summary(final_mod$mod)

