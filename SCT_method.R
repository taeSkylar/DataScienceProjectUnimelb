library(tidyverse)
library(data.table)
library(qqman)
library(bigsnpr)

setwd('C:/Users/olladapu/Downloads/plink_win64_20220402')
eigenValues <- read_tsv('H2_dataset.eigenval',col_names = T)
eigenVectors <- read_tsv('H2_dataset.eigenvec',col_names = T)

bc <- read_delim('H2_dataset_final_1000G_filtered.fam',delim=" ",col_names = F)


eigen_percent <- round((eigenValues / (sum(eigenValues))*100), 2)

# PCA plot
ggplot(data = eigenVectors) +
  geom_point(mapping = aes(x = bc$X1, y = bc$X3, color = as.factor(bc$X6),shape=as.factor(bc$X6)), size = 3, show.legend = TRUE ) +
  geom_hline(yintercept = 0, linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  labs(title = "PCA of selected 10 components",
       x = paste0("Principal component 1 (",eigen_percent[3,1]," %)"),
       y = paste0("Principal component 2 (",eigen_percent[4,1]," %)"),
       colour = "breast cancer", shape = "breast cancers") +
  theme_minimal()


mgwas2 <- read.table('sample.txt',sep = " ")

big_snp <- snp_attach('h6.rds')

geno <- big_snp$genotypes

gwas_results <- fread('H2_dataset.assoc')


snp_annotation <- read_csv('SNP_annotation.csv')

snp_map_dict={}

snp_map_dict_rev={}

snp_map_dict[snp_annotation$RS_dbSNP137]=snp_annotation$rsID

snp_map_dict_rev[snp_annotation$rsID]=snp_annotation$RS_dbSNP137

weights_dict <- {}

weights_dict[mgwas2$ID]<- mgwas2$effect_weight

shortlisted_snps=c()

weights=c()

for(i in cols12){
  if(!is.na(snp_map_dict_rev[i])){
  weights <- c(weights,weights_dict[snp_map_dict_rev[i]]) }
}

geno$short_snps_dt$shortlisted_snps



library(dbplyr)

head(geno[,short_snps_dt$shortlisted_snps])

snp_subset(big_snp,ind.col=short_snps_dt$shortlisted_snps)

write.table(list(short_snps_dt$shortlisted_snps)[[1]],'snps.txt',row.names = F)

snp_readBed('h13.bed')

short_snp <- snp_attach('h5.rds')

short_snp_geno <- short_snp$genotypes


short_snp_geno


results <- read.table('plink.profile',header = T)

mean(results$SCORESUM)
