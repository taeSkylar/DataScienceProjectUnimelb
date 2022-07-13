library(bigsnpr)

setwd('C:/Users/olladapu/Downloads/plink_win64_20220402')

snp_readBed('h14.bed')

snp_data <- snp_attach('h14.rds')

G <- snp_data$genotypes

CHR <- snp_data$map$chromosome

POS <- snp_data$map$physical.pos

y <- snp_data$fam$affection-1

pop <- snp_data$fam$family.ID

NCORES <- nb_cores()

set.seed(1)

ind.train <- sample(nrow(G),5000)

ind.test <- setdiff(rows_along(G),ind.train)


svd <- big_randomSVD(G,big_scale(),ncores = NCORES)

plot(svd)