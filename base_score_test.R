#load libraries
library(tidyverse)
library(data.table)
library(qqman)
library(bigsnpr)

#set working directory
setwd('C:/Users/olladapu/Downloads/plink_win64_20220402')

bed <- 'H2_dataset_final_non_1000G_filtered'

#create .rds format from bed file
snp_readBed(paste0(bed,".bed"))

#read the .rds file to read genotype data
bed_data <- snp_attach(paste0(bed,".rds"))

#read genotype data
bed2 <- bed_data$genotypes

#DO NOT RUN it
#head(bed2)

#import weights from the Mavaddat paper

snp_annotation <- read_csv('SNP_annotation.csv')

#create a SNP and DBSNP mapping dictionary
snp_map_dict={}
snp_map_dict_rev={}
snp_map_dict[snp_annotation$RS_dbSNP137]=snp_annotation$rsID
snp_map_dict_rev[snp_annotation$rsID]=snp_annotation$RS_dbSNP137

#read weights from mavaddata file
mgwas2 <- read.table('sample.txt',sep = " ")

#create a mavaddat mapping dictionary map
weights_dict <- {}
weights_dict[mgwas2$ID]<- mgwas2$effect_weight

#get all the SNPs from the file
cols12 <- bed_data$map$marker.ID

weights=c()
shortlisted_snps=c()

shortlisted_snps_2=c()

for(i in mgwas2$ID){
  if(!is.na(snp_map_dict_rev[i])){
    shortlisted_snps_2 <- c(shortlisted_snps_2,snp_map_dict_rev[i]) }
}


#write shortlisted SNPs to a file 
write.table(shortlisted_snps,'short_snp_list.txt',sep='\n')


#read the data from the shortlisted rds files
#RUN this below command in a terminal to generate the h_new3 dataset with the shortlisted SNps

#NOTE: please note that after every run you need to generate a fresh h_new bed file and 

#run line 69 and 71,73 to update the dataset

#./plink --bfile H2_dataset_final_non_1000G_filtered --extract short_snp_list.txt  --make-bed --out h_new3

snp_readBed('h_new3.bed')

bed_data <- snp_attach('h_new3.rds')

bed2 <- bed_data$genotypes

cols_list <- bed_data$map$marker.ID


#import the weights for the shortlisted snps
for(i in cols_list){
  if(!is.na(snp_map_dict_rev[i])){
    weights <- c(weights,weights_dict[snp_map_dict_rev[i]]) }
}


bed2 <- big_apply(bed2,a.FUN = function(X,ind){
  X[,ind] <- X[,ind]*weights[ind]
},a.combine = 'plus')

bed2 <- as_FBM(bed2)

prs_scores <-big_apply(bed2,a.FUN = function(X,ind){rowSums(X[ind,])},ind = rows_along(bed2),a.combine = 'plus')

prs_scores_2 <- replace(prs_scores,is.na(prs_scores),0)

AUC(prs_scores_2,bed_data$fam$affection-1)
