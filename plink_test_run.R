#system command to run plink in R scripts
#--bfile binary input file to convert it to .ped and .map format data
#--cow to mention the type of chromosomes we are dealing with (here it is cow)
#--nonfounders don't include column about parent information
#--allow-no-sex don't include column about sex information
#--recode parameter to change the datatype format to .ped and .map format
#--out output file name

system("./plink --bfile ../ADAPTmap_genotypeTOP_20160222_full/ADAPTmap_genotypeTOP_20160222_full \
       --cow --nonfounders \
       --allow-no-sex --recode \
       --out ADAPTmap_TOP")
