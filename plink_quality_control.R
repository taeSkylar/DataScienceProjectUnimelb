
#perform Quality control on the genotypes using PLink

#Quality control criteria

#1.Missingness per SN: --geno 0.1
#2.Missingness per individual: --mind 0.1
#3.Minor allele frequency: --maf .05
#4.Hardy-Weinberg threshold: 0.0000001 --hwe

system("./plink --bfile ../ADAPTmap_genotypeTOP_20160222_full/ADAPTmap_genotypeTOP_20160222_full \\
       --cow --nonfounders \\
       --geno 0.1 \\
       --mind 0.1 \\
       --maf 0.05 \\
       --hwe 0.0000001 \\
       --allow-no-sex --recode \\
       --out AfterQC")
