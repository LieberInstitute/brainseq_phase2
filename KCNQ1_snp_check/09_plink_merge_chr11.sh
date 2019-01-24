#

############################
## make plink file sets, impute without filter to get rs34097980

mkdir -p imputedplink

# Since the PLINK 1 binary format cannot represent genotype probabilities, calls with uncertainty greater than 0.1 are normally treated as missing, and the rest are treated as hard calls.

### make plink files from imputed chunks
for chr in $(seq 11 12); do
	qsub -V -l bluejay,mf=30G,h_vmem=35G,h_stack=256M,h_fsize=100G -cwd -b y /users/ajaffe/bin/plink --bgen /dcs01/ajaffe/Imputation/Merged/GEN_files/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Onmi2pt5_Macrogen_chr$chr.bgen --sample /dcs01/ajaffe/Imputation/Merged/GEN_files/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Onmi2pt5_Macrogen_chr$chr.sample --make-bed --out imputedplink/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_chr$chr.imputed ;
done

# ### fix chrX chr field
# awk '$1==0 {$1=23} {print}' imputed_plinkFiles/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_chr23.imputed.bim > imputed_plinkFiles/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_chr23.imputed.new.bim
# mv imputed_plinkFiles/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_chr23.imputed.new.bim imputed_plinkFiles/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_chr23.imputed.bim

# ### merge autosomal imputed
# /users/ajaffe/bin/plink --bfile imputed_plinkFiles/LIBD_merged_h650_1M_Omni5M_Onmi2pt5_Macrogen_chr1.imputed --merge-list binary_to_merge_withMacrogen.txt --make-bed --out LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2 --memory 250000

# ### loose filter for most datasets
# plink --bfile LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2 --maf 0.005 --geno 0.1 --hwe 1e-10 --make-bed --out LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_Macrogen_imputed_run2_maf005_hwe10_geno10 

# ### get MDS components 
# plink --bfile LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_imputed_run2_maf005_hwe10_geno10 --indep 100 10 1.25 --maf 0.05 --geno 0.1 --hwe 0.000001 --out LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_imputed_run2
# plink --bfile LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_imputed_run2_maf005_hwe10_geno10 --extract LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_imputed_run2.prune.in --cluster --mds-plot 10 --out LIBD_Brain_Illumina_h650_1M_Omni5M_Omni2pt5_imputed_run2

