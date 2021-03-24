# module load plink

library("readxl")
library("here")
library("sessioninfo")

subset_snps <- read_xlsx(
    here(
        "eQTL_full",
        "subset_hancock",
        "LIBDeQTLs_OPRM1.xlsx"
    )
)$snp

snps_file <- here(
    "eQTL_full",
    "subset_hancock",
    "snps_to_extract.txt"
)

writeLines(subset_snps, snps_file)

bfile <- here(
    "genotype_data",
    "BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6"
)
newbfile <- here(
    "eQTL_full",
    "subset_hancock",
    "BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6_hancock"
)
system(paste("plink --bfile", bfile, "--extract", snps_file, " --make-bed --out", newbfile))
# PLINK v1.90b6.6 64-bit (10 Oct 2018)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to /dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/subset_hancock/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6_hancock.log.
# Options in effect:
#   --bfile /dcl01/lieber/ajaffe/lab/brainseq_phase2/genotype_data/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6
#   --extract /dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/subset_hancock/snps_to_extract.txt
#   --make-bed
#   --out /dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/subset_hancock/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6_hancock
#
# 257850 MB RAM detected; reserving 128925 MB for main workspace.
# Allocated 96693 MB successfully, after larger attempt(s) failed.
# 7023860 variants loaded from .bim file.
# 551 people (365 males, 186 females) loaded from .fam.
# --extract: 33 variants remaining.
# Warning: At least 18 duplicate IDs in --extract file.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 551 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.997305.
# 33 variants and 551 people pass filters and QC.
# Note: No phenotypes present.
# --make-bed to
# /dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/subset_hancock/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6_hancock.bed
# +
# /dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/subset_hancock/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6_hancock.bim
# +
# /dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/subset_hancock/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6_hancock.fam
# ... done.


system(paste("plink --bfile", newbfile, "--recode vcf --out", paste0(newbfile, "_vcf")))
# PLINK v1.90b6.6 64-bit (10 Oct 2018)           www.cog-genomics.org/plink/1.9/
# (C) 2005-2018 Shaun Purcell, Christopher Chang   GNU General Public License v3
# Logging to /dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/subset_hancock/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6_hancock_vcf.log.
# Options in effect:
#   --bfile /dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/subset_hancock/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6_hancock
#   --out /dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/subset_hancock/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6_hancock_vcf
#   --recode vcf
#
# 257850 MB RAM detected; reserving 128925 MB for main workspace.
# Allocated 96693 MB successfully, after larger attempt(s) failed.
# 33 variants loaded from .bim file.
# 551 people (365 males, 186 females) loaded from .fam.
# Using 1 thread (no multithreaded calculations invoked).
# Before main variant filters, 551 founders and 0 nonfounders present.
# Calculating allele frequencies... done.
# Total genotyping rate is 0.997305.
# 33 variants and 551 people pass filters and QC.
# Note: No phenotypes present.
# --recode vcf to
# /dcl01/lieber/ajaffe/lab/brainseq_phase2/eQTL_full/subset_hancock/BrainSeq_Phase2_RiboZero_Genotypes_n551_maf05_geno10_hwe1e6_hancock_vcf.vcf
# ... done.
