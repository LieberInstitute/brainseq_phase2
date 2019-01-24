#!/bin/bash

wget http://walters.psycm.cf.ac.uk/clozuk_pgc2.meta.sumstats.txt.gz
wget http://walters.psycm.cf.ac.uk/clozuk_pgc2.meta.sumstats.info9.snplist.txt.gz

gunzip clozuk_pgc2.meta.sumstats.txt.gz
gunzip clozuk_pgc2.meta.sumstats.info9.snplist.txt.gz


git clone https://github.com/bulik/ldsc.git /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/ldsc
cd /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/ldsc
conda env create --file environment.yml
source activate ldsc
cd /dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/psycm

/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/ldsc/munge_sumstats.py -h

/dcl01/lieber/ajaffe/lab/brainseq_phase2/twas/ldsc/munge_sumstats.py --sumstats clozuk_pgc2.meta.sumstats.txt --N-cas 40675 --N-con 64643 --out clozuk_pgc2.meta.reformatted --frq "Freq.A1" --maf-min -0.1
# *********************************************************************
# * LD Score Regression (LDSC)
# * Version 1.0.0
# * (C) 2014-2015 Brendan Bulik-Sullivan and Hilary Finucane
# * Broad Institute of MIT and Harvard / MIT Department of Mathematics
# * GNU General Public License v3
# *********************************************************************
# Call:
# ./munge_sumstats.py \
# --out clozuk_pgc2.meta.reformatted \
# --frq Freq.A1 \
# --N-con 64643.0 \
# --maf-min -0.1 \
# --N-cas 40675.0 \
# --sumstats clozuk_pgc2.meta.sumstats.txt
#
# Interpreting column names as follows:
# Freq.A1:    Allele frequency
# A1:    Allele 1, interpreted as ref allele for signed sumstat.
# P:    p-Value
# A2:    Allele 2, interpreted as non-ref allele for signed sumstat.
# SNP:    Variant ID (e.g., rs number)
# OR:    Odds ratio (1 --> no effect; above 1 --> A1 is risk increasing)
#
# Reading sumstats from clozuk_pgc2.meta.sumstats.txt into memory 5000000 SNPs at a time.
# .. done
# Read 8171061 SNPs from --sumstats file.
# Removed 0 SNPs with missing values.
# Removed 0 SNPs with INFO <= 0.9.
# Removed 0 SNPs with MAF <= -0.1.
# Removed 0 SNPs with out-of-bounds p-values.
# Removed 1670587 variants that were not SNPs or were strand-ambiguous.
# 6500474 SNPs remain.
# Removed 0 SNPs with duplicated rs numbers (6500474 SNPs remain).
# Median value of OR was 0.99957, which seems sensible.
# Writing summary statistics for 6500474 SNPs (6500474 with nonmissing beta) to clozuk_pgc2.meta.reformatted.sumstats.gz.
#
# Metadata:
# Mean chi^2 = 1.867
# Lambda GC = 1.582
# Max chi^2 = 194.697
# 15178 Genome-wide significant SNPs (some may have been removed by filtering).
#
# Conversion finished at Thu Jan 24 10:58:27 2019
# Total time elapsed: 3.0m:7.52s

## Uncompress for using 'head'
gunzip clozuk_pgc2.meta.reformatted.sumstats.gz

## Check that it indeed removed non-SNPs
R

x <- readr::read_tsv('clozuk_pgc2.meta.reformatted.sumstats')
table(nchar(x$A2))
table(nchar(x$A1))
addmargins(table('A2' = nchar(x$A2) == 1, 'A1' = nchar(x$A1) == 1))
#       A1
# A2        TRUE     Sum
#   TRUE 6500474 6500474
#   Sum  6500474 6500474

y <- readr::read_tsv('clozuk_pgc2.meta.sumstats.txt')
table(nchar(y$A2))
table(nchar(y$A1))
# > table(nchar(y$A2))
#
#       1       2       3       4       5       6       7       8       9      10
# 7807339  202949   45833   33565   39590   13987    6904    3627    3839    2477
#      11      12      13      14      15      16      17      18      19      20
#    2359    1644    1701    1005     941     685     676     428     373     282
#      21      22      23      24      25      26      27      28      29      30
#     227     150     143      73      83      44      41      18      21      12
#      31      32      33      34      35      36      37      38      39      40
#       9      10       9       4       4       2       2       1       1       2
#      44
#       1
# > table(nchar(y$A1))
#
#       1       2       3       4       5       6       7       8       9      10
# 8042407   62321   30310   10058   14065    4200    1940     984    1099     710
#      11      12      13      14      15      16      17      18      19      20
#     654     444     505     266     255     166     160     119      98      62
#      21      22      23      24      25      26      27      28      29      30
#      66      39      29      23      19      22      11       7       7       1
#      31      32      33      34      36      38      42
#       6       2       1       1       1       2       1

table(nchar(y$A2) == 1)
table(nchar(y$A1) == 1)
addmargins(table('A2' = nchar(y$A2) == 1, 'A1' = nchar(y$A1) == 1))
#        A1
# A2        FALSE    TRUE     Sum
#   FALSE       0  363722  363722
#   TRUE   128654 7678685 7807339
#   Sum    128654 8042407 8171061


## Looks like there are 7678685 - 6500474 = 1178211 strand ambiguous snps?

# > 7678685 - 6500474
# [1] 1178211
# > 8171061 - 6500474
# [1] 1670587
