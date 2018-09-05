##
library(jaffelab)


hippo = read.csv("results/raggr_179_snps_hippo_eqtls_fdr01.csv")
dlp = read.csv("results/raggr_179_snps_dlp_eqtls_fdr01.csv")


## can you guys also tabulate how many of the SNPs only associate with 1 gene
h = hippo[order(hippo$SNP),]
h_snps = data.frame(SNP = unique(h$SNP), ngenes = NA)
for (i in 1:nrow(h_snps)) {
	s = h_snps$SNP[i]
	h_snps$ngenes[i] = length(unique(h$EnsemblGeneID[which(h$SNP == s)]))
}
table(h_snps$ngenes)
   # 1    2    3    4    5    6    7    8    9   10   11   12	13
# 1982 1353  653  580  317   91  229   28  106   47   58   59    7

d = dlp[order(dlp$SNP),]
d_snps = data.frame(SNP = unique(d$SNP), ngenes = NA)
for (i in 1:nrow(d_snps)) {
	s = d_snps$SNP[i]
	d_snps$ngenes[i] = length(unique(d$EnsemblGeneID[which(d$SNP == s)]))
}
table(d_snps$ngenes)
   # 1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
# 2447 1790  679  686  329  292  107  136   50   49   70   73   36   28    7    1





#### Index-SNP level

## can you guys also tabulate how many of the SNPs only associate with 1 gene
h = hippo[order(hippo$IndexSNP),]
h_snps = data.frame(IndexSNP = unique(h$IndexSNP), ngenes = NA)
for (i in 1:nrow(h_snps)) {
	s = h_snps$IndexSNP[i]
	h_snps$ngenes[i] = length(unique(h$EnsemblGeneID[which(h$IndexSNP == s)]))
}
table(h_snps$ngenes)
 # 1  2  3  4  5  6  7  8  9 10 11 12 14 18 22 25
# 38 18 13  9  4  1  4  2  5  1  1  1  3  1  1  1

d = dlp[order(dlp$IndexSNP),]
d_snps = data.frame(IndexSNP = unique(d$IndexSNP), ngenes = NA)
for (i in 1:nrow(d_snps)) {
	s = d_snps$IndexSNP[i]
	d_snps$ngenes[i] = length(unique(d$EnsemblGeneID[which(d$IndexSNP == s)]))
}
table(d_snps$ngenes)
 # 1  2  3  4  5  7  8  9 10 11 12 13 14 15 17 20 23 24
# 37 19 20  8  7  2  3  4  4  2  3  1  1  1  1  1  1  1






library(RColorBrewer)
library(VennDiagram)

pal = c('skyblue3','dark orange')
venn.diagram(list(Hippo = unique(hippo$IndexSNP), 
				  DLPFC = unique(dlp$IndexSNP) ), 
	fill = pal[1:2], main="", main.pos = c(.5, .05), cat.cex = 1.9, cex=3,
	margin = .1, imagetype="png",  filename = "venn_unique_IndexSNP_UPDATED.png")



	
	
	
	

