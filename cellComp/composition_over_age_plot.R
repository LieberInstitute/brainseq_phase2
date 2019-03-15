## Original: https://github.com/LieberInstitute/DNAm_Hippo/blob/master/brainseq_phase2_composition_over_age_plot_01.R

## script for brainseq phase2 methylation results
library(ggplot2)
library(ggrepel)
library(jaffelab)
library('RColorBrewer')
library('sessioninfo')

theme_set(theme_bw(base_size=30) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
##
#load('/dcl01/lieber/ajaffe/Steve/Hippo_meQTL/rdas/cleanSamples_n694_Mset_SQN_postfiltered.rda')
load('/dcl01/lieber/ajaffe/Steve/Hippo_meQTL/rdas/cleanSamples_n694_processed_data_postfiltered.rda')
pd$ageGroup = cut(pd$Age, breaks = c(-0.5,0,0.6,10,20,50,100))
levels(pd$ageGroup) = c("Fetal","Infant","Child","Teens","Adult","50+")
pd$ageStage = ifelse(pd$Age<0, "Fetal", ifelse(pd$Age>17, "Adult", NA ))
pd$Dx <- factor(pd$Dx, levels=c("Control","Schizo") )
pd[pd$Brain.Region=="DLPFC",'Sample_Plate']=jaffelab::ss(pd[pd$Brain.Region=="DLPFC",'BasePath'],"\\/",9)

## Save for use in limma_dev_adjNeunProp.R
save(pd, file = 'rda/methprop_pd.Rdata')

##### Cell Type over Development Plot #########
dat <- pd
dat$Brain.Region <- toupper(dat$Brain.Region)
dat$Glial = dat$NeuN_neg
dat$NeuN <- dat$NeuN_pos
dat$DA_Neuron <- dat$DA_NEURON
dat = dat[,c('BrNum','Brain.Region','ageGroup', 'Glial', 'NeuN', 'ES', 'NPC', 'DA_Neuron')]
dat = tidyr::gather(dat, key="CellType", value="Proportion", Glial, NeuN, ES, NPC, DA_Neuron)
dat$CellType = factor(dat$CellType,levels=c('NPC', 'ES', 'Glial', 'NeuN', 'DA_Neuron') )
dat$Region=dat$`Brain.Region`

df3 = data.frame(ageGroup = as.numeric(c(NA, NA)),
                 Region = c("DLPFC", "HIPPO"),
                 Proportion = as.numeric(c(NA, NA)))
				 
cell_type = ggplot(data=dat, aes(x=ageGroup,y=Proportion,fill=Region ))  + 
geom_boxplot(col='black',show.legend=FALSE) + 
facet_wrap(~`CellType`,scales='free',nrow=1)  + 
labs(x='Age Group', y = "Cell Type Proportion",fill='Region') +   
geom_point(data = df3, aes(x = ageGroup, y = Proportion, col = Region), size=8, shape=15) +
  scale_fill_manual(values = c(DLPFC = 'darkgoldenrod2', HIPPO = 'steelblue1'), breaks = c("DLPFC","HIPPO")) +
  scale_color_manual(values = c(DLPFC = 'darkgoldenrod2', HIPPO = 'steelblue1'), breaks = c("DLPFC","HIPPO")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = c(0.1, 0.8),legend.background = element_rect(fill = "white", colour = 'black',linetype='solid'), legend.key = element_blank(),axis.title.x=element_blank()) 

ggsave(cell_type, filename='pdf/bothRegions_estimated_cellType_proportions_over_lifespan_n347.pdf',height=8,width=20)
ggsave(cell_type, filename='pdf/bothRegions_estimated_cellType_proportions_over_lifespan_n347.png',height=8,width=20)

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.5.1 Patched (2018-10-29 r75535)
#  os       Red Hat Enterprise Linux Server release 6.9 (Santiago)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-03-11
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package          * version   date       lib source
#  assertthat         0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  bindr              0.1.1     2018-03-13 [1] CRAN (R 3.5.0)
#  bindrcpp           0.2.2     2018-03-29 [1] CRAN (R 3.5.0)
#  BiocGenerics       0.28.0    2018-10-30 [1] Bioconductor
#  bitops             1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  cli                1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout         * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace         1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon             1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  digest             0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr              0.7.8     2018-11-10 [1] CRAN (R 3.5.1)
#  GenomeInfoDb       1.18.1    2018-11-12 [1] Bioconductor
#  GenomeInfoDbData   1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges      1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2          * 3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  ggrepel          * 0.8.0     2018-05-09 [1] CRAN (R 3.5.0)
#  glue               1.3.0     2018-07-17 [1] CRAN (R 3.5.1)
#  gtable             0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools          0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets        1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv             1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges            2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab         * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  labeling           0.3       2014-08-23 [2] CRAN (R 3.5.0)
#  later              0.7.5     2018-09-18 [2] CRAN (R 3.5.1)
#  lattice            0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval           0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma              3.38.3    2018-12-02 [1] Bioconductor
#  magrittr           1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  munsell            0.5.0     2018-06-12 [2] CRAN (R 3.5.0)
#  pillar             1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig          2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr               1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises           1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr              0.2.5     2018-05-29 [2] CRAN (R 3.5.0)
#  R6                 2.3.0     2018-10-04 [2] CRAN (R 3.5.1)
#  rafalib          * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer     * 1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp               1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl              1.95-4.11 2018-07-15 [2] CRAN (R 3.5.1)
#  rlang              0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote            * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors          0.20.1    2018-11-09 [1] Bioconductor
#  scales             1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented          0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr              0.11      2018-10-23 [1] CRAN (R 3.5.1)
#  sessioninfo      * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  tibble             2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyr              0.8.2     2018-10-28 [2] CRAN (R 3.5.1)
#  tidyselect         0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr              2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun               0.4       2018-10-23 [1] CRAN (R 3.5.1)
#  XVector            0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc           1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
