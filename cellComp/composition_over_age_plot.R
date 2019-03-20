## Original: https://github.com/LieberInstitute/DNAm_Hippo/blob/master/brainseq_phase2_composition_over_age_plot_01.R

library('sessioninfo')
library('SummarizedExperiment')
library('ggplot2')
library('jaffelab')
library('limma')

source('../development/load_funs.R')

## Load the data
load("RNA_cell_proportions_brainSeq_phase2.rda")
load("../expr_cutoff/rse_gene.Rdata")

## Re-weight so the sum is 1? Yes
summary(rowSums(propEsts))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1.205   1.548   1.579   1.575   1.617   1.871
propEsts <- sweep(propEsts, 1, rowSums(propEsts), '/')
summary(rowSums(propEsts))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1       1       1       1       1       1
   
## Add age groups
colData(rse_gene)$ageGroup <- factor(with(colData(rse_gene), dplyr::case_when(
        Age < 0 ~ 'Prenatal',
        Age >= 0 & Age < 1 ~ 'Infant',
        Age >= 1 & Age < 10 ~ 'Child',
        Age >= 10 & Age < 20 ~ 'Teen',
        Age >= 20 & Age < 50 ~ 'Adult',
        Age >= 50 ~ '50+'
)), levels = c('Prenatal', 'Infant', 'Child', 'Teen', 'Adult', '50+'))

colData(rse_gene)$ageStage <- factor(with(colData(rse_gene), dplyr::case_when(
    Age < 0 ~ 'Prenatal',
    Age > 17 ~ 'Adult'
)), levels = c('Prenatal', 'Adult'))

## Add cell type info
m_cell <- match(colnames(rse_gene), rownames(propEsts))
colData(rse_gene) <- cbind(colData(rse_gene), propEsts[m_cell, ])

## Add age linear splines
fetal <- ifelse(colData(rse_gene)$Age < 0, 1,0)
birth <- colData(rse_gene)$Age
birth[birth < 0] <- 0 # linear spline
infant <- colData(rse_gene)$Age - 1
infant[infant < 0] <- 0 # linear spline
child <- colData(rse_gene)$Age - 10
child[child < 0] <- 0 # linear spline
teen <- colData(rse_gene)$Age - 20
teen[teen < 0] <- 0 # linear spline
adult <- colData(rse_gene)$Age - 50
adult[adult < 0] <- 0 # linear spline

colData(rse_gene)$fetal <- fetal
colData(rse_gene)$birth <- birth
colData(rse_gene)$infant <- infant
colData(rse_gene)$child <- child
colData(rse_gene)$teen <- teen
colData(rse_gene)$adult <- adult

## Add means
colData(rse_gene)$mean_mitoRate <- mean(colData(rse_gene)$mitoRate)
colData(rse_gene)$mean_totalAssignedGene <- mean(colData(rse_gene)$totalAssignedGene)
colData(rse_gene)$mean_RIN <- mean(colData(rse_gene)$RIN)

## Save for later use (if necessary)
pd <- colData(rse_gene)
save(pd, file = 'methprop_pd.Rdata')


##### Cell Type over Development Plot #########

theme_set(theme_bw(base_size=30) + 
		  theme(panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				 plot.title = element_text(hjust = 0.5),
				 legend.position="none"))
                 
                 
dat <- pd
dat$Region <- toupper(dat$Region)
dat$FPN <- dat$Fetal_replicating
dat$FQN <- dat$Fetal_quiescent
dat = as.data.frame(dat[,c('BrNum','Region','ageGroup', 'FPN', 'FQN', 'OPC', 'Neurons', 'Astrocytes', 'Oligodendrocytes', 'Microglia', 'Endothelial')])
dat = tidyr::gather(dat, key="CellType", value="Proportion", FPN, FQN, OPC, Neurons, Astrocytes, Oligodendrocytes, Microglia, Endothelial)
dat$CellType = factor(dat$CellType,levels=c('FPN', 'FQN', 'Neurons', 'Microglia', 'OPC', 'Astrocytes', 'Oligodendrocytes', 'Endothelial') )

df3 = data.frame(ageGroup = as.numeric(c(NA, NA)),
                 Region = c("DLPFC", "HIPPO"),
                 Proportion = as.numeric(c(NA, NA)))
                 
plot_code <- function(dat, legend_pos = c(0.053, 0.8)) {
    ggplot(data=dat, aes(x=ageGroup,y=Proportion,fill=Region ))  + 
        geom_boxplot(col='black',show.legend=FALSE) + 
        facet_wrap(~`CellType`,scales='free',nrow=1)  + 
        labs(x='Age Group', y = "Cell Type Proportion",fill='Region') +   
        geom_point(data = df3, aes(x = ageGroup, y = Proportion, col = Region), size=8, shape=15) +
        scale_fill_manual(values = c(DLPFC = 'darkgoldenrod2', HIPPO = 'steelblue1'), breaks = c("DLPFC","HIPPO")) +
        scale_color_manual(values = c(DLPFC = 'darkgoldenrod2', HIPPO = 'steelblue1'), breaks = c("DLPFC","HIPPO")) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = legend_pos,legend.background = element_rect(fill = "white", colour = 'black',linetype='solid'), legend.key = element_blank(),axis.title.x=element_blank()) 
}
				 
cell_type <- plot_code(dat)
ggsave(cell_type, filename='bothRegions_estimated_cellType_proportions_over_lifespan.pdf',height=8,width=30)
ggsave(cell_type, filename='bothRegions_estimated_cellType_proportions_over_lifespan.png',height=8,width=30)


cell_type_main <- plot_code(subset(dat, CellType %in% c('FQN', 'Neurons', 'Microglia')), legend_pos = c(0.2, 0.8))
ggsave(cell_type_main, filename='bothRegions_estimated_cellType_proportions_over_lifespan_main.pdf',height=8,width=15)


## Make plots over age
## Based on https://github.com/LieberInstitute/brainseq_phase2/blob/89e48d9e612381d857aef358595933c0b9adbf73/development/explore_limma_dev_top.R#L96-L115
p_cols <- ifelse(colData(rse_gene)$Region == 'HIPPO', 'steelblue1', 'darkgoldenrod2')
l_cols <- c('lightgoldenrod', 'light blue')
age_brks <- c(-1, 0, 1, 10, 20, 50, 100)

design <- get_mods( colData(rse_gene), int = TRUE)$mod

plot_age_mod <-
    design[, c(
        '(Intercept)',
        'Age',
        'RegionHIPPO',
        'fetal',
        'birth',
        'infant',
        'child',
        'teen',
        'adult'
    )]

pdf('proportions_over_age_by_cell_type.pdf', width = 14, useDingbats = FALSE)
for(i in seq_along(propEsts)) {
    set.seed(20190315)
    agePlotter(
        y = pd[, colnames(propEsts)[i]],
        age = pd$Age,
        pointColor = p_cols,
        ageBreaks = age_brks,
        ## If we wanted the abbreviated names:
        # mainText = c('FQN', 'FQN', colnames(propEsts)[-c(1:2)])[i],
        mainText = gsub('_', ' ', colnames(propEsts))[i],
        lineColor = l_cols,
        ylab = 'Estimated Proportion',
        mod = plot_age_mod
    )
    legend('top', c('DLPFC', 'HIPPO'), col = l_cols, lwd = 3, bty = 'n', ncol = 1, cex = 1.5)
}
dev.off()

## Use limma like we did at
## https://github.com/LieberInstitute/brainseq_phase2/blob/master/development/limma_dev.R
brnum <- pd$BrNum
props <- t(as.matrix(pd[, colnames(propEsts)]))
corfit <- duplicateCorrelation(props,
        design[, c('(Intercept)', 'RegionHIPPO')], block=brnum)
        
## Only adjust for the variables plotted
mod_lm <- cbind(plot_age_mod, design[, grep(':', colnames(design))])
colnames(mod_lm)
#  [1] "(Intercept)"        "Age"                "RegionHIPPO"
#  [4] "fetal"              "birth"              "infant"
#  [7] "child"              "teen"               "adult"
# [10] "Age:RegionHIPPO"    "RegionHIPPO:fetal"  "RegionHIPPO:birth"
# [13] "RegionHIPPO:infant" "RegionHIPPO:child"  "RegionHIPPO:teen"
# [16] "RegionHIPPO:adult"

## Now fit with limma
fit <- lmFit(props, design = mod_lm, block=brnum,
        correlation = corfit$consensus.correlation)
fit <- eBayes(fit)
top <- topTable(fit, coef = grep(':', colnames(mod_lm)), n = nrow(props),
    sort.by = 'none')
    
## Add some values for convinience
top$P.Bonf <- p.adjust(top$P.Value, 'bonf')
top$P.Bonf.sig <- top$P.Bonf < 0.05
# top$P.Bonf.Less.1perc <- top$P.Bonf < 0.01
top$FDR.sig <- top$adj.P.Val < 0.05

## Use a better name and print
top_cell_prop <- top
options(width = 300)
top_cell_prop


## Save the results
write.table(top_cell_prop, file = 'top_cell_prop.txt', sep = '\t', quote = FALSE)
save(top_cell_prop, file = 'top_cell_prop.Rdata')

## Manually with lm() although that doesn't take into account the repeated measures
# x <- do.call(rbind, mapply(function(i, cell) {
#     tmp <- cbind(data.frame(y = pd[, colnames(propEsts)[i]]), pd[, c('Region', colnames(plot_age_mod)[-c(1, 3)])])
#     fit_lm <- lm(y ~ Region * Age + Region * fetal + Region * birth + Region * infant + Region * infant + Region * child + Region * teen + Region * adult, data = tmp)
#     fit_lm0 <- lm(y ~ Region + Age + fetal + birth + infant + infant + child + teen + adult, data = tmp)
#     res <- as.data.frame(anova(fit_lm, fit_lm0)[2, 5:6])
#     res$cell <- cell
#     rownames(res) <- NULL
#     return(res)
# }, seq_len(ncol(propEsts)), colnames(propEsts), SIMPLIFY = FALSE))
# x$pbonf <- p.adjust(x[, 2], 'bonf')
# x$pbonf < 0.05




top_by_age_group <- function(agegroup) {
    index <- pd$ageGroup == agegroup
    pd <- pd[index, ]
    props <- props[, index]
    brnum <- pd$BrNum
    mod <- with(pd, model.matrix( ~ Region * Age))
    corfit <- duplicateCorrelation(props,
            mod[, c('(Intercept)', 'RegionHIPPO')], block=brnum)
            
    fit <- lmFit(props, design = mod, block=brnum,
            correlation = corfit$consensus.correlation)
    fit <- eBayes(fit)
    top <- topTable(fit, coef = grep(':', colnames(mod)), n = nrow(props),
        sort.by = 'none')
    top$ageGroup <- agegroup
    top$cell <- rownames(top)
    
    ## Compute these globally
    ## Leads to basically the same results (like the same one is P.Bonf < 5%)
    # top$P.Bonf <- p.adjust(top$P.Value, 'bonf')
    # top$P.Bonf.sig <- top$P.Bonf < 0.05
    # top$FDR.sig <- top$adj.P.Val < 0.05
    return(top)        
}


top_cell_by_agegroup <- do.call(rbind, lapply(levels(pd$ageGroup), top_by_age_group))
top_cell_by_agegroup$P.Bonf <- p.adjust(top_cell_by_agegroup$P.Value, 'bonf')
top_cell_by_agegroup$P.Bonf.sig <- top_cell_by_agegroup$P.Bonf < 0.05
top_cell_by_agegroup$FDR.sig <- top_cell_by_agegroup$adj.P.Val < 0.05
rownames(top_cell_by_agegroup) <- NULL
options(width = 150)
subset(top_cell_by_agegroup, P.Bonf.sig)
#           logFC   AveExpr         t      P.Value   adj.P.Val         B ageGroup       cell      P.Bonf P.Bonf.sig FDR.sig
# 37 -0.001422805 0.2399906 -3.755022 0.0002008077 0.001606461 -5.820941    Adult Astrocytes 0.009638768       TRUE    TRUE

## Save the results
write.table(top_cell_by_agegroup, file = 'top_cell_by_agegroup.txt', sep = '\t', quote = FALSE, row.names = FALSE)
save(top_cell_by_agegroup, file = 'top_cell_by_agegroup.Rdata')


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
#  date     2019-03-20
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.0     2017-04-11 [2] CRAN (R 3.5.0)
#  Biobase              * 2.42.0    2018-10-30 [2] Bioconductor
#  BiocGenerics         * 0.28.0    2018-10-30 [1] Bioconductor
#  BiocParallel         * 1.16.6    2019-02-10 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.5.0)
#  cli                    1.0.1     2018-09-25 [1] CRAN (R 3.5.1)
#  colorout             * 1.2-0     2018-05-02 [1] Github (jalvesaq/colorout@c42088d)
#  colorspace             1.4-0     2019-01-13 [2] CRAN (R 3.5.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.5.0)
#  DelayedArray         * 0.8.0     2018-10-30 [2] Bioconductor
#  digest                 0.6.18    2018-10-10 [1] CRAN (R 3.5.1)
#  dplyr                  0.8.0.1   2019-02-15 [1] CRAN (R 3.5.1)
#  GenomeInfoDb         * 1.18.2    2019-02-12 [1] Bioconductor
#  GenomeInfoDbData       1.2.0     2018-11-02 [2] Bioconductor
#  GenomicRanges        * 1.34.0    2018-10-30 [1] Bioconductor
#  ggplot2              * 3.1.0     2018-10-25 [1] CRAN (R 3.5.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.5.1)
#  gtable                 0.2.0     2016-02-26 [2] CRAN (R 3.5.0)
#  htmltools              0.3.6     2017-04-28 [2] CRAN (R 3.5.0)
#  htmlwidgets            1.3       2018-09-30 [1] CRAN (R 3.5.1)
#  httpuv                 1.4.5.1   2018-12-18 [2] CRAN (R 3.5.1)
#  IRanges              * 2.16.0    2018-10-30 [1] Bioconductor
#  jaffelab             * 0.99.21   2018-05-03 [1] Github (LieberInstitute/jaffelab@7ed0ab7)
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.5.1)
#  later                  0.8.0     2019-02-11 [2] CRAN (R 3.5.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.5.1)
#  lazyeval               0.2.1     2017-10-29 [2] CRAN (R 3.5.0)
#  limma                * 3.38.3    2018-12-02 [1] Bioconductor
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.5.0)
#  Matrix                 1.2-15    2018-11-01 [3] CRAN (R 3.5.1)
#  matrixStats          * 0.54.0    2018-07-23 [1] CRAN (R 3.5.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.5.1)
#  pillar                 1.3.1     2018-12-15 [1] CRAN (R 3.5.1)
#  pkgconfig              2.0.2     2018-08-16 [1] CRAN (R 3.5.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.5.0)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.5.0)
#  promises               1.0.1     2018-04-13 [2] CRAN (R 3.5.0)
#  purrr                  0.3.1     2019-03-03 [2] CRAN (R 3.5.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.5.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.5.0)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.5.0)
#  Rcpp                   1.0.0     2018-11-07 [1] CRAN (R 3.5.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.5.1)
#  rlang                  0.3.1     2019-01-08 [1] CRAN (R 3.5.1)
#  rmote                * 0.3.4     2018-05-02 [1] deltarho (R 3.5.0)
#  S4Vectors            * 0.20.1    2018-11-09 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.5.1)
#  segmented              0.5-3.0   2017-11-30 [2] CRAN (R 3.5.0)
#  servr                  0.13      2019-03-04 [1] CRAN (R 3.5.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.5.1)
#  statmod                1.4.30    2017-06-18 [2] CRAN (R 3.5.0)
#  SummarizedExperiment * 1.12.0    2018-10-30 [1] Bioconductor
#  tibble                 2.0.1     2019-01-12 [1] CRAN (R 3.5.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.5.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.5.0)
#  xfun                   0.5       2019-02-20 [1] CRAN (R 3.5.1)
#  XVector                0.22.0    2018-10-30 [1] Bioconductor
#  zlibbioc               1.28.0    2018-10-30 [2] Bioconductor
#
# [1] /users/lcollado/R/x86_64-pc-linux-gnu-library/3.5.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda-3/envs/svnR-3.5.x/R/3.5.x/lib64/R/library
