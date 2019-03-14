library('clusterProfiler')
library('gplots')
library('GenomicRanges')
library('devtools')
library('VennDiagram')
library('RColorBrewer')
library('ggplot2')
library('jaffelab')
library('limma')
library('edgeR')
library('SummarizedExperiment')

## Load case-control results
files <- c(
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_geneLevel.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_geneLevel.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_geneLevel_noHGoldQSV_matchHIPPO.rda',
    '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_geneLevel_noHGoldQSV_matchDLPFC.rda'
)

outGene <- lapply(files, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(outGene)
})

outGene0 <- lapply(files, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(outGene0)
})

outGeneNoAdj <- lapply(files, function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    return(outGeneNoAdj)
})


names(outGeneNoAdj) <- names(outGene0) <- names(outGene) <- c('HIPPO_allQSV', 'DLPFC_allQSV', 'HIPPO_noHGoldQSV', 'DLPFC_noHGoldQSV', 'HIPPO_matchQSV', 'DLPFC_matchQSV')

## Load BrainSeq Phase 1 and Common Mind results
load('/users/ajaffe/Lieber/Projects/RNAseq/firstRnaSeqPaper/caseControl/rdas/expressed_de_features.rda', verbose = TRUE)

prev <- list(
    'BSP1' = data.frame(
        ensemblID = names(outStatsExprs$Gene),
        adj.P.Val = outStatsExprs$Gene$fdr_qsva,
        logFC = outStatsExprs$Gene$log2FC_qsva,
        t = outStatsExprs$Gene$tstat_qsva
        ),
    'CMC' = data.frame(
        ensemblID = names(outStatsExprs$Gene),
        adj.P.Val = p.adjust(outStatsExprs$Gene$CMC_pval_qsva, method = 'fdr'),
        logFC = outStatsExprs$Gene$CMC_log2FC_qsva,
        t = outStatsExprs$Gene$CMC_tstat_qsva
        )
)

outGene <- c(outGene, prev)

## Load exon/jx/tx results
outFeat <- lapply(c('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda', '/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_hippo_filtered_qSVA_noHGoldQSV_matchHIPPO.rda'), function(f) {
    message(paste(Sys.time(), 'loading', f))
    load(f, verbose = TRUE)
    outTx$ensemblID <- gsub('\\..*', '', outTx$gene_id)
    return(list('exon' = outExon, 'jxn' = outJxn, 'tx' = outTx))
})
names(outFeat) <- c('DLPFC', 'HIPPO')
# > sapply(outFeat, function(x) { sapply(x, nrow) })
#       DLPFC  HIPPO
# exon 396583 396583
# jxn  297181 297181
# tx    92732  92732


## Get number of DE genes at different FDR cutoffs
n_de <- do.call(rbind, lapply(c(0.05, 0.1, 0.15, 0.2), function(cut) {
    xx <- sapply(c(outGene, do.call(c, outFeat)), function(x) {
        table(factor(x$adj.P.Val < cut, levels = c('FALSE', 'TRUE')))
    })
    cbind(xx, cutoff = cut)
}))

options(width = 200)
n_de
#       HIPPO_allQSV DLPFC_allQSV HIPPO_noHGoldQSV DLPFC_noHGoldQSV HIPPO_matchQSV DLPFC_matchQSV  BSP1   CMC DLPFC.exon DLPFC.jxn DLPFC.tx HIPPO.exon HIPPO.jxn HIPPO.tx cutoff
# FALSE        24652        24307            24635            24190          24604          24407 23939 23746     396143    297144    92726     396386    297140    92732   0.05
# TRUE             0          345               17              462             48            245   183   376        440        37        6        197        41        0   0.05
# FALSE        24648        23923            24592            23800          24551          24020 23616 23220     394775    297091    92712     396007    297116    92731   0.10
# TRUE             4          729               60              852            101            632   506   902       1808        90       20        576        65        1   0.10
# FALSE        24601        23543            24502            23360          24489          23618 23205 22661     393038    296960    92699     395442    297025    92709   0.15
# TRUE            51         1109              150             1292            163           1034   917  1461       3545       221       33       1141       156       23   0.15
# FALSE        24524        23112            24369            22857          24320          23103 22736 22002     390412    296762    92671     394224    296905    92643   0.20
# TRUE           128         1540              283             1795            332           1549  1386  2120       6171       419       61       2359       276       89   0.20

data.frame(colnames(n_de))
#     colnames.n_de.
# 1     HIPPO_allQSV
# 2     DLPFC_allQSV
# 3 HIPPO_noHGoldQSV
# 4 DLPFC_noHGoldQSV
# 5   HIPPO_matchQSV
# 6   DLPFC_matchQSV
# 7             BSP1
# 8              CMC
# 9           cutoff

n_de_sign <- do.call(rbind, lapply(c(0.05, 0.1, 0.15, 0.2), function(cut) {
    xx <- lapply(c(outGene[5:8], do.call(c, outFeat)), function(x) {
        y <- table(factor(x$adj.P.Val < cut, levels = c('FALSE', 'TRUE')), factor(sign(x$logFC), levels = c(-1, 0, 1)))
        data.frame(de_status = rownames(y), n = as.vector(y), sign = rep(colnames(y), each = 2), cutoff = cut)
    })
    for(i in seq_len(length(xx))) { xx[[i]]$model = names(xx)[i] }
    names(xx) <- NULL
    do.call(rbind, xx)
}))
n_de_sign$group <- ifelse(n_de_sign$sign == 0, 'none', ifelse(n_de_sign$sign == -1, 'Control', 'SCZD'))
# n_de_sign
#     de_status      n sign cutoff          model   group
# 1       FALSE  12825   -1   0.05 HIPPO_matchQSV Control
# 2        TRUE     27   -1   0.05 HIPPO_matchQSV Control
# 3       FALSE      0    0   0.05 HIPPO_matchQSV    none
# 4        TRUE      0    0   0.05 HIPPO_matchQSV    none
# 5       FALSE  11779    1   0.05 HIPPO_matchQSV    SCZD
# 6        TRUE     21    1   0.05 HIPPO_matchQSV    SCZD
# 7       FALSE  13065   -1   0.05 DLPFC_matchQSV Control
# 8        TRUE    142   -1   0.05 DLPFC_matchQSV Control
# 9       FALSE      0    0   0.05 DLPFC_matchQSV    none
# 10       TRUE      0    0   0.05 DLPFC_matchQSV    none
# 11      FALSE  11342    1   0.05 DLPFC_matchQSV    SCZD
# 12       TRUE    103    1   0.05 DLPFC_matchQSV    SCZD
# 13      FALSE  11415   -1   0.05           BSP1 Control
# 14       TRUE     50   -1   0.05           BSP1 Control
# 15      FALSE      0    0   0.05           BSP1    none
# 16       TRUE      0    0   0.05           BSP1    none
# 17      FALSE  12524    1   0.05           BSP1    SCZD
# 18       TRUE    133    1   0.05           BSP1    SCZD
# 19      FALSE  12101   -1   0.05            CMC Control
# 20       TRUE    243   -1   0.05            CMC Control
# 21      FALSE     21    0   0.05            CMC    none
# 22       TRUE      0    0   0.05            CMC    none
# 23      FALSE  11624    1   0.05            CMC    SCZD
# 24       TRUE    133    1   0.05            CMC    SCZD
# 25      FALSE 206717   -1   0.05     DLPFC.exon Control
# 26       TRUE    215   -1   0.05     DLPFC.exon Control
# 27      FALSE      0    0   0.05     DLPFC.exon    none
# 28       TRUE      0    0   0.05     DLPFC.exon    none
# 29      FALSE 189426    1   0.05     DLPFC.exon    SCZD
# 30       TRUE    225    1   0.05     DLPFC.exon    SCZD
# 31      FALSE 143266   -1   0.05      DLPFC.jxn Control
# 32       TRUE     20   -1   0.05      DLPFC.jxn Control
# 33      FALSE      0    0   0.05      DLPFC.jxn    none
# 34       TRUE      0    0   0.05      DLPFC.jxn    none
# 35      FALSE 153878    1   0.05      DLPFC.jxn    SCZD
# 36       TRUE     17    1   0.05      DLPFC.jxn    SCZD
# 37      FALSE  44746   -1   0.05       DLPFC.tx Control
# 38       TRUE      4   -1   0.05       DLPFC.tx Control
# 39      FALSE      0    0   0.05       DLPFC.tx    none
# 40       TRUE      0    0   0.05       DLPFC.tx    none
# 41      FALSE  47980    1   0.05       DLPFC.tx    SCZD
# 42       TRUE      2    1   0.05       DLPFC.tx    SCZD
# 43      FALSE 194911   -1   0.05     HIPPO.exon Control
# 44       TRUE     56   -1   0.05     HIPPO.exon Control
# 45      FALSE      0    0   0.05     HIPPO.exon    none
# 46       TRUE      0    0   0.05     HIPPO.exon    none
# 47      FALSE 201475    1   0.05     HIPPO.exon    SCZD
# 48       TRUE    141    1   0.05     HIPPO.exon    SCZD
# 49      FALSE 152287   -1   0.05      HIPPO.jxn Control
# 50       TRUE      9   -1   0.05      HIPPO.jxn Control
# 51      FALSE      0    0   0.05      HIPPO.jxn    none
# 52       TRUE      0    0   0.05      HIPPO.jxn    none
# 53      FALSE 144853    1   0.05      HIPPO.jxn    SCZD
# 54       TRUE     32    1   0.05      HIPPO.jxn    SCZD
# 55      FALSE  39554   -1   0.05       HIPPO.tx Control
# 56       TRUE      0   -1   0.05       HIPPO.tx Control
# 57      FALSE      0    0   0.05       HIPPO.tx    none
# 58       TRUE      0    0   0.05       HIPPO.tx    none
# 59      FALSE  53178    1   0.05       HIPPO.tx    SCZD
# 60       TRUE      0    1   0.05       HIPPO.tx    SCZD
# 61      FALSE  12800   -1   0.10 HIPPO_matchQSV Control
# 62       TRUE     52   -1   0.10 HIPPO_matchQSV Control
# 63      FALSE      0    0   0.10 HIPPO_matchQSV    none
# 64       TRUE      0    0   0.10 HIPPO_matchQSV    none
# 65      FALSE  11751    1   0.10 HIPPO_matchQSV    SCZD
# 66       TRUE     49    1   0.10 HIPPO_matchQSV    SCZD
# 67      FALSE  12828   -1   0.10 DLPFC_matchQSV Control
# 68       TRUE    379   -1   0.10 DLPFC_matchQSV Control
# 69      FALSE      0    0   0.10 DLPFC_matchQSV    none
# 70       TRUE      0    0   0.10 DLPFC_matchQSV    none
# 71      FALSE  11192    1   0.10 DLPFC_matchQSV    SCZD
# 72       TRUE    253    1   0.10 DLPFC_matchQSV    SCZD
# 73      FALSE  11288   -1   0.10           BSP1 Control
# 74       TRUE    177   -1   0.10           BSP1 Control
# 75      FALSE      0    0   0.10           BSP1    none
# 76       TRUE      0    0   0.10           BSP1    none
# 77      FALSE  12328    1   0.10           BSP1    SCZD
# 78       TRUE    329    1   0.10           BSP1    SCZD
# 79      FALSE  11796   -1   0.10            CMC Control
# 80       TRUE    548   -1   0.10            CMC Control
# 81      FALSE     21    0   0.10            CMC    none
# 82       TRUE      0    0   0.10            CMC    none
# 83      FALSE  11403    1   0.10            CMC    SCZD
# 84       TRUE    354    1   0.10            CMC    SCZD
# 85      FALSE 205872   -1   0.10     DLPFC.exon Control
# 86       TRUE   1060   -1   0.10     DLPFC.exon Control
# 87      FALSE      0    0   0.10     DLPFC.exon    none
# 88       TRUE      0    0   0.10     DLPFC.exon    none
# 89      FALSE 188903    1   0.10     DLPFC.exon    SCZD
# 90       TRUE    748    1   0.10     DLPFC.exon    SCZD
# 91      FALSE 143242   -1   0.10      DLPFC.jxn Control
# 92       TRUE     44   -1   0.10      DLPFC.jxn Control
# 93      FALSE      0    0   0.10      DLPFC.jxn    none
# 94       TRUE      0    0   0.10      DLPFC.jxn    none
# 95      FALSE 153849    1   0.10      DLPFC.jxn    SCZD
# 96       TRUE     46    1   0.10      DLPFC.jxn    SCZD
# 97      FALSE  44740   -1   0.10       DLPFC.tx Control
# 98       TRUE     10   -1   0.10       DLPFC.tx Control
# 99      FALSE      0    0   0.10       DLPFC.tx    none
# 100      TRUE      0    0   0.10       DLPFC.tx    none
# 101     FALSE  47972    1   0.10       DLPFC.tx    SCZD
# 102      TRUE     10    1   0.10       DLPFC.tx    SCZD
# 103     FALSE 194789   -1   0.10     HIPPO.exon Control
# 104      TRUE    178   -1   0.10     HIPPO.exon Control
# 105     FALSE      0    0   0.10     HIPPO.exon    none
# 106      TRUE      0    0   0.10     HIPPO.exon    none
# 107     FALSE 201218    1   0.10     HIPPO.exon    SCZD
# 108      TRUE    398    1   0.10     HIPPO.exon    SCZD
# 109     FALSE 152282   -1   0.10      HIPPO.jxn Control
# 110      TRUE     14   -1   0.10      HIPPO.jxn Control
# 111     FALSE      0    0   0.10      HIPPO.jxn    none
# 112      TRUE      0    0   0.10      HIPPO.jxn    none
# 113     FALSE 144834    1   0.10      HIPPO.jxn    SCZD
# 114      TRUE     51    1   0.10      HIPPO.jxn    SCZD
# 115     FALSE  39554   -1   0.10       HIPPO.tx Control
# 116      TRUE      0   -1   0.10       HIPPO.tx Control
# 117     FALSE      0    0   0.10       HIPPO.tx    none
# 118      TRUE      0    0   0.10       HIPPO.tx    none
# 119     FALSE  53177    1   0.10       HIPPO.tx    SCZD
# 120      TRUE      1    1   0.10       HIPPO.tx    SCZD
# 121     FALSE  12768   -1   0.15 HIPPO_matchQSV Control
# 122      TRUE     84   -1   0.15 HIPPO_matchQSV Control
# 123     FALSE      0    0   0.15 HIPPO_matchQSV    none
# 124      TRUE      0    0   0.15 HIPPO_matchQSV    none
# 125     FALSE  11721    1   0.15 HIPPO_matchQSV    SCZD
# 126      TRUE     79    1   0.15 HIPPO_matchQSV    SCZD
# 127     FALSE  12583   -1   0.15 DLPFC_matchQSV Control
# 128      TRUE    624   -1   0.15 DLPFC_matchQSV Control
# 129     FALSE      0    0   0.15 DLPFC_matchQSV    none
# 130      TRUE      0    0   0.15 DLPFC_matchQSV    none
# 131     FALSE  11035    1   0.15 DLPFC_matchQSV    SCZD
# 132      TRUE    410    1   0.15 DLPFC_matchQSV    SCZD
# 133     FALSE  11126   -1   0.15           BSP1 Control
# 134      TRUE    339   -1   0.15           BSP1 Control
# 135     FALSE      0    0   0.15           BSP1    none
# 136      TRUE      0    0   0.15           BSP1    none
# 137     FALSE  12079    1   0.15           BSP1    SCZD
# 138      TRUE    578    1   0.15           BSP1    SCZD
# 139     FALSE  11491   -1   0.15            CMC Control
# 140      TRUE    853   -1   0.15            CMC Control
# 141     FALSE     21    0   0.15            CMC    none
# 142      TRUE      0    0   0.15            CMC    none
# 143     FALSE  11149    1   0.15            CMC    SCZD
# 144      TRUE    608    1   0.15            CMC    SCZD
# 145     FALSE 204863   -1   0.15     DLPFC.exon Control
# 146      TRUE   2069   -1   0.15     DLPFC.exon Control
# 147     FALSE      0    0   0.15     DLPFC.exon    none
# 148      TRUE      0    0   0.15     DLPFC.exon    none
# 149     FALSE 188175    1   0.15     DLPFC.exon    SCZD
# 150      TRUE   1476    1   0.15     DLPFC.exon    SCZD
# 151     FALSE 143181   -1   0.15      DLPFC.jxn Control
# 152      TRUE    105   -1   0.15      DLPFC.jxn Control
# 153     FALSE      0    0   0.15      DLPFC.jxn    none
# 154      TRUE      0    0   0.15      DLPFC.jxn    none
# 155     FALSE 153779    1   0.15      DLPFC.jxn    SCZD
# 156      TRUE    116    1   0.15      DLPFC.jxn    SCZD
# 157     FALSE  44734   -1   0.15       DLPFC.tx Control
# 158      TRUE     16   -1   0.15       DLPFC.tx Control
# 159     FALSE      0    0   0.15       DLPFC.tx    none
# 160      TRUE      0    0   0.15       DLPFC.tx    none
# 161     FALSE  47965    1   0.15       DLPFC.tx    SCZD
# 162      TRUE     17    1   0.15       DLPFC.tx    SCZD
# 163     FALSE 194588   -1   0.15     HIPPO.exon Control
# 164      TRUE    379   -1   0.15     HIPPO.exon Control
# 165     FALSE      0    0   0.15     HIPPO.exon    none
# 166      TRUE      0    0   0.15     HIPPO.exon    none
# 167     FALSE 200854    1   0.15     HIPPO.exon    SCZD
# 168      TRUE    762    1   0.15     HIPPO.exon    SCZD
# 169     FALSE 152256   -1   0.15      HIPPO.jxn Control
# 170      TRUE     40   -1   0.15      HIPPO.jxn Control
# 171     FALSE      0    0   0.15      HIPPO.jxn    none
# 172      TRUE      0    0   0.15      HIPPO.jxn    none
# 173     FALSE 144769    1   0.15      HIPPO.jxn    SCZD
# 174      TRUE    116    1   0.15      HIPPO.jxn    SCZD
# 175     FALSE  39550   -1   0.15       HIPPO.tx Control
# 176      TRUE      4   -1   0.15       HIPPO.tx Control
# 177     FALSE      0    0   0.15       HIPPO.tx    none
# 178      TRUE      0    0   0.15       HIPPO.tx    none
# 179     FALSE  53159    1   0.15       HIPPO.tx    SCZD
# 180      TRUE     19    1   0.15       HIPPO.tx    SCZD
# 181     FALSE  12681   -1   0.20 HIPPO_matchQSV Control
# 182      TRUE    171   -1   0.20 HIPPO_matchQSV Control
# 183     FALSE      0    0   0.20 HIPPO_matchQSV    none
# 184      TRUE      0    0   0.20 HIPPO_matchQSV    none
# 185     FALSE  11639    1   0.20 HIPPO_matchQSV    SCZD
# 186      TRUE    161    1   0.20 HIPPO_matchQSV    SCZD
# 187     FALSE  12286   -1   0.20 DLPFC_matchQSV Control
# 188      TRUE    921   -1   0.20 DLPFC_matchQSV Control
# 189     FALSE      0    0   0.20 DLPFC_matchQSV    none
# 190      TRUE      0    0   0.20 DLPFC_matchQSV    none
# 191     FALSE  10817    1   0.20 DLPFC_matchQSV    SCZD
# 192      TRUE    628    1   0.20 DLPFC_matchQSV    SCZD
# 193     FALSE  10949   -1   0.20           BSP1 Control
# 194      TRUE    516   -1   0.20           BSP1 Control
# 195     FALSE      0    0   0.20           BSP1    none
# 196      TRUE      0    0   0.20           BSP1    none
# 197     FALSE  11787    1   0.20           BSP1    SCZD
# 198      TRUE    870    1   0.20           BSP1    SCZD
# 199     FALSE  11132   -1   0.20            CMC Control
# 200      TRUE   1212   -1   0.20            CMC Control
# 201     FALSE     21    0   0.20            CMC    none
# 202      TRUE      0    0   0.20            CMC    none
# 203     FALSE  10849    1   0.20            CMC    SCZD
# 204      TRUE    908    1   0.20            CMC    SCZD
# 205     FALSE 203413   -1   0.20     DLPFC.exon Control
# 206      TRUE   3519   -1   0.20     DLPFC.exon Control
# 207     FALSE      0    0   0.20     DLPFC.exon    none
# 208      TRUE      0    0   0.20     DLPFC.exon    none
# 209     FALSE 186999    1   0.20     DLPFC.exon    SCZD
# 210      TRUE   2652    1   0.20     DLPFC.exon    SCZD
# 211     FALSE 143068   -1   0.20      DLPFC.jxn Control
# 212      TRUE    218   -1   0.20      DLPFC.jxn Control
# 213     FALSE      0    0   0.20      DLPFC.jxn    none
# 214      TRUE      0    0   0.20      DLPFC.jxn    none
# 215     FALSE 153694    1   0.20      DLPFC.jxn    SCZD
# 216      TRUE    201    1   0.20      DLPFC.jxn    SCZD
# 217     FALSE  44721   -1   0.20       DLPFC.tx Control
# 218      TRUE     29   -1   0.20       DLPFC.tx Control
# 219     FALSE      0    0   0.20       DLPFC.tx    none
# 220      TRUE      0    0   0.20       DLPFC.tx    none
# 221     FALSE  47950    1   0.20       DLPFC.tx    SCZD
# 222      TRUE     32    1   0.20       DLPFC.tx    SCZD
# 223     FALSE 194097   -1   0.20     HIPPO.exon Control
# 224      TRUE    870   -1   0.20     HIPPO.exon Control
# 225     FALSE      0    0   0.20     HIPPO.exon    none
# 226      TRUE      0    0   0.20     HIPPO.exon    none
# 227     FALSE 200127    1   0.20     HIPPO.exon    SCZD
# 228      TRUE   1489    1   0.20     HIPPO.exon    SCZD
# 229     FALSE 152209   -1   0.20      HIPPO.jxn Control
# 230      TRUE     87   -1   0.20      HIPPO.jxn Control
# 231     FALSE      0    0   0.20      HIPPO.jxn    none
# 232      TRUE      0    0   0.20      HIPPO.jxn    none
# 233     FALSE 144696    1   0.20      HIPPO.jxn    SCZD
# 234      TRUE    189    1   0.20      HIPPO.jxn    SCZD
# 235     FALSE  39535   -1   0.20       HIPPO.tx Control
# 236      TRUE     19   -1   0.20       HIPPO.tx Control
# 237     FALSE      0    0   0.20       HIPPO.tx    none
# 238      TRUE      0    0   0.20       HIPPO.tx    none
# 239     FALSE  53108    1   0.20       HIPPO.tx    SCZD
# 240      TRUE     70    1   0.20       HIPPO.tx    SCZD

pdf('pdf/n_de_sign.pdf', useDingbats = FALSE, width = 10, height = 18)
ggplot(n_de_sign, aes(x = de_status, y = n, fill = group)) + geom_bar(stat = 'identity', width = 0.5, position = 'dodge') + facet_grid(model ~ cutoff) + theme_bw(base_size = 18)
ggplot(subset(n_de_sign, de_status == 'TRUE'), aes(x = de_status, y = n, fill = group)) + geom_bar(stat = 'identity', width = 0.5, position = 'dodge') + facet_grid(model ~ cutoff) + theme_bw(base_size = 18)
ggplot(subset(n_de_sign, de_status == 'TRUE' & sign != 0), aes(x = de_status, y = n, fill = group)) + geom_bar(stat = 'identity', width = 0.5, position = 'dodge') + facet_grid(model ~ cutoff) + theme_bw(base_size = 18)
dev.off()

## Explore number of DE genes across models
make_venn <- function(i, txt = 'HIPPO_', cut = 0.1) {
    vinfo <- lapply(c(outGene, do.call(c, outFeat))[i], function(x) {
        x$ensemblID[x$adj.P.Val < cut]
    })
    names(vinfo) <- gsub('_match', '.gene', gsub(txt, '', names(vinfo)))
    venn(vinfo) + title(paste('FDR cutoff:', cut))
}

make_venn2 <- function(i, txt = 'QSV|IPPO|LPFC') {
    make_venn(i, txt = txt, cut = 0.05)
    make_venn(i, txt = txt)
}

data.frame(name = names(c(outGene, do.call(c, outFeat))))
#                name
# 1      HIPPO_allQSV
# 2      DLPFC_allQSV
# 3  HIPPO_noHGoldQSV
# 4  DLPFC_noHGoldQSV
# 5    HIPPO_matchQSV
# 6    DLPFC_matchQSV
# 7              BSP1
# 8               CMC
# 9        DLPFC.exon
# 10        DLPFC.jxn
# 11         DLPFC.tx
# 12       HIPPO.exon
# 13        HIPPO.jxn
# 14         HIPPO.tx

pdf('pdf/venn_across_models.pdf', useDingbats = FALSE)
make_venn2(c(1, 3, 5))
make_venn2(c(2, 4, 6))
make_venn2(3:6)
make_venn2(c(4, 6, 7, 8))
make_venn2(5:6)
make_venn2(c(5, 6, 7, 8))
make_venn2(c(6, 9:10))
make_venn2(c(6, 9:11))
make_venn2(c(5, 12:13))
make_venn2(c(5, 12:14))
make_venn2(c(5:6, 9:10))
make_venn2(c(5:6, 9, 12))
make_venn2(c(5:6, 10, 13))
make_venn2(c(9, 12, 10, 13))
dev.off()

## DE genes by sign
de_genes_sign <- mapply(function(x, cut, sign) {
    x$ensemblID[x$adj.P.Val < cut & sign(x$logFC) == sign]
}, c(outGene[rep(5:6, each = 2)], do.call(c, outFeat)[rep(1:6, each = 2)]), c(rep(c(0.2, 0.1), each = 2), rep(c(0.1, 0.2), each = 6)), sign = c(-1, 1))
## Sign -1 corresponds to control, +1 to SCZD
# load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata', verbose = TRUE)
# > colnames(modQsva)[grep('Dx', colnames(modQsva))]
# [1] "DxSchizo"
names(de_genes_sign) <- paste0(gsub('HIPPO', 'H', gsub('DLPFC', 'D', gsub('_matchQSV', '.gene', names(de_genes_sign)))), '_', c('Control', 'SCZD'))

xx <- as.data.frame(sapply(de_genes_sign, function(x) length(unique(x[!is.na(x)]))))
colnames(xx) <- 'unique_gene_id'
xx$i <- seq_len(nrow(xx))
xx$n_noNA <- sapply(de_genes_sign, function(x) length(x[!is.na(x)]))
xx
#                unique_gene_id  i n_noNA
# H.gene_Control            171  1    171
# H.gene_SCZD               161  2    161
# D.gene_Control            379  3    379
# D.gene_SCZD               253  4    253
# D.exon_Control            412  5   1060
# D.exon_SCZD               292  6    748
# D.jxn_Control              28  7     36
# D.jxn_SCZD                 29  8     34
# D.tx_Control               10  9     10
# D.tx_SCZD                  10 10     10
# H.exon_Control            346 11    870
# H.exon_SCZD               629 12   1489
# H.jxn_Control              55 13     68
# H.jxn_SCZD                111 14    155
# H.tx_Control               19 15     19
# H.tx_SCZD                  70 16     70


## Top DE genes by sign
de_genes_sign_top <- mapply(function(x, sign) {
    x <- x[sign(x$logFC) == sign, ]
    head(x$ensemblID[order(x$adj.P.Val, decreasing = FALSE)], 400)
}, c(outGene[c(5, 5, 6, 6)], do.call(c, outFeat)[rep(1:6, each = 2)]), sign = rep(c(-1, 1), 8), SIMPLIFY = FALSE)
names(de_genes_sign_top) <- names(de_genes_sign)
as.data.frame(sapply(de_genes_sign_top, length))

## Pretty venn code
#venn_cols <- brewer.pal('Set1', n = 4)
venn_cols <- c('skyblue3', 'dark orange', 'red', 'purple')
names(venn_cols) <- paste0(rep(c('DLPFC', 'HIPPO'), each = 2), '_', c('control', 'SCZD'))
make_venn_pretty <- function(genes, title = 'DE features grouped by gene id') {
    genes <- lapply(genes, function(x) x[!is.na(x)])
    v <- venn.diagram(genes, filename = NULL,
        main = title,
        col = 'transparent', fill = venn_cols[seq_len(length(genes))],
        alpha = 0.5, margin = 0,
        main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
        cat.col = venn_cols[seq_len(length(genes))], scaled = FALSE, cat.pos = 0)
    grid.newpage()
    grid.draw(v)
}

pdf('pdf/venn_de_genes_by_sign.pdf', useDingbats = FALSE)
## gene
make_venn_pretty(de_genes_sign[1:4], 'Gene DLPFC FDR10%, HIPPO FDR20%')
make_venn_pretty(list('HIPPO' = de_genes_sign[[1]], 'DLPFC' = de_genes_sign[[3]]), 'Gene Control (DLPFC FDR10%, HIPPO FDR20%)')
make_venn_pretty(list('HIPPO' = de_genes_sign[[2]], 'DLPFC' = de_genes_sign[[4]]), 'Gene SCZD (DLPFC FDR10%, HIPPO FDR20%)')

## exon
make_venn_pretty(de_genes_sign[c(11:12, 5:6)], 'Exon DLPFC FDR10%, HIPPO FDR20%')
make_venn_pretty(list('HIPPO' = de_genes_sign[[11]], 'DLPFC' = de_genes_sign[[5]]), 'Exon Control (DLPFC FDR10%, HIPPO FDR20%)')
make_venn_pretty(list('HIPPO' = de_genes_sign[[12]], 'DLPFC' = de_genes_sign[[6]]), 'Exon SCZD (DLPFC FDR10%, HIPPO FDR20%)')

## jxn
make_venn_pretty(de_genes_sign[c(13:14, 7:8)], 'JXN DLPFC FDR10%, HIPPO FDR20%')
make_venn_pretty(list('HIPPO' = de_genes_sign[[13]], 'DLPFC' = de_genes_sign[[7]]), 'JXN Control (DLPFC FDR10%, HIPPO FDR20%)')
make_venn_pretty(list('HIPPO' = de_genes_sign[[14]], 'DLPFC' = de_genes_sign[[8]]), 'JXN SCZD (DLPFC FDR10%, HIPPO FDR20%)')

## tx
make_venn_pretty(de_genes_sign[c(15:16, 9:10)], 'TX DLPFC FDR10%, HIPPO FDR20%')
make_venn_pretty(list('HIPPO' = de_genes_sign[[15]], 'DLPFC' = de_genes_sign[[9]]), 'TX Control (DLPFC FDR10%, HIPPO FDR20%)')
make_venn_pretty(list('HIPPO' = de_genes_sign[[16]], 'DLPFC' = de_genes_sign[[10]]), 'TX SCZD (DLPFC FDR10%, HIPPO FDR20%)')

## Top 400 gene
make_venn_pretty(de_genes_sign_top[1:4], 'gene: top 400 in each group')
make_venn_pretty(list('HIPPO' = de_genes_sign_top[[1]], 'DLPFC' = de_genes_sign_top[[3]]), 'gene: Control (top 400)')
make_venn_pretty(list('HIPPO' = de_genes_sign_top[[2]], 'DLPFC' = de_genes_sign_top[[4]]), 'gene: SCZD (top 400)')

## exon
make_venn_pretty(de_genes_sign_top[c(11:12, 5:6)], 'exon: top 400 in each group')
make_venn_pretty(list('HIPPO' = de_genes_sign_top[[11]], 'DLPFC' = de_genes_sign_top[[5]]), 'exon: Control (top 400)')
make_venn_pretty(list('HIPPO' = de_genes_sign[[12]], 'DLPFC' = de_genes_sign_top[[6]]), 'exon: SCZD (top 400)')

## jxn
make_venn_pretty(de_genes_sign_top[c(13:14, 7:8)], 'jxn: top 400 in each group')
make_venn_pretty(list('HIPPO' = de_genes_sign_top[[13]], 'DLPFC' = de_genes_sign_top[[7]]), 'jxn: Control (top 400)')
make_venn_pretty(list('HIPPO' = de_genes_sign_top[[14]], 'DLPFC' = de_genes_sign_top[[8]]), 'jxn: SCZD (top 400)')

## tx
make_venn_pretty(de_genes_sign_top[c(15:16, 9:10)], 'tx: top 400 in each group')
make_venn_pretty(list('HIPPO' = de_genes_sign_top[[15]], 'DLPFC' = de_genes_sign_top[[9]]), 'tx: Control (top 400)')
make_venn_pretty(list('HIPPO' = de_genes_sign_top[[16]], 'DLPFC' = de_genes_sign_top[[10]]), 'tx: SCZD (top 400)')

## top 200 gene
make_venn_pretty(lapply(de_genes_sign_top[1:4], head, n = 200), 'gene: top 200 in each group')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[1]], 200), 'DLPFC' = head(de_genes_sign_top[[3]], 200)), 'gene: Control (top 200)')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[2]], 200), 'DLPFC' = head(de_genes_sign_top[[4]], 200)), 'gene: SCZD (top 200)')

make_venn_pretty(lapply(de_genes_sign_top[1:4], head, n = 150), 'gene: top 150 in each group')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[1]], 150), 'DLPFC' = head(de_genes_sign_top[[3]], 150)), 'gene: Control (top 150)')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[2]], 150), 'DLPFC' = head(de_genes_sign_top[[4]], 150)), 'gene: SCZD (top 150)')

make_venn_pretty(lapply(de_genes_sign_top[1:4], head, n = 100), 'gene: top 100 in each group')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[1]], 100), 'DLPFC' = head(de_genes_sign_top[[3]], 100)), 'gene: Control (top 100)')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[2]], 100), 'DLPFC' = head(de_genes_sign_top[[4]], 100)), 'gene: SCZD (top 100)')

make_venn_pretty(lapply(de_genes_sign_top[1:4], head, n = 50), 'gene: top 50 in each group')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[1]], 50), 'DLPFC' = head(de_genes_sign_top[[3]], 50)), 'gene: Control (top 50)')
make_venn_pretty(list('HIPPO' = head(de_genes_sign_top[[2]], 50), 'DLPFC' = head(de_genes_sign_top[[4]], 50)), 'gene: SCZD (top 50)')
dev.off()
system('rm VennDiagram*.log')


## Degradation results
load("/dcl01/ajaffe/data/lab/qsva_brain/ERs/rdas/DLPFC_Plus_HIPPO_RiboZero_geneLevel_degradationStats_forDEqual_hg38.rda", verbose = TRUE)
load("/dcl01/ajaffe/data/lab/qsva_brain/ERs/rdas/DLPFC_HIPPO_degradationStats_hg38.rda", verbose = TRUE)

## For DE_qual plots
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
                     function(x)
                       rgb(x[1], x[2], x[3], alpha=alpha))
}

plot_dequal <- function(out_input, degrade_input, var = 't', xlabtxt = 'case-control', main, ylabtxt = '') {
	both <- intersect(rownames(out_input), rownames(degrade_input))
	degrade <- degrade_input[both, ]
	interest <- out_input[both, ]

	stopifnot(identical(rownames(degrade), rownames(interest)))
	corr = signif(cor( degrade[, var],  interest[, var]), 3)
	plot(y = degrade[, var], x = interest[, var], xlab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), xlabtxt), ylab = paste(ylabtxt, var, 'degradation'), main = main, col = add.alpha('black', 1/10), pch = 16)
	legend('topleft', legend = paste('r =', corr))
}


## Make de_qual plots
plot_dequal2 <- function(i) {
    geneinfo <- list(outGeneNoAdj[[i]], outGeneNoAdj[[i]], outGene0[[i]], outGene0[[i]], outGene[[i]], outGene[[i]])
    makedequal <- function(v, gene, main, xlab) {
        plot_dequal(gene, degradeStats, var = v, xlabtxt = xlab, main = main, ylabtxt = 'Combined')
        plot_dequal(gene, degradeStats_HIPPO, var = v, xlabtxt = xlab, main = main, ylabtxt = 'HIPPO')
        plot_dequal(gene, degradeStats_DLPFC, var = v, xlabtxt = xlab, main = main, ylabtxt = 'DLPFC')
        plot_dequal(gene, degradeStatsInt, var = v, xlabtxt = xlab, main = main, ylabtxt = 'Combined adj interaction')
    }
    par(mfcol = c(4, 2))
    mapply(makedequal,
           v = rep(c('t', 'logFC'), 3),
           gene = geneinfo,
           main = rep(names(outGene)[i], 6),
           xlab = rep(c('case-control (Dx only)', 'case-control (without qSVs)', 'case-control (with qSVs)'), each = 2))
}


pdf('pdf/dequal_plots.pdf', useDingbats = FALSE, width = 8, height = 16)
for(i in 1:6) plot_dequal2(i)
dev.off()


## Gene ontology analysis
uni <- c(outGene[[3]]$ensemblID, outFeat[[1]][[1]]$ensemblID, outFeat[[1]][[2]]$ensemblID)
uni <- unique(uni[!is.na(uni)])
length(uni)
# [1] 27514

run_go <- function(genes, ont = c('BP', 'MF', 'CC')) {
    ## Change to ENSEMBL ids and remove NAs
    genes_ens <- lapply(lapply(genes, function(x) { gsub('\\..*', '', x) }), function(y) y[!is.na(y)])

    #genes_venn <- venn(genes_ens, show.plot = FALSE)

    ## Run GO analysis
    go_cluster <- lapply(ont, function(bp) {
        message(paste(Sys.time(), 'running GO analysis for', bp))
        tryCatch(compareCluster(genes_ens, fun = "enrichGO",
            universe = uni, OrgDb = 'org.Hs.eg.db',
            ont = bp, pAdjustMethod = "BH",
            pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
            readable = TRUE, keyType = 'ENSEMBL'),
            error = function(e) { return(NULL) })
    })
    names(go_cluster) <- ont
    
    message(paste(Sys.time(), 'running GO analysis for KEGG'))
    genes_ncbi <- lapply(lapply(genes_ens, bitr, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db'), function(x) x$ENTREZID)
    
    uni_ncbi <- bitr(uni, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID
    
    go_cluster$KEGG <- tryCatch(compareCluster(genes_ncbi, fun = 'enrichKEGG',
        universe = uni_ncbi, organism = 'hsa', pAdjustMethod = 'BH',
        pvalueCutoff = 0.1, qvalueCutoff = 0.05, keyType = 'ncbi-geneid'),
        error = function(e) { return(NULL) })
    
    return(go_cluster)
}


## Development DE genes
if(!file.exists('rdas/go_de_genes.Rdata')) {
    system.time( go_de_genes <- run_go(de_genes_sign[c(1:2, 11:14, 3:8)]) )
    message(paste(Sys.time(), 'saving rdas/go_de_genes.Rdata'))
    save(go_de_genes, file = 'rdas/go_de_genes.Rdata')
} else {
    message(paste(Sys.time(), 'loading rdas/go_de_genes.Rdata'))
    load('rdas/go_de_genes.Rdata', verbose = TRUE)
}
sapply(go_de_genes, class)

if(!file.exists('rdas/go_de_genes_top.Rdata')) {
    go_de_genes_top <- lapply(c(50, 100, 150, 200), function(n) {
        run_go(lapply(de_genes_sign_top[c(1:2, 11:14, 3:8)], head, n = n))
    })
    names(go_de_genes_top) <- c(50, 100, 150, 200)
    message(paste(Sys.time(), 'saving rdas/go_de_genes_top.Rdata'))
    save(go_de_genes_top, file = 'rdas/go_de_genes_top.Rdata')
} else {
    message(paste(Sys.time(), 'loading rdas/go_de_genes_top.Rdata'))
    load('rdas/go_de_genes_top.Rdata', verbose = TRUE)
}
lapply(go_de_genes_top, function(x) sapply(x, class))


simplify_go <- function(x) {
    #gsub('QSV|IPPO|LPFC', '', x)
    #gsub('_matchQSV', '', x)
    gsub('IPPO|LPFC|ontrol|chizo|CZD|xon_|xn_|x_|ene_|\\.', '', x)
}

plot_go <- function(go_cluster, cat = 10) {
    lapply(names(go_cluster), function(bp) {
        go <- go_cluster[[bp]]
        if(is.null(go)) {
            message(paste(Sys.time(), 'found no results for', bp))
            return(NULL)
        }

        ## Simplify names
        go@compareClusterResult$Cluster <- simplify_go(go@compareClusterResult$Cluster)
        names(go@geneClusters) <- simplify_go(names(go@geneClusters))

        print(plot(go, title = paste('ontology:', bp), font.size = 18, showCategory = cat, includeAll = TRUE))
        return(NULL)
    })
}


pdf('pdf/go_de_genes.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_go(go_de_genes)
dev.off()

pdf('pdf/go_all_de_genes.pdf', width = 16, height = 70, useDingbats = FALSE)
plot_go(go_de_genes, cat = NULL)
dev.off()

for(i in names(go_de_genes_top)) {
    pdf(paste0('pdf/go_de_genes_top', i, '.pdf'), width = 14, height = 9, useDingbats = FALSE)
    plot_go(go_de_genes_top[[i]])
    dev.off()

    pdf(paste0('pdf/go_all_de_genes_top', i, '.pdf'), width = 16, height = 25, useDingbats = FALSE)
    plot_go(go_de_genes_top[[i]], cat = NULL)
    dev.off()
}

run_gse <-  function(region, ont = c('BP', 'MF', 'CC')) {

    genes <- outGene[[paste0(region, '_matchQSV')]]$t
    names(genes) <- outGene[[paste0(region, '_matchQSV')]]$ensemblID
    genes <- genes[order(genes, decreasing = TRUE)]

    ## Run gseGO analysis
    go_cluster <- lapply(ont, function(bp) {
        message(paste(Sys.time(), 'running gseGO analysis for', bp))
        tryCatch(gseGO(genes, OrgDb = 'org.Hs.eg.db',
            ont = bp, pAdjustMethod = "BH",
            keyType = 'ENSEMBL', verbose = TRUE),
            error = function(e) { return(NULL) })
    })
    names(go_cluster) <- ont
    
    
    
    genes_tab <- bitr(names(genes), fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
    genes_tab <- genes_tab[!duplicated(genes_tab$ENTREZID), ]
        
    genes_m <- match(genes_tab$ENSEMBL, names(genes))
    genes_ncbi <- genes[genes_m]
    names(genes_ncbi) <- genes_tab$ENTREZID

    go_cluster$KEGG <- tryCatch(gseKEGG(genes_ncbi, organism = 'hsa',
        pAdjustMethod = "BH",
        keyType = 'ncbi-geneid', verbose = TRUE),
        error = function(e) { return(NULL) })
        
    return(go_cluster)
}


system.time( gse_hippo <- run_gse('HIPPO') )
system.time( gse_dlpfc <- run_gse('DLPFC') )
save(gse_hippo, gse_dlpfc, file = 'rdas/gse.Rdata')


plot_gse <- function(gse) {
    lapply(names(gse), function(bp) {
        go <- gse[[bp]]
        if(is.null(go)) {
            message(paste(Sys.time(), 'found no results for', bp))
            return(NULL)
        }

        print(dotplot(go, title = paste('ontology:', bp), font.size = 18))
        return(NULL)
    })
}


pdf('pdf/gse_hippo.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_gse(gse_hippo)
dev.off()

pdf('pdf/gse_dlpfc.pdf', width = 14, height = 9, useDingbats = FALSE)
plot_gse(gse_dlpfc)
dev.off()



## Scatter of logFC
comp_log <- function(x, y, xlab, ylab, var = 'logFC', de = FALSE, n = 150, onlyx = FALSE) {
    if(de) {
        if(!onlyx) {
            common <- unique(c(
                head(x$ensemblID[order(x$adj.P.Val, decreasing = FALSE)], n),
                head(y$ensemblID[order(y$adj.P.Val, decreasing = FALSE)], n)
            ))
        } else {
            common <- unique(c(
                head(x$ensemblID[order(x$adj.P.Val, decreasing = FALSE)], n)
            ))
        }
        
    } else {
        common <- intersect(x$ensemblID, y$ensemblID)
    }
    x <- x[match(common, x$ensemblID), ]
    y <- y[match(common, y$ensemblID), ]
    corr = signif(cor(x[, var], y[, var], use = 'pairwise.complete.obs'), 3)

    plot(x = x[, var], y = y[, var],
         xlab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), xlab),
         ylab = paste(ifelse(var == 't', 't-statistic', 'log2 FC'), ylab),
         col = add.alpha('black', ifelse(de, 1/2, 1/10)), pch = 16)
    legend('topleft', legend = paste('r =', corr))
    lines(loess.smooth(y = y[, var], x = x[, var]), col = 'red')
    abline(lm(y[, var] ~ x[, var]), col = 'blue')
    abline(h = 0, col = 'grey20')
    abline(v = 0, col = 'grey20')
}


pdf('pdf/scatter_models.pdf', useDingbats = FALSE)
comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC')
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1')
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC')
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1')
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC')

comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', var = 't')
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', var = 't')
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', var = 't')
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', var = 't')
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', var = 't')
dev.off()

pdf('pdf/scatter_models_top150de.pdf', useDingbats = FALSE)
comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', de = TRUE)
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', de = TRUE)
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', de = TRUE)
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', de = TRUE)
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', de = TRUE)

comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', var = 't', de = TRUE)
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', var = 't', de = TRUE)
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', var = 't', de = TRUE)
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', var = 't', de = TRUE)
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', var = 't', de = TRUE)
dev.off()

pdf('pdf/scatter_models_top400de.pdf', useDingbats = FALSE)
comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', de = TRUE, n = 400)
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', de = TRUE, n = 400)
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', de = TRUE, n = 400)
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', de = TRUE, n = 400)
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', de = TRUE, n = 400)

comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', var = 't', de = TRUE, n = 400)
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', var = 't', de = TRUE, n = 400)
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', var = 't', de = TRUE, n = 400)
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', var = 't', de = TRUE, n = 400)
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', var = 't', de = TRUE, n = 400)
dev.off()


pdf('pdf/scatter_models_top400de_onlyX.pdf', useDingbats = FALSE)
comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', de = TRUE, n = 400, onlyx = TRUE)

comp_log(outGene[[5]], outGene[[6]], 'HIPPO', 'DLPFC', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[5]], outGene[[7]], 'HIPPO', 'BSP1', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[5]], outGene[[8]], 'HIPPO', 'CMC', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[6]], outGene[[7]], 'DLPFC', 'BSP1', var = 't', de = TRUE, n = 400, onlyx = TRUE)
comp_log(outGene[[6]], outGene[[8]], 'DLPFC', 'CMC', var = 't', de = TRUE, n = 400, onlyx = TRUE)
dev.off()




## Load expression data
load("/dcl01/lieber/ajaffe/lab/brainseq_phase2/expr_cutoff/rse_gene.Rdata", verbose = TRUE)
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_HIPPO.Rdata', verbose = TRUE)
rse_gene_HIPPO <- rse_gene[, keepIndex]
mod_HIPPO <- mod
modQsva_HIPPO <- modQsva
load('/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/brainseq_phase2_qsvs_age17_noHGold_DLPFC.Rdata', verbose = TRUE)
rse_gene_DLPFC <- rse_gene[, keepIndex]
mod_DLPFC <- mod
modQsva_DLPFC <- modQsva

## Get normalized expression and clean it
cleaned <- mapply(function(expr, model) {
    dge = DGEList(counts = assays(expr)$counts,
	genes = rowData(expr))
    #calculate library-size adjustment
    dge = calcNormFactors(dge)
    vGene = voom(dge, model, plot=FALSE)
    list('norm' = vGene$E, 'cleaned' = cleaningY(vGene$E, model, P = 2))
},
    list('HIPPO' = rse_gene_HIPPO, 'DLPFC' = rse_gene_DLPFC),
    list('HIPPO' = modQsva_HIPPO, 'DLPFC' = modQsva_DLPFC),
    SIMPLIFY = FALSE
)

## Adapted some functions from
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/region_specific/explore_reg_specific.R
# and
# https://github.com/LieberInstitute/brainseq_phase2/blob/master/region_specific/explore_reg_specific_top.R
get_ylab <- function(type, cleaned = FALSE) {
    if(type %in% c('gene', 'exon', 'jxn')) {
        res <- 'log2(CPM + 0.5)'
    } else {
        res <- 'log2(TPM + 0.5)'
    }

    if(cleaned) res <- paste(res, '- covariate effects removed')
    return(res)
}

get_main <- function(i, region, group = 'schizo', type = 'gene') {
    if(type %in% c('gene', 'exon')) {
        var <- 'gencodeID'
        vars <- 'Symbol'
    } else if (type == 'jxn') {
        var <- 'gencodeGeneID'
        vars <- 'Symbol'
    } else {
        var <- 'gene_id'
        vars <- 'gene_name'
    }

    rse <- if(region == 'HIPPO') rse_gene_HIPPO else rse_gene_DLPFC

    topnow <- outGene[[paste0(region, '_matchQSV')]]

    j <- which(topnow$ensemblID == de_genes_sign_top[[paste0(region, '_', group)]][i])
    k <- which(names(rowRanges(rse)) == rownames(topnow)[j])

    paste(if(type != 'jxn') mcols(rowRanges(rse))[, var][k] else rownames(topnow)[j],
          if(is.na(mcols(rowRanges(rse))[, vars][k])) '' else mcols(rowRanges(rse))[, vars][k],
          'FDR',
          signif(topnow$adj.P.Val[j], 3),
          group)
}

get_ylim_mult <- function(rang) {
    c(
        ifelse(sign(rang[1]) == 1, 0.95, 1.05),
        ifelse(sign(rang[2]) == 1, 1.05, 0.95)
    )
}


plot_top_sign <- function(i, region, group = 'schizo', normtype = 'norm') {
    g <- de_genes_sign_top[[paste0(region, '_', group)]][i]
    j <- which(gsub('\\..*', '', rownames(cleaned[[region]][[normtype]])) == g)

    set.seed(20180426)
    dx <- if(region == 'HIPPO') colData(rse_gene_HIPPO)$Dx else colData(rse_gene_DLPFC)$Dx
    dx <- factor(dx, levels = c('Control', 'Schizo'))

    y <- cleaned[[region]][[normtype]][j, ]
    boxplot(y ~ dx,
            ylab = get_ylab('gene', cleaned = normtype == 'cleaned'),
            main = get_main(i, region, group),
            col = c('orchid1', 'aquamarine1'),
            ylim = abs(range(y)) * get_ylim_mult(range(y)) * sign(range(y)),
            outline = FALSE)
    points(y ~ jitter(as.integer(dx), 1), pch = 21,
           bg = c('orchid4', 'aquamarine4')[as.integer(dx)])
}

for(reg in c('HIPPO', 'DLPFC')) {
    pdf(paste0('pdf/top_50each_', reg, '_norm.pdf'))
    for(i in 1:50) {
        plot_top_sign(i, reg)
        plot_top_sign(i, reg, 'control')
    }
    dev.off()
    pdf(paste0('pdf/top_50each_', reg, '_cleaned.pdf'))
    for(i in 1:50) {
        plot_top_sign(i, reg, normtype = 'cleaned')
        plot_top_sign(i, reg, 'control', normtype = 'cleaned')
    }
    dev.off()
}


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
