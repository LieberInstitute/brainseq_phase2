get_proxy_info <- function(pos_hg19, prefix) {
    status <- rep('Other', length(pos_hg19))
    status[pos_hg19 %in% riskLoci$hg19POS2] <- 'Proxy'
    status[pos_hg19 %in% indexLoci$hg19POS] <- 'Index'
    print(table(status))
    
    indexSNP <- rep(NA, length(pos_hg19))
    m <- match(pos_hg19, riskLoci$hg19POS2)
    stopifnot(sum(is.na(m)) == sum(status == 'Other'))
    indexSNP[!is.na(m)] <- riskLoci$SNP1_Name[m[!is.na(m)]]
    
    indexSNP_pos_hg19 <- rep(NA, length(pos_hg19))
    indexSNP_pos_hg19[!is.na(m)] <- riskLoci$hg19POS1[m[!is.na(m)]]
    
    distance <- rep(NA, length(pos_hg19))
    distance[!is.na(m)] <- riskLoci$Distance[m[!is.na(m)]]
    
    
    
    res <- tibble(
        status = status,
        indexSNP = indexSNP,
        indexSNP_pos_hg19 = indexSNP_pos_hg19,
        indexSNP_distance = distance
    )
    
    print(length(unique(res$indexSNP[!is.na(res$indexSNP)])))
    print(summary(res$indexSNP_distance[!is.na(res$indexSNP_distance)]))
    
    colnames(res) <- paste0(prefix, colnames(res))
    return(res)
}

get_variable_by_region <- function(var, NAs_0 = FALSE) {
    
    m <- match(ttReg_map$ID, tt$ID)
    
    res <- data.frame(
        ID = ttReg_map$ID,
        feature = tt$feature[m],
        geneid = tt$geneid[m],
        genesymbol = tt$genesymbol[m],
        DLPFC = tt[ttReg_map$i_DLPFC, var, drop = TRUE],
        HIPPO = tt[ttReg_map$i_HIPPO, var, drop = TRUE],
        in_both = !is.na(ttReg_map$i_DLPFC) & !is.na(ttReg_map$i_HIPPO),
        TWAS.FDR_DLPFC = tt$TWAS.FDR[ttReg_map$i_DLPFC],
        TWAS.FDR_HIPPO = tt$TWAS.FDR[ttReg_map$i_HIPPO],
        TWAS.Bonf_DLPFC = tt$TWAS.Bonf[ttReg_map$i_DLPFC],
        TWAS.Bonf_HIPPO = tt$TWAS.Bonf[ttReg_map$i_HIPPO],
        BEST.GWAS.status_DLPFC = tt$BEST.GWAS.status[ttReg_map$i_DLPFC],
        BEST.GWAS.status_HIPPO = tt$BEST.GWAS.status[ttReg_map$i_HIPPO],
        stringsAsFactors = FALSE
    )
    
    res$FDR.5perc <- 'None'
    res$FDR.5perc[res$TWAS.FDR_DLPFC < 0.05] <- 'DLPFC'
    res$FDR.5perc[res$TWAS.FDR_HIPPO < 0.05] <- 'HIPPO'
    res$FDR.5perc[res$TWAS.FDR_DLPFC < 0.05 & res$TWAS.FDR_HIPPO < 0.05] <- 'Both'
    
    res$Bonf.5perc <- 'None'
    res$Bonf.5perc[res$TWAS.Bonf_DLPFC < 0.05] <- 'DLPFC'
    res$Bonf.5perc[res$TWAS.Bonf_HIPPO < 0.05] <- 'HIPPO'
    res$Bonf.5perc[res$TWAS.Bonf_DLPFC < 0.05 & res$TWAS.Bonf_HIPPO < 0.05] <- 'Both'

    res$BEST.GWAS.status <- 'Other'
    res$BEST.GWAS.status[c(which(res$BEST.GWAS.status_DLPFC != 'Other'),  which(res$BEST.GWAS.status_HIPPO != 'Other'))] <- 'Risk Locus'
    
    if(NAs_0 == TRUE) {
        res$DLPFC[is.na(res$DLPFC)] <- 0
        res$HIPPO[is.na(res$HIPPO)] <- 0
    }
    
    ## Make the features as factor, so its looks ok when plotting
    res$feature <- factor(res$feature, levels = c('gene', 'exon', 'jxn', 'tx'))
    res$FDR.5perc <- factor(res$FDR.5perc, levels = c('None', 'DLPFC', 'HIPPO', 'Both'))
    res$Bonf.5perc <- factor(res$Bonf.5perc, levels = c('None', 'DLPFC', 'HIPPO', 'Both'))
    res$BEST.GWAS.status <- factor(res$BEST.GWAS.status, levels = c('Other', 'Risk Locus'))
    
    return(res)
}

check_cor <- function(x, y) {
    cor(-log10(x), -log10(y))
}

create_gwas_or_eqtl <- function(tt_sigonly, filename = 'pdf/twas_fdr5perc_vs_gwas_or_eqtl.pdf', titleslug = 'FDR') {
    pdf(filename, useDingbats = FALSE, width = 21, height = 14)
    print(ggplot(tt_sigonly, aes(
        x = -log10(TWAS.P),
        y = -log10(EQTL.P.computed),
        color = BEST.GWAS.P.computed < 5e-08
    )) + geom_point() +
        facet_grid(region * 
            ifelse(EQTL.status == 'Other', 'Other', 'Risk Locus') ~
            factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
        ) +
        theme_bw(base_size = 30) +
        ggtitle(paste0('TWAS (', titleslug, '<5%) vs EQTL p-values')) +
        labs(caption = 'Risk Loci by EQTL')
    )

    print(ggplot(tt_sigonly, aes(
        x = -log10(TWAS.P),
        y = -log10(EQTL.P.computed),
        color = BEST.GWAS.P.computed < 5e-08
    )) + geom_point() +
        facet_grid(region * 
            ifelse(BEST.GWAS.status == 'Other', 'Other', 'Risk Locus') ~
            factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
        ) +
        theme_bw(base_size = 30) +
        ggtitle(paste0('TWAS (', titleslug, '<5%) vs EQTL p-values')) +
        labs(caption = 'Risk Loci by BEST GWAS')
    )
    
    print(ggplot(tt_sigonly, aes(
        x = TWAS.Z,
        y = EQTL.Z,
        color = BEST.GWAS.P.computed < 5e-08
    )) + geom_point() +
        facet_grid(region * 
            ifelse(EQTL.status == 'Other', 'Other', 'Risk Locus') ~
            factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
        ) +
        theme_bw(base_size = 30) +
        ggtitle(paste0('TWAS (', titleslug, '<5%) vs EQTL z-scores')) +
        labs(caption = 'Risk Loci by EQTL')
    )

    print(ggplot(tt_sigonly, aes(
        x = TWAS.Z,
        y = EQTL.Z,
        color = BEST.GWAS.P.computed < 5e-08
    )) + geom_point() +
        facet_grid(region * 
            ifelse(BEST.GWAS.status == 'Other', 'Other', 'Risk Locus') ~
            factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
        ) +
        theme_bw(base_size = 30) +
        ggtitle(paste0('TWAS (', titleslug, '<5%) vs EQTL z-scores')) +
        labs(caption = 'Risk Loci by BEST GWAS')
    )
    
    dev.off()
}

create_by_status <- function(tt_sigonly, filename = 'pdf/twas_fdr5perc_by_status.pdf', titleslug = 'FDR') {
    pdf(filename, useDingbats = FALSE, width = 28, height = 14)
    print(ggplot(tt_sigonly, aes(
        y = -log10(TWAS.P),
        x = ifelse(EQTL.status == 'Other', 'Other', 'Risk Locus'),
        fill = BEST.GWAS.P.computed < 5e-08
    )) + geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_point(aes(fill = BEST.GWAS.P.computed < 5e-08), shape = 21, position = position_jitterdodge(jitter.width = 0.2)) +
        facet_grid(region ~ factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
        ) +
        theme_bw(base_size = 30) +
        ggtitle(paste0('TWAS (', titleslug, '<5%) by locus')) +
        xlab('Risk Loci assignment by EQTL SNP') +
        ylim(c(0, max(-log10(tt_sigonly$TWAS.P))))
    )

    print(ggplot(tt_sigonly, aes(
        y = -log10(TWAS.P),
        x = ifelse(BEST.GWAS.status == 'Other', 'Other', 'Risk Locus'),
        fill = BEST.GWAS.P.computed < 5e-08
    )) + geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_point(aes(fill = BEST.GWAS.P.computed < 5e-08), shape = 21, position = position_jitterdodge(jitter.width = 0.2)) +
        facet_grid(region ~ factor(feature, levels = c('gene', 'exon', 'jxn', 'tx'))
        ) +
        theme_bw(base_size = 30) +
        ggtitle(paste0('TWAS (', titleslug, '<5%) by locus')) +
        xlab('Risk Loci assignment by BEST GWAS SNP') +
        ylim(c(0, max(-log10(tt_sigonly$TWAS.P))))
    )
    dev.off()
    
}

check_by_locus <- function(rag, ref) {
    by_loc <- split(rag$SNP, rag$IndexSNP)
    map_dbl(by_loc, ~ sum(.x %in% ref))
}


clean_by_state <- function(x) {
    r <- map_dfr(x, ~ .x)
    r$state <- c(FALSE, TRUE)
    return(r)
}

by_locus <- function(cut, list = FALSE, var = 'TWAS.FDR') {
    by_locus <- map2(
        raggr_clean,
        map(names(raggr_clean), ~ tt$BEST.GWAS.ID[tt[, var, drop = TRUE] < cut & tt$region == .x]),
        check_by_locus
    )
    if(list) return(by_locus)
    clean_by_state(map(by_locus, ~ table(.x > 0)))
}
perc_locus <- function(cut, var = 'TWAS.FDR') {
    x <- by_locus(cut, var = var)
    x[2, 1:2] / colSums(x[, 1:2]) * 100
}

get_matrix <- function(x) {
    matrix(x, ncol = ncol(x), dimnames = attr(x, 'dimnames'))
}

get_venn_info <- function(cut, var = 'TWAS.FDR') {
    map(
        by_locus(cut, list = TRUE, var = var),
        ~ names(which(.x > 0))
    )
}
venn_by_locus <- function(cut, var = 'TWAS.FDR') {
    venn(
        get_venn_info(cut, var = var),
        show.plot = FALSE
    )
}

shared_by_locus <- function(cut, var = 'TWAS.FDR') {
    get_matrix(
        venn_by_locus(cut, var = var)
    )
}

make_pretty_venn <- function(cut, title = '', var = 'TWAS.FDR') {
    info <- get_venn_info(cut, var = var)
    cols <- c('DLPFC' = 'dark orange', 'HIPPO' = 'skyblue3')
    v <- venn.diagram(info, filename = NULL,
        main = title,
        col = 'transparent', fill = rev(cols),
        alpha = 0.5, margin = 0,
        main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
        cat.col = rev(cols))
    grid.newpage()
    grid.draw(v)
}

gene_by_locus <- function(rag, ref) {
    by_loc <- split(rag$gene, rag$IndexSNP)
    map_dbl(by_loc, ~ sum(.x %in% ref))
}

features <- c('gene', 'exon', 'jxn', 'tx')

by_feature <- function(cut, list = FALSE, var = 'TWAS.FDR') {
    g_by_locus <- map(features, function(feature) {
        map2(
            map(raggr_clean, ~ subset(.x, tolower(Type) == feature)),
            map(names(raggr_clean), ~ tt$ID[tt$feature == feature & tt$region == .x & tt[, var, drop = TRUE] < cut]),
            gene_by_locus
        )
    })
    names(g_by_locus) <- features
    if(list) return(g_by_locus)
    map(g_by_locus, function(x) {
        r <- map_dfr(x, ~ table(.x > 0))
        r$state <- c(FALSE, TRUE)
        return(r)
    })
}

clean_tabs <- function(l) {
    map2_dfr(l, names(l), function(x, y) {
        x$feature <- y
        return(x)
    })
}

perc_feature <- function(cut, var = 'TWAS.FDR') {
    y <- clean_tabs(by_feature(cut, var = var))
    clean_tabs(
        map(
            split(y, factor(y$feature, levels = features)), 
            ~ .x[2, 1:2] / colSums(.x[, 1:2]) * 100
    ))
}

get_venn_info_by_feature <- function(cut, var = 'TWAS.FDR') {
    map(
        by_feature(cut, list = TRUE, var = var),
        ~ map(.x, 
            ~ names(which(.x > 0))
        )
    )
}
venn_by_locus_by_feature <- function(cut, var = 'TWAS.FDR') {
    map(get_venn_info_by_feature(cut, var = var), venn, show.plot = FALSE)
}

shared_by_locus_by_feature <- function(cut, var = 'TWAS.FDR') {
    map(venn_by_locus_by_feature(cut, var = var), get_matrix)
}

make_pretty_venn_by_feature <- function(cut, title = '', var = 'TWAS.FDR') {
    info_all <- get_venn_info_by_feature(cut, var = var)
    cols <- c('DLPFC' = 'dark orange', 'HIPPO' = 'skyblue3')
    map2(
        info_all,
        names(info_all),
        function(info, feature) {
        v <- venn.diagram(info, filename = NULL,
            main = paste0(title, ' - ', feature),
            col = 'transparent', fill = rev(cols),
            alpha = 0.5, margin = 0,
            main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
            cat.col = rev(cols))
        grid.newpage()
        grid.draw(v)
        }
    )
}

get_feat <- function(cut, var = 'TWAS.FDR') {
    x <- tt[ tt[, var, drop = TRUE]< cut, ]
    map(
        split(x, factor(x$feature, levels = features)),
        ~ map(split(.x, .x$region), ~.x$ID)
    )
}

shared_by_feature <- function(cut, var = 'TWAS.FDR') {
    map(
        get_feat(cut, var = var),
        ~ get_matrix(venn(.x, show.plot = FALSE))
    )
}

make_pretty_venn_shared_by_feature <- function(cut, title = '', var = 'TWAS.FDR') {
    info_all <- get_feat(cut, var = var)
    cols <- c('DLPFC' = 'dark orange', 'HIPPO' = 'skyblue3')
    map2(
        info_all,
        names(info_all),
        function(info, feature) {
        v <- venn.diagram(info, filename = NULL,
            main = paste0(title, ' - ', feature),
            col = 'transparent', fill = cols,
            alpha = 0.5, margin = 0,
            main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
            cat.col = cols)
        grid.newpage()
        grid.draw(v)
        }
    )
}

get_feat_geneid <- function(cut, var = 'TWAS.FDR') {
    x <- tt[ tt[, var, drop = TRUE]< cut, ]
    map(
        split(x, factor(x$feature, levels = features)),
        ~ map(split(.x, .x$region), ~ .x$geneid[!is.na(.x$geneid)])
    )
}

shared_by_geneid <- function(cut, var = 'TWAS.FDR') {
    map(
        get_feat_geneid(cut, var = var),
        ~ get_matrix(venn(.x, show.plot = FALSE))
    )
}

make_pretty_venn_shared_by_geneid <- function(cut, title = '', var = 'TWAS.FDR') {
    info_all <- get_feat_geneid(cut, var = var)
    cols <- c('DLPFC' = 'dark orange', 'HIPPO' = 'skyblue3')
    map2(
        info_all,
        names(info_all),
        function(info, feature) {
        v <- venn.diagram(info, filename = NULL,
            main = paste0(title, ' - ', feature),
            col = 'transparent', fill = cols,
            alpha = 0.5, margin = 0,
            main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
            cat.col = cols)
        grid.newpage()
        grid.draw(v)
        }
    )
}

get_feat_geneid2 <- function(cut, var = 'TWAS.FDR') {
    x <- tt[ tt[, var, drop = TRUE]< cut, ]
    map(
        split(x, x$region),
        ~ map(split(.x, factor(.x$feature, levels = features)), ~ .x$geneid[!is.na(.x$geneid)])
    )
}

make_pretty_venn_shared_by_geneid2 <- function(cut, title = '', var = 'TWAS.FDR') {
    info_all <- get_feat_geneid2(cut, var = var)
    cols <- brewer.pal('Set1', n = 4)
    map2(
        info_all,
        names(info_all),
        function(info, feature) {
        v <- venn.diagram(info, filename = NULL,
            main = paste0(title, ' - ', feature),
            col = 'transparent', fill = cols,
            alpha = 0.5, margin = 0,
            main.cex = 2, cex = 2, cat.fontcase = 'bold', cat.cex = 2,
            cat.col = cols)
        grid.newpage()
        grid.draw(v)
        }
    )
}
