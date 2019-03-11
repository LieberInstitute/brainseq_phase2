SCZD and protein coding gene check
==================================

BrainSeq Phase II RNA-seq data was prepared with a RiboZero library preparation protocol that allows preservation of non-coding RNA. We found that features passing our expression cutoffs are enriched for protein coding genes compared to the other 44 annotation categories using Gencode v25. That's because 19,950 (34.4%) out of 58,037 genes are protein coding yet 14,653 (59.4%) out of the 24,652 genes passing our expression cutoff are protein coding; that is 14,653 (73.4%) out of 19,950 protein coding genes pass the expression cutoff. In contrast lincRNAs decrease from 13% to 8% after the expression cuts with 1,965 (26.1%) out of 7,539 lincRNAs passing the expression cutoff (see [check_protein_coding.R](https://github.com/LieberInstitute/brainseq_phase2/blob/master/check_protein_coding/check_protein_coding.R#L99-L144) for more details). Below is the protein coding enrichment information at the gene feature level.

![Expressed vs PC](protein_coding_checks_gene1.png)

Among genes used for our SCZD case-control differential expression analysis (FDR <5%) in each brain region, we observe a significant enrichment (p-value <5%) for protein coding genes in the differentially expressed genes. 

![SCZD DE vs PC](protein_coding_checks_gene2.png)

At [check_protein_coding.R](https://github.com/LieberInstitute/brainseq_phase2/blob/master/check_protein_coding/check_protein_coding.R#L310-L354) we included the list of the 33 and 7 SCZD differentially expressed genes that are __not__ protein coding in DLPFC and HIPPO respectively, including 11 lincRNAs in DLPFC.

Note that `gene_type` is a column we provided in [supplementary table 2](https://github.com/LieberInstitute/brainseq_phase2/tree/master/supp_tabs).