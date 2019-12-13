###
library(LIBDpheno)
pd = read.csv("BrainSeq_Phase2_phenotype_data_small_n900.csv")
br = unique(pd$BrNum)

## add most recent phenotype data from lims
pheno = toxicant[[1]]
pheno = pheno[match(br, paste0("Br", pheno$brnumerical)),]
pheno$BrNum = paste0("Br", pheno$brnumerical)

## just smoking info
vars = c("BrNum", "smoking", "cotinine", "nicotine", "nicotine_comments")
pheno = pheno[,vars]
colnames(pheno)[2] = "smoking_hx"

## based on steve's previous code

################################################################
###### Nicotine Brain Numbers
pheno$nicotine_brain_number <- NA
pheno$nicotine_brain_number <- ifelse( grepl('ng/g|pg/mg', pheno$nicotine_comments), pheno$nicotine_comments, NA)
pheno$nicotine_brain_number <- gsub(';.*', '', pheno$nicotine_brain_number, ignore.case = TRUE) #removing all text after ;
pheno$nicotine_brain_number <- gsub("[^0-9|.]", "", pheno$nicotine_brain_number)
pheno$nicotine_brain_number <- as.numeric(pheno$nicotine_brain_number)
#pheno$nicotine_brain_number[which(!pheno$nicotine & !df$cotinine & df$primarydx=='Control'& grepl('CB', df$nicotine_comments))] <- 0
######Cotinine Brain Number
pheno$cotinine_brain_number <- NA
pheno$cotinine_brain_number <- ifelse( grepl('ng/g|pg/mg', pheno$nicotine_comments), pheno$nicotine_comments, NA)
pheno$cotinine_brain_number <- gsub('.*cotinine', '', pheno$cotinine_brain_number, ignore.case = TRUE) #removing all text after ;
pheno$cotinine_brain_number <- gsub('nicotine.*', '', pheno$cotinine_brain_number, ignore.case = TRUE) #removing all text after ;
pheno$cotinine_brain_number <- gsub("[^0-9|.]", "", pheno$cotinine_brain_number)
pheno$cotinine_brain_number <- as.numeric(pheno$cotinine_brain_number)
#pheno$nicotine_brain_number[which(!pheno$nicotine & !df$cotinine & df$primarydx=='Control'& grepl('CB', df$nicotine_comments))] <- 0
###### Nicotine Blood Numbers
pheno$nicotine_blood_number <- NA
pheno$nicotine_blood_number <- ifelse(grepl('ng/mL', pheno$nicotine_comments), pheno$nicotine_comments, NA)
pheno$nicotine_blood_number <- gsub('cotinine.*', '', pheno$nicotine_blood_number, ignore.case = TRUE) #removing all text after ng/mL
pheno$nicotine_blood_number <- gsub('.*Nicotine', '', pheno$nicotine_blood_number, ignore.case = TRUE) #removing all text before Nicotine
pheno$nicotine_blood_number <- gsub("[^0-9|.]", "", pheno$nicotine_blood_number)
pheno$nicotine_blood_number <- as.numeric(pheno$nicotine_blood_number)
###### Cotinine Blood Concentrations
pheno$cotinine_blood_number <- NA
pheno$cotinine_blood_number <- ifelse(grepl('ng/mL', pheno$nicotine_comments), pheno$nicotine_comments, NA)
pheno$cotinine_blood_number <- gsub('.*cotinine', '', pheno$cotinine_blood_number, ignore.case = TRUE) #removing all text before ; (usually nicotine said first)
pheno$cotinine_blood_number <- gsub("[^0-9|.]", "", pheno$cotinine_blood_number)
pheno$cotinine_blood_number <- as.numeric(pheno$cotinine_blood_number)
################################################################

#pheno = pheno[pheno$Age <0,]
pheno$modelGroup <- pheno$cotinine| pheno$BrNum %in% c('Br1813','Br1826','Br2047', "Br2049")
pheno$modelGroup <- as.factor(pheno$modelGroup)
pheno$modelGroup <- plyr::revalue(pheno$modelGroup, c("TRUE" = "Smoker", "FALSE" = "Non-Smoker"))
pheno$modelGroup[which(pheno$nicotine)] <- "Smoker"

pheno <- pheno[!is.na(pheno$modelGroup),]
pheno <- pheno[pheno$agedeath<0|pheno$agedeath>16,]
pheno$AgeGroup = factor(ifelse(pheno$agedeath<0, "Fetal","Adult"))
pheno<- pheno[-which(pheno$AgeGroup=="Adult"&pheno$smoking & pheno$modelGroup =="Non-Smoker"),]


write.csv(pheno, file = "/dcl01/lieber/ajaffe/lab/Nicotine/smoking_info_bs2.csv",
	row.names=FALSE)