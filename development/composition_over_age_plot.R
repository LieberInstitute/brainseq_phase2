## Original: https://github.com/LieberInstitute/DNAm_Hippo/blob/master/brainseq_phase2_composition_over_age_plot_01.R

## script for brainseq phase2 methylation results
library(ggplot2)
library(ggrepel)
library(jaffelab)
library('RColorBrewer')
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


##### Cell Type over Development Plot #########
dat <- pd
dat$Pluripotency = dat$NPC + dat$ES
dat$Glial = dat$NeuN_neg
dat = dat[,c('BrNum','Brain.Region','ageGroup', 'Pluripotency','Glial')]
dat = tidyr::gather(dat, key="CellType", value="Proportion", Pluripotency, Glial)
dat$CellType = factor(dat$CellType,levels=c('Pluripotency','Glial') )
dat$Region=dat$`Brain.Region`

df3 = data.frame(ageGroup = as.numeric(c(NA, NA)),
                 Region = c("DLPFC", "Hippo"),
                 Proportion = as.numeric(c(NA, NA)))
				 
cell_type = ggplot(data=dat, aes(x=ageGroup,y=Proportion,fill=Region ))  + 
geom_boxplot(col='black',show.legend=FALSE) + 
facet_wrap(~`CellType`,scales='free',nrow=1)  + 
labs(x='Age Group', y = "Cell Type Proportion",fill='Region') +   
geom_point(data = df3, aes(x = ageGroup, y = Proportion, col = Region), size=8, shape=15) +
  scale_fill_manual(values = brewer.pal(5,'Set1')[1:2], breaks = c("DLPFC","Hippo")) +
  scale_color_manual(values = brewer.pal(5,'Set1')[1:2] , breaks = c("DLPFC","Hippo")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = c(0.375, 0.88),legend.background = element_rect(fill = "white", colour = 'black',linetype='solid'), legend.key = element_blank(),axis.title.x=element_blank()) 

ggsave(cell_type, filename='/dcl01/lieber/ajaffe/Steve/Hippo_meQTL/brainseq_phase2/plots/bothRegions_estimated_cellType_proportions_over_lifespan_n347_for_BoG_may2018_poster.pdf',height=8,width=12)
ggsave(cell_type, filename='/dcl01/lieber/ajaffe/Steve/Hippo_meQTL/brainseq_phase2/plots/bothRegions_estimated_cellType_proportions_over_lifespan_n347_for_BoG_may2018_poster.png',height=8,width=12)
