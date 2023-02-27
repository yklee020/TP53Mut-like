#2022.11.09
#GEP analysis using updated data: 401 BEAT AML and 178 TCGA

source("C:/Users/yklee/Desktop/Sachs Lab/R_project/TP53/TP53GEP_functions.r")
packages<-c('edgeR','limma','dplyr','tidyr','sva','tidyverse',
            'ggplot2','survival','survminer','HGNChelper','stringr','ComplexHeatmap')
for(p in packages){
  library(p,character.only = T)
}


#new Analyses to understand TP53Mut-like WT 
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version")

BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')

#below GEP data has genes were transformed into same reference
BEAT_RNA_counts<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_counts_mtmtlwt.rds")
TCGA_AML_counts<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/TCGA_AML_counts_mtmtlwt.rds")
# BEAT_TCGA_merge.corrected<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_TCGA_OL_data/BEAT461_TCGA_combat_count.rds")


keep.beat<-rowSums(BEAT_RNA_counts>0)>=5
keep.tcga<-rowSums(TCGA_AML_counts>0)>=5
#keep.comb<-rowSums(BEAT_TCGA_merge.corrected>0)>=5
y.beat<-DGEList(counts=BEAT_RNA_counts)
y.beat<-y.beat[keep.beat,]
y.beat<-calcNormFactors(y.beat)
y.beat.cpm<-as.data.frame(cpm(y.beat,log=T,normalized.lib.sizes = TRUE))
y.tcga<-DGEList(counts=TCGA_AML_counts)
colnames(y.tcga)<-gsub('.03A.*|.03B.*','',colnames(y.tcga))
y.tcga<-y.tcga[keep.tcga,]
y.tcga<-calcNormFactors(y.tcga)
y.tcga.cpm<-as.data.frame(cpm(y.tcga,log=T,normalized.lib.sizes = TRUE))

OL_genes<-intersect(rownames(BEAT_RNA_counts),rownames(TCGA_AML_counts))
BEAT_RNA_counts.OL<-BEAT_RNA_counts[OL_genes,]
TCGA_AML_counts.OL<-TCGA_AML_counts[OL_genes,]
BEAT_TCGA_counts.OL<-cbind(BEAT_RNA_counts.OL,TCGA_AML_counts.OL)
keep<-rowSums(BEAT_TCGA_counts.OL>0)>=5
BEAT_TCGA_counts.OL<-BEAT_TCGA_counts.OL[keep,]
beat_tcga_group<-c(rep('BEAT',dim(BEAT_RNA_counts.OL)[2]),
                   rep('TCGA',dim(TCGA_AML_counts.OL)[2]))
BEAT_TCGA_merge.corrected<-ComBat_seq(counts = as.matrix(BEAT_TCGA_counts.OL), batch = beat_tcga_group)
#saveRDS(BEAT_TCGA_merge.corrected,"C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_TCGA_OL_data/BEAT403_TCGA_combat_count.rds")
BEAT_TCGA_merge.corrected<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_TCGA_OL_data/BEAT461_TCGA_combat_count.rds")
y.comb<-DGEList(counts=BEAT_TCGA_merge.corrected)
y.comb<-calcNormFactors(y.comb)
y.comb.logcpm<-as.data.frame(log2(cpm(y.comb,log=FALSE,normalized.lib.sizes = TRUE)+1))
#rm(BEAT_TCGA_merge.corrected)

#DO PCA
###############
#BEAT_TCGA_merge.corrected<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_TCGA_OL_data/BEAT461_TCGA_combat_count.rds")
beat_tcga_group
BEAT_MUT_id<-BEAT_clinical_final$LabId[BEAT_clinical_final$cls.GLM=='TP53_MUT']
BEAT_MUT_like_id<-BEAT_clinical_final$LabId[BEAT_clinical_final$cls.GLM=='TP53_MUTlike']#actually just WT_MTlike, no unknown added
BEAT_WT_id<-BEAT_clinical_final$LabId[BEAT_clinical_final$cls.GLM=='TP53_WT']#actually just WT_MTlike, no unknown added

TCGA_MUT_id<-TCGA_clinical_final$LabId.y[TCGA_clinical_final$cls.GLM.y=='TP53_MUT']
TCGA_MUT_like_id<-TCGA_clinical_final$LabId.y[TCGA_clinical_final$cls.GLM.y=='TP53MUT_like']
TCGA_WT_id<-TCGA_clinical_final$LabId.y[TCGA_clinical_final$cls.GLM.y=='TP53_WT']
beat_tcga_tp53<-c(gsub('.03A.*|.03B.*','',colnames(y.comb)))
beat_tcga_tp53<-ifelse(beat_tcga_tp53 %in% c(BEAT_MUT_id,TCGA_MUT_id),'TP53MUT','TP53WT')

pca.corrected <- prcomp(t(y.comb.logcpm),scale. = F)
pcaData.corrected <- as.data.frame(pca.corrected$x)
pcaData.corrected$sampleID <- rownames(pcaData.corrected)
pcaData.corrected$group<-beat_tcga_group
pcaData.corrected$TP53<-beat_tcga_tp53

library(ggplot2)
p.comb.correct01<-ggplot(pcaData.corrected, aes(x = PC1, y = PC2,color=group))+geom_point()+
  #labs(caption = "pca:logcpm_allgenes, group=TP53Mut")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),
        axis.title = element_text(size=14),title = element_text(size=14),
        legend.text=element_text(size=14))+  #theme(legend.text=element_text(size=18))
   guides(colour = guide_legend(override.aes = list(size=3)))

p.comb.correct02<-ggplot(pcaData.corrected, aes(x = PC1, y = PC2,color=TP53))+geom_point()+
  #labs(caption = "pca:logcpm_allgenes, group=TP53Mut")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),
        axis.title = element_text(size=14),title = element_text(size=14),
        legend.text=element_text(size=14))+  #theme(legend.text=element_text(size=18))
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_color_manual(values=c('red','grey'))
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/GEP/Clustering/")
jpeg('PCA_BEAT_TCGA.jpeg',width=5.2,height=3.8,units='in',res=300)
p.comb.correct01
dev.off()

jpeg('PCA_BEAT_TCGA_tp53mut.jpeg',width=5.2,height=3.8,units='in',res=300)
p.comb.correct02
dev.off()


###############
#Do HC heat map

library(ComplexHeatmap)
library(matrixStats)
#subset using highly variable genes

var5K.beat<-head(order(rowVars(as.matrix(y.beat.cpm)),decreasing = T),5000)
var8K.beat<-head(order(rowVars(as.matrix(y.beat.cpm)),decreasing = T),8000)

var5K.tcga<-head(order(rowVars(as.matrix(y.tcga.cpm)),decreasing = T),5000)
var8K.tcga<-head(order(rowVars(as.matrix(y.tcga.cpm)),decreasing = T),8000)

# meangenes<-rowMeans(as.matrix(RNAseq_logcpm))
# sdgenes<-rowSds(as.matrix(RNAseq_logcpm))
# rowstat<-cbind(meangenes,sdgenes,vargenes)

Zmat_5K.beat<-t(scale(t(y.beat.cpm[var5K.beat,])))
Zmat_8K.beat<-t(scale(t(y.beat.cpm[var8K.beat,])))

Zmat_5K.tcga<-t(scale(t(y.tcga.cpm[var5K.tcga,])))
Zmat_8K.tcga<-t(scale(t(y.tcga.cpm[var8K.tcga,])))

#group_anno<-data.frame('TP53_group'=TP53_group,'ZS_anno_group'=ZS_anno_group,'LSC17_score'=LSC17_score)
#rownames(group_anno)<-colnames(RNAseq_logcpm)
group_anno.beat<-data.frame('LabId'=colnames(y.beat.cpm))
group_anno.beat<-right_join(group_anno.beat,BEAT_clinical_final,by='LabId') 
rownames(group_anno.beat)<-colnames(y.beat.cpm)

group_anno.beat.trim<-group_anno.beat[,c(110,111)]
group_anno.beat.trim$TP53_mut_stat<-ifelse(group_anno.beat.trim$TP53_mut_stat==1,'TP53_MUT','TP53_WT')
group_anno.tcga<-data.frame('patient_RNA_id'=gsub('.03A.*|.03B.*','',colnames(y.tcga.cpm)))
group_anno.tcga<-right_join(group_anno.tcga,TCGA_clinical_final,by='patient_RNA_id') 
rownames(group_anno.tcga)<-colnames(y.tcga.cpm)
group_anno.tcga.trim<-group_anno.tcga[,c(3,2102)]

library(circlize)
group_anno_h.beat<-HeatmapAnnotation(df=group_anno.beat.trim,
                                     col = list(TP53_mut_stat = c("TP53_MUT" = "red", "TP53_WT" = "blue", "unknown" = "white")))
group_anno_h.tcga<-HeatmapAnnotation(df=group_anno.tcga.trim,
                                     col = list(TP53_mut_stat = c("TP53_MUT" = "red", "TP53_WT" = "blue")))
HM_5K_pear.avg<-Heatmap(Zmat_5K.beat[1:10,], name = "Zmat_5K", clustering_distance_columns = "pearson",
                        clustering_method_columns = "average",#complete,ward.D2
                        column_title = "pre-defined distance method (1 - pearson),complete",
                        show_row_names = F,show_column_names = F,top_annotation = group_anno_h.beat)


HM_5K_pear.avg<-Heatmap(Zmat_5K.beat, name = "Zmat_5K", clustering_distance_columns = "pearson",
                        clustering_method_columns = "average",#complete,ward.D2
                        column_title = "pre-defined distance method (1 - pearson),complete",
                        show_row_names = F,show_column_names = F,top_annotation = group_anno_h.beat)
HM_8K_pear.avg<-Heatmap(Zmat_8K.beat, name = "Zmat_5K", clustering_distance_columns = "pearson",
                        clustering_method_columns = "average",#complete,ward.D2
                        column_title = "pre-defined distance method (1 - pearson),complete",
                        show_row_names = F,show_column_names = F,top_annotation = group_anno_h.beat)

HM_tcga5K_pear.avg<-Heatmap(Zmat_5K.tcga, name = "Zmat_5K", clustering_distance_columns = "pearson",
                            clustering_method_columns = "average",#complete,ward.D2
                            column_title = "pre-defined distance method (1 - pearson),complete",
                            show_row_names = F,show_column_names = F,top_annotation = group_anno_h.tcga)
HM_tcga8K_pear.avg<-Heatmap(Zmat_8K.tcga, name = "Zmat_5K", clustering_distance_columns = "pearson",
                            clustering_method_columns = "average",#complete,ward.D2
                            column_title = "pre-defined distance method (1 - pearson),complete",
                            show_row_names = F,show_column_names = F,top_annotation = group_anno_h.tcga)

pdf("Analysis/GEP/BEAT_AML_HC_pearson_var5_8k.pdf")
print(HM_5K_pear.avg)
print(HM_8K_pear.avg)
dev.off()
pdf("Analysis/GEP/TCGA_AML_HC_pearson_var5_8k.pdf")
print(HM_tcga5K_pear.avg)
print(HM_tcga8K_pear.avg)
dev.off()

###############
#do PCA for TP53mut-like

#new Analyses to understand TP53Mut-like WT 
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version")
BEAT_RNA_counts<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_counts_mtmtlwt.rds")
TCGA_AML_counts<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/TCGA_AML_counts_mtmtlwt.rds")

BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')


#saveRDS(BEAT_TCGA_merge.corrected,"C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_TCGA_OL_data/BEAT403_TCGA_combat_count.rds")
BEAT_TCGA_merge.corrected<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_TCGA_OL_data/BEAT461_TCGA_combat_count.rds")
y.comb<-DGEList(counts=BEAT_TCGA_merge.corrected)
y.comb<-calcNormFactors(y.comb)
y.comb.logcpm<-as.data.frame(log2(cpm(y.comb,log=FALSE,normalized.lib.sizes = TRUE)+1))


beat_tcga_group
BEAT_MUT_id<-BEAT_clinical_final$LabId[BEAT_clinical_final$cls.GLM=='TP53_MUT']
BEAT_MUT_like_id<-BEAT_clinical_final$LabId[BEAT_clinical_final$cls.GLM=='TP53_MUTlike']#actually just WT_MTlike, no unknown added
BEAT_WT_id<-BEAT_clinical_final$LabId[BEAT_clinical_final$cls.GLM=='TP53_WT']#actually just WT_MTlike, no unknown added

TCGA_MUT_id<-TCGA_clinical_final$LabId.y[TCGA_clinical_final$cls.GLM.y=='TP53_MUT']
TCGA_MUT_like_id<-TCGA_clinical_final$LabId.y[TCGA_clinical_final$cls.GLM.y=='TP53MUT_like']
TCGA_WT_id<-TCGA_clinical_final$LabId.y[TCGA_clinical_final$cls.GLM.y=='TP53_WT']
# beat_tcga_tp53<-c(gsub('.03A.*|.03B.*','',colnames(y.comb)))
beat_tcga_tp53<-colnames(y.comb)
beat_tcga_tp53<-ifelse(beat_tcga_tp53 %in% c(BEAT_MUT_id,TCGA_MUT_id),'TP53Mut','TP53Wt')
#beat_tcga_tp53mtl<-c(gsub('.03A.*|.03B.*','',colnames(y.comb)))
beat_tcga_tp53mtl<-colnames(y.comb)
beat_tcga_tp53mtl<-ifelse(beat_tcga_tp53mtl %in% c(BEAT_MUT_like_id,TCGA_MUT_like_id),'TP53Mut-like','TP53Wt')
beat_tcga_group<-c(rep('BEAT',461),
                   rep('TCGA',178))

pca.corrected <- prcomp(t(y.comb.logcpm),scale. = F)
pcaData.corrected <- as.data.frame(pca.corrected$x)
pcaData.corrected$sampleID <- rownames(pcaData.corrected)
pcaData.corrected$group<-beat_tcga_group
pcaData.corrected$TP53mt<-beat_tcga_tp53
pcaData.corrected$TP53mtl<-beat_tcga_tp53mtl

library(ggplot2)
p.comb.correct01<-ggplot(pcaData.corrected, aes(x = PC1, y = PC2,color=group))+geom_point()+
  #labs(caption = "pca:logcpm_allgenes, group=TP53Mut")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),
        axis.title = element_text(size=14),title = element_text(size=14),
        legend.text=element_text(size=14))+  #theme(legend.text=element_text(size=18))
  guides(colour = guide_legend(override.aes = list(size=3)))

p.comb.correct02<-ggplot(pcaData.corrected, aes(x = PC1, y = PC2,color=TP53mt))+geom_point()+
  #labs(caption = "pca:logcpm_allgenes, group=TP53Mut")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),
        axis.title = element_text(size=14),title = element_text(size=14),
        legend.text=element_text(size=14))+  #theme(legend.text=element_text(size=18))
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_color_manual(values=c('red','grey'))
p.comb.correct03<-ggplot(pcaData.corrected, aes(x = PC1, y = PC2,color=TP53mtl))+geom_point()+
  #labs(caption = "pca:logcpm_allgenes, group=TP53Mut")+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),
        axis.title = element_text(size=14),title = element_text(size=14),
        legend.text=element_text(size=14))+  #theme(legend.text=element_text(size=18))
  guides(colour = guide_legend(override.aes = list(size=3)))+
  scale_color_manual(values=c('purple','grey'))


setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/GEP/Clustering/")
jpeg('PCA_BEAT_TCGA.jpeg',width=5.2,height=3.8,units='in',res=300)
p.comb.correct01
dev.off()

jpeg('PCA_BEAT_TCGA_tp53mut.jpeg',width=5.2,height=3.8,units='in',res=300)
p.comb.correct02
dev.off()

jpeg('PCA_BEAT_TCGA_tp53mutlike.jpeg',width=5.2,height=3.8,units='in',res=300)
p.comb.correct03
dev.off()




###############
#DE gene first and GSEA


setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version")

BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')

#below GEP data has genes were transformed into same reference
BEAT_RNA_counts<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_counts_mtmtlwt.rds")
TCGA_AML_counts<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/TCGA_AML_counts_mtmtlwt.rds")
#BEAT_TCGA_merge.corrected<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_TCGA_OL_data/BEAT461_TCGA_combat_count.rds")
msigdb_filt<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/GSEA/Compiled data/Gene set db/MSIGDB_ZS_gs_filtered_long_complete_gene_conv2022.rds')



BEAT_MUT_id<-BEAT_clinical_final$LabId[BEAT_clinical_final$cls.GLM=='TP53_MUT']
BEAT_MUT_like_id<-BEAT_clinical_final$LabId[BEAT_clinical_final$cls.GLM=='TP53_MUTlike']#actually just WT_MTlike, no unknown added
BEAT_WT_id<-BEAT_clinical_final$LabId[BEAT_clinical_final$cls.GLM=='TP53_WT']#actually just WT_MTlike, no unknown added

TCGA_MUT_id<-TCGA_clinical_final$patient_RNA_id[TCGA_clinical_final$cls.GLM.y=='TP53_MUT']
TCGA_MUT_like_id<-TCGA_clinical_final$patient_RNA_id[TCGA_clinical_final$cls.GLM.y=='TP53MUT_like']
TCGA_WT_id<-TCGA_clinical_final$patient_RNA_id[TCGA_clinical_final$cls.GLM.y=='TP53_WT']




# BEAT_TCGA_count<-cbind(BEAT_RNA_counts.461.clinic[OL_genes, ],
#                        TCGA_AML_counts[OL_genes, ])

keep.beat<-rowSums(BEAT_RNA_counts>0)>=5
keep.tcga<-rowSums(TCGA_AML_counts>0)>=5
#keep.comb<-rowSums(BEAT_TCGA_merge.corrected>0)>=5
y.beat<-DGEList(counts=BEAT_RNA_counts)
y.beat<-y.beat[keep.beat,]
y.beat<-calcNormFactors(y.beat)
y.tcga<-DGEList(counts=TCGA_AML_counts)
y.tcga<-y.tcga[keep.tcga,]
y.tcga<-calcNormFactors(y.tcga)
colnames(y.tcga)<-gsub('.03A.*|.03B.*','',colnames(y.tcga))

all.equal(colnames(y.beat),BEAT_clinical_final$LabId)#TRUE
all.equal(colnames(y.tcga),TCGA_clinical_final$patient_RNA_id)
# TP53WT_ix<-which(BEAT_clinical_final$TP53_mut_stat != 'unknown' & BEAT_clinical_final$TP53_mut_stat != 'TP53_MUT')
# TP53WT_ix.tcga<-which(TCGA_clinical_final$TP53_mut_stat != 'TP53_MUT')
TP53WT_id<-which(BEAT_clinical_final$TP53_mut_stat != 'unknown' & BEAT_clinical_final$TP53_mut_stat != 'TP53_MUT')


#need to have 2 diferent data: MT vs. WT and MTlike vs. WTlike
y.beat.mtwt<-y.beat[,c(BEAT_MUT_id,BEAT_WT_id)]
y.beat.mtlwt<-y.beat[,c(BEAT_MUT_like_id,BEAT_WT_id)]

#y.tcga
y.tcga.mtwt<-y.tcga[,c(TCGA_MUT_id,TCGA_WT_id)]
y.tcga.mtlwt<-y.tcga[,c(TCGA_MUT_like_id,TCGA_WT_id)]

#make group annotation
BEAT_cls.mtwt<-rep('WT',dim(y.beat.mtwt)[2]);BEAT_cls.mtwt[colnames(y.beat.mtwt)%in%BEAT_MUT_id]<-'TP53_MUT'
BEAT_cls.mtlwt<-rep('WT',dim(y.beat.mtlwt)[2]);BEAT_cls.mtlwt[colnames(y.beat.mtlwt)%in%BEAT_MUT_like_id]<-'TP53_MUTlike'
TCGA_cls.mtwt<-rep('WT',dim(y.tcga.mtwt)[2]);TCGA_cls.mtwt[gsub('.03A.*|.03B.*','',colnames(y.tcga.mtwt))%in%TCGA_MUT_id]<-'TP53_MUT'
TCGA_cls.mtlwt<-rep('WT',dim(y.tcga.mtlwt)[2]);TCGA_cls.mtlwt[gsub('.03A.*|.03B.*','',colnames(y.tcga.mtlwt))%in%TCGA_MUT_like_id]<-'TP53_MUTlike'

#data to use
y.beat.mtwt
y.beat.mtlwt
y.tcga.mtwt
y.tcga.mtlwt


beat.mtwt.group<-factor(ifelse(BEAT_cls.mtwt=='TP53_MUT',1,0))
beat.mtlwt.group<-factor(ifelse(BEAT_cls.mtlwt=='TP53_MUTlike',1,0))

tcga.mtwt.group<-factor(ifelse(TCGA_cls.mtwt=='TP53_MUT',1,0))
tcga.mtlwt.group<-factor(ifelse(TCGA_cls.mtlwt=='TP53_MUTlike',1,0))

design.beat.mtwt<-model.matrix(~beat.mtwt.group)
design.beat.mtlwt<-model.matrix(~beat.mtlwt.group)
design.tcga.mtwt<-model.matrix(~tcga.mtwt.group)
design.tcga.mtlwt<-model.matrix(~tcga.mtlwt.group)

y.beat.mtwt<-estimateDisp(y.beat.mtwt,design = design.beat.mtwt)
y.beat.mtlwt<-estimateDisp(y.beat.mtlwt,design = design.beat.mtlwt)
y.tcga.mtwt<-estimateDisp(y.tcga.mtwt,design = design.tcga.mtwt)
y.tcga.mtlwt<-estimateDisp(y.tcga.mtlwt,design = design.tcga.mtlwt)

# win.graph()
# plotBCV(y.beat.mtwt.n)
# plotBCV(y.beat.mtlwt)
# plotBCV(y.tcga.mtwt)
# plotBCV(y.tcga.mtlwt)
# plotBCV(y.comb.mtwt)
# plotBCV(y.comb.mtlwt)

fit.beat.mtwt<-glmQLFit(y.beat.mtwt,design = design.beat.mtwt,robust=T)
fit.beat.mtlwt<-glmQLFit(y.beat.mtlwt,design = design.beat.mtlwt,robust=T)
fit.tcga.mtwt<-glmQLFit(y.tcga.mtwt,design = design.tcga.mtwt,robust=T)
fit.tcga.mtlwt<-glmQLFit(y.tcga.mtlwt,design = design.tcga.mtlwt,robust=T)

qlf.beat.mtwt.qlf <- glmQLFTest(fit.beat.mtwt,coef=2)
qlf.beat.mtwt.treat <- glmTreat(fit.beat.mtwt,coef=2,lfc=1)
qlf.beat.mtlwt.qlf <- glmQLFTest(fit.beat.mtlwt,coef=2)
qlf.beat.mtlwt.treat <- glmTreat(fit.beat.mtlwt,coef=2,lfc=1)

qlf.tcga.mtwt.qlf <- glmQLFTest(fit.tcga.mtwt,coef=2)
qlf.tcga.mtwt.treat <- glmTreat(fit.tcga.mtwt,coef=2,lfc=1)
qlf.tcga.mtlwt.qlf <- glmQLFTest(fit.tcga.mtlwt,coef=2)
qlf.tcga.mtlwt.treat <- glmTreat(fit.tcga.mtlwt,coef=2,lfc=1)

summary(decideTestsDGE(qlf.beat.mtwt.qlf))
summary(decideTestsDGE(qlf.beat.mtwt.treat))
summary(decideTestsDGE(qlf.beat.mtlwt.qlf))
summary(decideTestsDGE(qlf.beat.mtlwt.treat))

summary(decideTestsDGE(qlf.tcga.mtwt.qlf))
summary(decideTestsDGE(qlf.tcga.mtwt.treat))
summary(decideTestsDGE(qlf.tcga.mtlwt.qlf))
summary(decideTestsDGE(qlf.tcga.mtlwt.treat))

DEG.beat.mtwt.qlf<-(topTags(qlf.beat.mtwt.qlf,Inf)$table)
DEG.beat.mtwt.treat<-(topTags(qlf.beat.mtwt.treat,Inf)$table)
DEG.beat.mtlwt.qlf<-(topTags(qlf.beat.mtlwt.qlf,Inf)$table)
DEG.beat.mtlwt.treat<-(topTags(qlf.beat.mtlwt.treat,Inf)$table)

DEG.tcga.mtwt.qlf<-(topTags(qlf.tcga.mtwt.qlf,Inf)$table)
DEG.tcga.mtwt.treat<-(topTags(qlf.tcga.mtwt.treat,Inf)$table)
DEG.tcga.mtlwt.qlf<-(topTags(qlf.tcga.mtlwt.qlf,Inf)$table)
DEG.tcga.mtlwt.treat<-(topTags(qlf.tcga.mtlwt.treat,Inf)$table)
# 
write.csv(DEG.beat.mtwt.qlf,'Analysis/GEP/DEG_GSEA/DEG/BEAT_MTWT_DEG_qlf_fin.csv')
write.csv(DEG.beat.mtwt.treat,'Analysis/GEP/DEG_GSEA/DEG/BEAT_MTWT_DEG_treat_fin.csv')
write.csv(DEG.beat.mtlwt.qlf,'Analysis/GEP/DEG_GSEA/DEG/BEAT_MTlikeWT_DEG_qlf_fin.csv')
write.csv(DEG.beat.mtlwt.treat,'Analysis/GEP/DEG_GSEA/DEG/BEAT_MTlikeWT_DEG_treat_fin.csv')

write.csv(DEG.tcga.mtwt.qlf,'Analysis/GEP/DEG_GSEA/DEG/TCGA_MTWT_DEG_qlf_fin.csv')
write.csv(DEG.tcga.mtwt.treat,'Analysis/GEP/DEG_GSEA/DEG/TCGA_MTWT_DEG_treat_fin.csv')
write.csv(DEG.tcga.mtlwt.qlf,'Analysis/GEP/DEG_GSEA/DEG/TCGA_MTlikeWT_DEG_qlf_fin.csv')
write.csv(DEG.tcga.mtlwt.treat,'Analysis/GEP/DEG_GSEA/DEG//TCGA_MTlikeWT_DEG_treat_fin.csv')
shared_DEG_updn_sig<-function(DEG_01.df,DEG_02.df,pval_cutoff=0.05){
  #pval_cutoff
  DEG01<-DEG_01.df
  DEG02<-DEG_02.df
  OL_DEG<-intersect(rownames(DEG01)[DEG01$FDR<pval_cutoff],
                    rownames(DEG02)[DEG02$FDR<pval_cutoff])
  shared.all<-data.frame('DEG'=OL_DEG,'DEG01_FC'=rep(NA,length(OL_DEG)),'DEG02_FC'=rep(NA,length(OL_DEG)),
                         'DEG01_FDR'=rep(NA,length(OL_DEG)),'DEG02_FDR'=rep(NA,length(OL_DEG)))
  for(i in 1:length(shared.all$DEG)){
    DEG01_ix<-which(rownames(DEG01) %in% shared.all$DEG[i])
    DEG02_ix<-which(rownames(DEG02) %in% shared.all$DEG[i])
    shared.all$DEG01_FC[i]<-DEG01$logFC[DEG01_ix]
    shared.all$DEG02_FC[i]<-DEG02$logFC[DEG02_ix]
    shared.all$DEG01_FDR[i]<-DEG01$FDR[DEG01_ix]
    shared.all$DEG02_FDR[i]<-DEG02$FDR[DEG02_ix]
  }
  shared.all.clean<-shared.all[c(which(shared.all$DEG01_FC>0 & shared.all$DEG02_FC>0),
                                 which(shared.all$DEG01_FC<0 & shared.all$DEG02_FC<0)),]
  return(shared.all.clean)
}
OL_genes<-intersect(rownames(y.beat),rownames(y.tcga))

DEG.beat.mtwt.qlf.OL<-DEG.beat.mtwt.qlf[match(OL_genes,rownames(DEG.beat.mtwt.qlf)),]%>%na.omit()
DEG.tcga.mtwt.qlf.OL<-DEG.tcga.mtwt.qlf[match(OL_genes,rownames(DEG.tcga.mtwt.qlf)),]%>%na.omit()
DEG.beat.mtwt.qlf.OL$FDR<-p.adjust(DEG.beat.mtwt.qlf.OL$PValue,method='BH')
DEG.tcga.mtwt.qlf.OL$FDR<-p.adjust(DEG.tcga.mtwt.qlf.OL$PValue,method='BH')

DEG.beat.mtlwt.qlf.OL<-DEG.beat.mtlwt.qlf[match(OL_genes,rownames(DEG.beat.mtlwt.qlf)),]%>%na.omit()
DEG.tcga.mtlwt.qlf.OL<-DEG.tcga.mtlwt.qlf[match(OL_genes,rownames(DEG.tcga.mtlwt.qlf)),]%>%na.omit()
DEG.beat.mtlwt.qlf.OL$FDR<-p.adjust(DEG.beat.mtlwt.qlf.OL$PValue,method='BH')
DEG.tcga.mtlwt.qlf.OL$FDR<-p.adjust(DEG.tcga.mtlwt.qlf.OL$PValue,method='BH')


DEG.BEAT.TCGA.mtwt.shared<-shared_DEG_updn_sig(DEG.beat.mtwt.qlf.OL,DEG.tcga.mtwt.qlf.OL,0.1) 
DEG.BEAT.TCGA.mtlwt.shared<-shared_DEG_updn_sig(DEG.beat.mtlwt.qlf.OL,DEG.tcga.mtlwt.qlf.OL,0.1)

DEG.BEAT.TCGA.mtwt.shared.filt<-DEG.BEAT.TCGA.mtwt.shared%>% filter(DEG01_FDR<0.05&DEG02_FDR<0.05)
DEG.BEAT.TCGA.mtlwt.shared.filt<-DEG.BEAT.TCGA.mtlwt.shared%>% filter(DEG01_FDR<0.05&DEG02_FDR<0.05)
# write.csv(DEG.BEAT.TCGA.mtwt.shared,'Analysis/GEP/DEG_GSEA/DEG/BEAT_TCGA_OL_MTWT_DEG_qlf_fin.csv')
# write.csv(DEG.BEAT.TCGA.mtlwt.shared,'Analysis/GEP/DEG_GSEA/DEG/BEAT_TCGA_OL_MTlWT_DEG_qlf_fin.csv')
# write.csv(DEG.BEAT.TCGA.mtwt.shared.filt,'Analysis/GEP/DEG_GSEA/DEG/BEAT_TCGA_OL_MTWT_DEG_qlf_fin_010.csv')
# write.csv(DEG.BEAT.TCGA.mtlwt.shared.filt,'Analysis/GEP/DEG_GSEA/DEG/BEAT_TCGA_OL_MTlWT_DEG_qlf_fin_010.csv')

# write.csv(DEG.BEAT.TCGA.mtwt.shared.filt,'DEG_GSEA/DEG/Reanalysis_aug02/BEAT_TCGA_OL_MTWT_DEG_qlf_OL_re_010.csv')
# write.csv(DEG.BEAT.TCGA.mtlwt.shared.filt,'DEG_GSEA/DEG/Reanalysis_aug02/BEAT_TCGA_OL_MTlWT_DEG_qlf_OL_re_010.csv')




##########################
#CD markers:BEAT and TCGA

#CD markers:BEAT and TCGA
CD_genes_mtwt_ix<-grep('^(CD\\d+)',(DEG.BEAT.TCGA.mtwt.shared$DEG))
CD_genes_mtlwt_ix<-grep('^(CD\\d+)',(DEG.BEAT.TCGA.mtlwt.shared$DEG))

DEG.BEAT.TCGA.mtwt.shared.CD<-DEG.BEAT.TCGA.mtwt.shared[CD_genes_mtwt_ix,]
DEG.BEAT.TCGA.mtlwt.shared.CD<-DEG.BEAT.TCGA.mtlwt.shared[CD_genes_mtlwt_ix,]
write.csv(DEG.BEAT.TCGA.mtwt.shared.CD,'Analysis/GEP/DEG_GSEA/DEG/BEAT_TCGA_OL_MTWT_DEG_qlf_CD_fin.csv')
write.csv(DEG.BEAT.TCGA.mtlwt.shared.CD,'Analysis/GEP/DEG_GSEA/DEG/BEAT_TCGA_OL_MTlikeWT_DEG_qlf_CD_fin.csv')

#make 2 verions of heatmap (FDR<0.1 on one side FDR<0.05 in both sides)
BEAT.TCGA.mtwt.CDid.FDR005<-DEG.BEAT.TCGA.mtwt.shared.CD%>% filter(DEG01_FDR<0.05 &DEG02_FDR<0.05)
BEAT.TCGA.mtlwt.CDid.FDR005<-DEG.BEAT.TCGA.mtlwt.shared.CD%>% filter(DEG01_FDR<0.05 &DEG02_FDR<0.05)

shCD_mtmtl<-intersect(BEAT.TCGA.mtwt.CDid.FDR005$DEG,BEAT.TCGA.mtlwt.CDid.FDR005$DEG)
rownames(BEAT.TCGA.mtwt.CDid.FDR005)<-BEAT.TCGA.mtwt.CDid.FDR005$DEG
rownames(BEAT.TCGA.mtlwt.CDid.FDR005)<-BEAT.TCGA.mtlwt.CDid.FDR005$DEG
BEAT.TCGA.mtwt.CDid.FDR005[shCD_mtmtl,]
BEAT.TCGA.mtlwt.CDid.FDR005[shCD_mtmtl,]
BEAT.TCGA.mtwt.CDid.FDR01<-DEG.BEAT.TCGA.mtwt.shared.CD
BEAT.TCGA.mtlwt.CDid.FDR01<-DEG.BEAT.TCGA.mtlwt.shared.CD

#generate heatmap: mt vs. WT CD marker 
y.beat.cpm #log norm
y.tcga.cpm #log norm
mt_ix03<-which(colnames(y.beat.cpm) %in% BEAT_MUT_id)
mtl_ix03<-which(colnames(y.beat.cpm) %in% BEAT_MUT_like_id)
wt_ix03<-setdiff(1:dim(y.beat.cpm)[2],c(mt_ix03,mtl_ix03))
#cls.beat<-rep('TP53_WT',dim(y.beat.cpm)[2]);cls.beat[mt_ix03]<-'TP53_MUT';cls.beat[mtl_ix03]<-'TP53MUT_like'
cls.beat<-rep('TP53_WT',length(mt_ix03)+length(wt_ix03));cls.beat[1:36]<-'TP53_MUT'
#cls.beat.df<-data.frame('ID'=colnames(y.beat.cpm),'cls'= cls.beat ) %>% column_to_rownames(var='ID')
cls.beat.df<-data.frame('ID'=c(BEAT_MUT_id,BEAT_WT_id),'cls'= cls.beat ) %>% column_to_rownames(var='ID')
anno_beat_gep<-HeatmapAnnotation(df=cls.beat.df,
                                 col = list(cls=c('TP53_MUT'='red','TP53_WT'='blue1')))
beat.mtwt.cd.fdr005<-Heatmap(t(scale(t(y.beat.cpm[BEAT.TCGA.mtwt.CDid.FDR005$DEG,c(BEAT_MUT_id,BEAT_WT_id)]))), 
        show_row_names = T,show_column_names = F,top_annotation = anno_beat_gep,
        column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(8, "cm"),
        row_names_gp = gpar(fontsize = 16),
        show_row_dend = F,cluster_rows = F,cluster_columns =F, 
       heatmap_legend_param  = list(title=NULL,legend_height = unit(4, "cm")))
beat.mtwt.cd.fdr01<-Heatmap(t(scale(t(y.beat.cpm[BEAT.TCGA.mtwt.CDid.FDR01$DEG,c(BEAT_MUT_id,BEAT_WT_id)]))), 
        show_row_names = T,show_column_names = F,top_annotation = anno_beat_gep,
        column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 16),
        show_row_dend = F,cluster_rows = F,cluster_columns =F,
        heatmap_legend_param  = list(title=NULL,legend_height = unit(4, "cm")))

TCGA_MUT_id
TCGA_MUT_id<-gsub('.03A.*|.03B.*','',TCGA_clinical_final$LabId.y[TCGA_clinical_final$cls.GLM.y=='TP53_MUT'])
TCGA_MUT_like_id<-gsub('.03A.*|.03B.*','',TCGA_clinical_final$LabId.y[TCGA_clinical_final$cls.GLM.y=='TP53MUT_like'])
TCGA_WT_id<-gsub('.03A.*|.03B.*','',TCGA_clinical_final$LabId.y[TCGA_clinical_final$cls.GLM.y=='TP53_WT'])

mt_ix03.tcga<-which(colnames(y.tcga.cpm) %in% TCGA_MUT_id)
mtl_ix03.tcga<-which(colnames(y.tcga.cpm) %in% TCGA_MUT_like_id)
wt_ix03.tcga<-setdiff(1:dim(y.tcga.cpm)[2],c(mt_ix03.tcga,mtl_ix03.tcga))
#cls.tcga<-rep('TP53_WT',dim(y.tcga.cpm)[2]);cls.tcga[mt_ix03.tcga]<-'TP53_MUT';cls.tcga[mtl_ix03.tcga]<-'TP53MUT_like'
cls.tcga<-rep('TP53_WT',length(mt_ix03.tcga)+length(wt_ix03.tcga));cls.tcga[1:15]<-'TP53_MUT'#;cls.tcga[mtl_ix03.tcga]<-'TP53MUT_like'
#cls.tcga.df<-data.frame('ID'=colnames(y.tcga.cpm),'cls'= cls.tcga ) %>% column_to_rownames(var='ID')
cls.tcga.df<-data.frame('ID'=c(TCGA_MUT_id,TCGA_WT_id),'cls'= cls.tcga ) %>% column_to_rownames(var='ID')
anno_tcga_gep<-HeatmapAnnotation(df=cls.tcga.df,
                                 col = list(cls=c('TP53_MUT'='red','TP53_WT'='blue1')))
tcga.mtwt.cd.fdr005<-Heatmap(scale_GEP(y.tcga.cpm[BEAT.TCGA.mtwt.CDid.FDR005$DEG,c(TCGA_MUT_id,TCGA_WT_id)]), 
        show_row_names = T,show_column_names = F,top_annotation = anno_tcga_gep,
        column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 16),
        show_row_dend = F,cluster_rows = F,cluster_columns =F,
        heatmap_legend_param  = list(title=NULL,legend_height = unit(4, "cm")))

tcga.mtwt.cd.fdr01<-Heatmap(scale_GEP(y.tcga.cpm[BEAT.TCGA.mtwt.CDid.FDR01$DEG,c(TCGA_MUT_id,TCGA_WT_id)]), 
        show_row_names = T,show_column_names = F,top_annotation = anno_tcga_gep,
        column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 16),
        show_row_dend = F,cluster_rows = F,cluster_columns =F,
        heatmap_legend_param  = list(title=NULL,legend_height = unit(4, "cm")))


#mtlwt
cls.beat<-rep('TP53_WT',length(mtl_ix03)+length(wt_ix03));cls.beat[1:40]<-'TP53_MUTlike'
#cls.beat.df<-data.frame('ID'=colnames(y.beat.cpm),'cls'= cls.beat ) %>% column_to_rownames(var='ID')
cls.beat.df<-data.frame('ID'=c(BEAT_MUT_like_id,BEAT_WT_id),'cls'= cls.beat ) %>% column_to_rownames(var='ID')
anno_beat_gep<-HeatmapAnnotation(df=cls.beat.df,
                                 col = list(cls=c('TP53_MUTlike'='purple','TP53_WT'='blue1')))
beat.mtlwt.cd.fdr005<-Heatmap(t(scale(t(y.beat.cpm[BEAT.TCGA.mtlwt.CDid.FDR005$DEG,c(BEAT_MUT_like_id,BEAT_WT_id)]))), 
                             show_row_names = T,show_column_names = F,top_annotation = anno_beat_gep,
                             column_names_max_height = unit(4, "cm"),
                             column_names_gp = gpar(fontsize = 6),
                             row_names_max_width = unit(6, "cm"),
                             row_names_gp = gpar(fontsize = 16),
                             show_row_dend = F,cluster_rows = F,cluster_columns =F,
                             heatmap_legend_param  = list(title=NULL,legend_height = unit(4, "cm")))
beat.mtlwt.cd.fdr01<-Heatmap(t(scale(t(y.beat.cpm[BEAT.TCGA.mtlwt.CDid.FDR01$DEG,c(BEAT_MUT_like_id,BEAT_WT_id)]))), 
                            show_row_names = T,show_column_names = F,top_annotation = anno_beat_gep,
                            column_names_max_height = unit(4, "cm"),
                            column_names_gp = gpar(fontsize = 6),
                            row_names_max_width = unit(6, "cm"),
                            row_names_gp = gpar(fontsize = 16),
                            show_row_dend = F,cluster_rows = F,cluster_columns =F,
                            heatmap_legend_param  = list(title=NULL,legend_height = unit(4, "cm")))

TCGA_MUT_id

# TCGA_MUT_id<-gsub('.03A.*|.03B.*','',TCGA_clinical_final$LabId.y[TCGA_clinical_final$cls.GLM.y=='TP53_MUT'])
# TCGA_MUT_like_id<-gsub('.03A.*|.03B.*','',TCGA_clinical_final$LabId.y[TCGA_clinical_final$cls.GLM.y=='TP53MUT_like'])
# TCGA_WT_id<-gsub('.03A.*|.03B.*','',TCGA_clinical_final$LabId.y[TCGA_clinical_final$cls.GLM.y=='TP53_WT'])
# 
# mt_ix03.tcga<-which(colnames(y.tcga.cpm) %in% TCGA_MUT_id)
# mtl_ix03.tcga<-which(colnames(y.tcga.cpm) %in% TCGA_MUT_like_id)
# wt_ix03.tcga<-setdiff(1:dim(y.tcga.cpm)[2],c(mt_ix03.tcga,mtl_ix03.tcga))

#cls.tcga<-rep('TP53_WT',dim(y.tcga.cpm)[2]);cls.tcga[mt_ix03.tcga]<-'TP53_MUT';cls.tcga[mtl_ix03.tcga]<-'TP53MUT_like'
cls.tcga<-rep('TP53_WT',length(mtl_ix03.tcga)+length(wt_ix03.tcga));cls.tcga[1:23]<-'TP53_MUTlike'#;cls.tcga[mtl_ix03.tcga]<-'TP53MUT_like'
#cls.tcga.df<-data.frame('ID'=colnames(y.tcga.cpm),'cls'= cls.tcga ) %>% column_to_rownames(var='ID')
cls.tcga.df<-data.frame('ID'=c(TCGA_MUT_like_id,TCGA_WT_id),'cls'= cls.tcga ) %>% column_to_rownames(var='ID')
anno_tcga_gep<-HeatmapAnnotation(df=cls.tcga.df,
                                 col = list(cls=c('TP53_MUTlike'='purple','TP53_WT'='blue1')))
tcga.mtlwt.cd.fdr005<-Heatmap(scale_GEP(y.tcga.cpm[BEAT.TCGA.mtlwt.CDid.FDR005$DEG,c(TCGA_MUT_like_id,TCGA_WT_id)]), 
                             show_row_names = T,show_column_names = F,top_annotation = anno_tcga_gep,
                             column_names_max_height = unit(4, "cm"),
                             column_names_gp = gpar(fontsize = 6),
                             row_names_max_width = unit(6, "cm"),
                             row_names_gp = gpar(fontsize = 16),
                             show_row_dend = F,cluster_rows = F,cluster_columns =F,
                             heatmap_legend_param  = list(title=NULL,legend_height = unit(4, "cm")))

tcga.mtlwt.cd.fdr01<-Heatmap(scale_GEP(y.tcga.cpm[BEAT.TCGA.mtlwt.CDid.FDR01$DEG,c(TCGA_MUT_like_id,TCGA_WT_id)]), 
                            show_row_names = T,show_column_names = F,top_annotation = anno_tcga_gep,
                            column_names_max_height = unit(4, "cm"),
                            column_names_gp = gpar(fontsize = 6),
                            row_names_max_width = unit(6, "cm"),
                            row_names_gp = gpar(fontsize = 16),
                            show_row_dend = F,cluster_rows = F,cluster_columns =F,
                            heatmap_legend_param  = list(title=NULL,legend_height = unit(4, "cm")))


jpeg('Analysis/GEP/DEG_GSEA/DEG/Heatmap/beat_mtwt_cd_fdr005.jpeg',height = 7, width = 7, units = 'in', res=300)
print(beat.mtwt.cd.fdr005)
dev.off()

jpeg('Analysis/GEP/DEG_GSEA/DEG/Heatmap/beat_mtwt_cd_fdr01.jpeg',height = 7, width = 7, units = 'in', res=300)
print(beat.mtwt.cd.fdr01)
dev.off()

jpeg('Analysis/GEP/DEG_GSEA/DEG/Heatmap/tcga_mtwt_cd_fdr005.jpeg',height = 7, width = 7, units = 'in', res=300)
print(tcga.mtwt.cd.fdr005)
dev.off()

jpeg('Analysis/GEP/DEG_GSEA/DEG/Heatmap/tcga_mtwt_cd_fdr01.jpeg',height = 7, width = 7, units = 'in', res=300)
print(tcga.mtwt.cd.fdr01)
dev.off()

##
jpeg('Analysis/GEP/DEG_GSEA/DEG/Heatmap/beat_mtlwt_cd_fdr005.jpeg',height = 7, width = 7, units = 'in', res=300)
print(beat.mtlwt.cd.fdr005)
dev.off()

jpeg('Analysis/GEP/DEG_GSEA/DEG/Heatmap/beat_mtlwt_cd_fdr01.jpeg',height = 7, width = 7, units = 'in', res=300)
print(beat.mtlwt.cd.fdr01)
dev.off()

jpeg('Analysis/GEP/DEG_GSEA/DEG/Heatmap/tcga_mtlwt_cd_fdr005.jpeg',height = 7, width = 7, units = 'in', res=300)
print(tcga.mtlwt.cd.fdr005)
dev.off()

jpeg('Analysis/GEP/DEG_GSEA/DEG/Heatmap/tcga_mtlwt_cd_fdr01.jpeg',height = 7, width = 7, units = 'in', res=300)
print(tcga.mtlwt.cd.fdr01)
dev.off()








#msigdb_filt<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/GSEA/Compiled data/Gene set db/MSIGDB_ZS_gs_filtered_long_complete_gene_converted.rds')
msigdb_filt<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/GSEA/Compiled data/Gene set db/MSIGDB_ZS_gs_filtered_long_complete_gene_conv2022.rds')

GSEA.beat.mtwt.qlf<-cprofiler_GSEA(DEG.beat.mtwt.qlf,msigdb_filt)
GSEA.beat.mtlwt.qlf<-cprofiler_GSEA(DEG.beat.mtlwt.qlf,msigdb_filt)

GSEA.tcga.mtwt.qlf<-cprofiler_GSEA(DEG.tcga.mtwt.qlf,msigdb_filt)
GSEA.tcga.mtlwt.qlf<-cprofiler_GSEA(DEG.tcga.mtlwt.qlf,msigdb_filt)


write.csv(GSEA.beat.mtwt.qlf,'Analysis/GEP/DEG_GSEA/GSEA/BEAT_MTWT_DEG_qlf_fin.csv')
write.csv(GSEA.beat.mtlwt.qlf,'Analysis/GEP/DEG_GSEA/GSEA/BEAT_MTlikeWT_DEG_qlf_fin.csv')

write.csv(GSEA.tcga.mtwt.qlf,'Analysis/GEP/DEG_GSEA/GSEA/TCGA_MTWT_DEG_qlf_fin.csv')
write.csv(GSEA.tcga.mtlwt.qlf,'Analysis/GEP/DEG_GSEA/GSEA/TCGA_MTlikeWT_DEG_qlf_fin.csv')

#save.image('Analysis/GEP/DEG_GSEA/TP53GEP_DEG_GSEA01.Rdata')
load('Analysis/GEP/DEG_GSEA/TP53GEP_DEG_GSEA01.Rdata')

source("C:/Users/yklee/Desktop/Sachs Lab/R_project/TP53/TP53GEP_functions.r")
packages<-c('edgeR','limma','dplyr','tidyr','sva','tidyverse',
            'ggplot2','survival','survminer','HGNChelper','stringr','ComplexHeatmap')
for(p in packages){
  library(p,character.only = T)
}
rm(fit.beat.mtwt,fit.beat.mtlwt,fit.tcga.mtwt,fit.tcga.mtlwt,
   qlf.beat.mtwt.qlf,qlf.beat.mtwt.treat,qlf.beat.mtlwt.qlf,qlf.beat.mtlwt.treat,
   qlf.tcga.mtwt.qlf,qlf.tcga.mtwt.treat,qlf.tcga.mtlwt.qlf,qlf.tcga.mtlwt.treat,
   BEAT_RNA_counts,TCGA_AML_counts)
gc()





shared_gs_updn_sig<-function(gs_01.df,gs_02.df,pval_cutoff){
  #pval_cutoff
  GS_id_OL<-intersect(gs_01.df$ID[gs_01.df$p.adjust<pval_cutoff],
                      gs_02.df$ID[gs_02.df$p.adjust<pval_cutoff])
  shared.all<-data.frame('GS_ID'=GS_id_OL,'gs_01_NES'=rep(NA,length(GS_id_OL)),'gs_01_pval'=rep(NA,length(GS_id_OL)),
                         'gs_02_NES'=rep(NA,length(GS_id_OL)),'gs_02_pval'=rep(NA,length(GS_id_OL)))
  for(i in 1:length(shared.all$GS_ID)){
    gs01_ix<-which(gs_01.df$ID %in% shared.all$GS_ID[i])
    gs02_ix<-which(gs_02.df$ID %in% shared.all$GS_ID[i])
    shared.all$gs_01_NES[i]<-gs_01.df$NES[gs01_ix]
    shared.all$gs_02_NES[i]<-gs_02.df$NES[gs02_ix]
    shared.all$gs_01_pval[i]<-gs_01.df$p.adjust[gs01_ix]
    shared.all$gs_02_pval[i]<-gs_02.df$p.adjust[gs02_ix]
  }
  shared.all.clean<-shared.all[c(which(shared.all$gs_01_NES>0 & shared.all$gs_02_NES>0),
                                 which(shared.all$gs_01_NES<0 & shared.all$gs_02_NES<0)),]
  # shared.gs.up
  # shared.gs.dn
  return(shared.all.clean)
}
#3 shared gene sets
#1.shared GSEA BEAT and TCGA
##a. MT vs. WT
##b. MTL vs. WT
#2. shared GSEA MTvsWT and MTLvsWT
##a. BEAT
##b. TCGA
##a. BEAT and TCGA

#1.shared GSEA BEAT and TCGA
##a. MT vs. WT
GSEA.beat.tcga.mtwt.qlf01<-shared_gs_updn_sig(GSEA.beat.mtwt.qlf,GSEA.tcga.mtwt.qlf,0.2)
GSEA.beat.tcga.mtwt.qlf005<-GSEA.beat.tcga.mtwt.qlf01%>% filter(gs_01_pval<0.05 & gs_02_pval<0.05)
##b. MTL vs. WT
GSEA.beat.tcga.mtlwt.qlf01<-shared_gs_updn_sig(GSEA.beat.mtlwt.qlf,GSEA.tcga.mtlwt.qlf,0.2)

#2. shared GSEA MTvsWT and MTLvsWT
##a. BEAT
GSEA.beat.mtwt_vs_mtlwt.qlf02<-shared_gs_updn_sig(GSEA.beat.mtwt.qlf,GSEA.beat.mtlwt.qlf,0.2)
##b. TCGA
GSEA.tcga.mtwt_vs_mtlwt.qlf02<-shared_gs_updn_sig(GSEA.tcga.mtwt.qlf,GSEA.tcga.mtlwt.qlf,0.2)
##a. BEAT and TCGA
rownames(GSEA.beat.mtwt_vs_mtlwt.qlf02)<-GSEA.beat.mtwt_vs_mtlwt.qlf02$GS_ID
rownames(GSEA.tcga.mtwt_vs_mtlwt.qlf02)<-GSEA.tcga.mtwt_vs_mtlwt.qlf02$GS_ID
ol_gs<-(intersect(GSEA.beat.mtwt_vs_mtlwt.qlf02$GS_ID,GSEA.tcga.mtwt_vs_mtlwt.qlf02$GS_ID))
View(GSEA.beat.mtwt_vs_mtlwt.qlf02[ol_gs,])
View(GSEA.tcga.mtwt_vs_mtlwt.qlf02[ol_gs,])

write.csv(GSEA.beat.tcga.mtwt.qlf01,'Analysis/GEP/DEG_GSEA/GSEA/shGSEA_BEAT_TCGA_mtwt_FDR02.csv')
write.csv(GSEA.beat.tcga.mtlwt.qlf01,'Analysis/GEP/DEG_GSEA/GSEA/shGSEA_BEAT_TCGA_mtlwt_FDR02.csv')

write.csv(GSEA.beat.mtwt_vs_mtlwt.qlf02,'Analysis/GEP/DEG_GSEA/GSEA/shGSEA_BEAT_mtwt_mtlwt_FDR02.csv')
write.csv(GSEA.tcga.mtwt_vs_mtlwt.qlf02,'Analysis/GEP/DEG_GSEA/GSEA/shGSEA_TCGA_mtwt_mtlwt_FDR02.csv')

write.csv(GSEA.beat.mtwt_vs_mtlwt.qlf02[ol_gs,],'Analysis/GEP/DEG_GSEA/GSEA/shGSEA_BEAT_TCGA_mtwt_mtlwt_FDR02_beatver.csv')
write.csv(GSEA.tcga.mtwt_vs_mtlwt.qlf02[ol_gs,],'Analysis/GEP/DEG_GSEA/GSEA/shGSEA_BEAT_TCGA_mtwt_mtlwt_FDR02_tcgaver.csv')

LSC_list<-c('HSC I','HSC II','HSC III','LGMP'	,'SELF RENEWAL UP',	'SELF RENEWAL DN',	'E2F1 TARGETS REPRESSED BY SERUM',
            'Somervaille LSC UP',	'Somervaille LSC Dn','UP when Myb OFF (lo FC)',	'DN when myb OFF (hi FC)',
            'UP when MLLAF9 OFF (Lo FC)',	'DN when MLLAF9 OFF (hi FC)',	'UP when Myb OFF (hi FC)',
            'DN when myb OFF (lo FC)',	'UP when MLLAF9 OFF (Hi FC)',	'DN when MLLAF9 OFF (lo FC)',
            'HSC_MATURE_FETAL',	'Hu LSC UP',	'Hu LSC Dn',	'Top HSC',	'Ng et al up',	'Ng et al dn',
            'Group1',	'Group3','TP53_proteosome_sig37','TP53_proteosome_72sig_up','TP53_proteosome_72sig_dn',
            'TP53_proteosome_205sig_up','TP53_proteosome_205sig_dn')

extract_subGS<-function(cprofile_res,custom_list,shared=1){
  if (shared==1){
    KEGG_gs_mtwt_ix<-grep('^(KEGG)',(cprofile_res$GS_ID))
    GOBP_gs_mtwt_ix<-grep('^(GOBP_)',(cprofile_res$GS_ID))
    WP_gs_mtwt_ix<-grep('^(WP)',(cprofile_res$GS_ID))
    LSC_gs_mtwt_ix<-which((cprofile_res$GS_ID)%in%LSC_list)
  }else{
    KEGG_gs_mtwt_ix<-grep('^(KEGG)',(cprofile_res$ID))
    GOBP_gs_mtwt_ix<-grep('^(GOBP_)',(cprofile_res$ID))
    WP_gs_mtwt_ix<-grep('^(WP)',(cprofile_res$ID))
    LSC_gs_mtwt_ix<-which((cprofile_res$ID)%in%LSC_list)
  }
  
  
  
  GSEA.KEGG<-cprofile_res[KEGG_gs_mtwt_ix,]
  
  GSEA.GOBP<-cprofile_res[GOBP_gs_mtwt_ix,]
  
  GSEA.WP<-cprofile_res[WP_gs_mtwt_ix,]
  
  GSEA.LSC<-cprofile_res[LSC_gs_mtwt_ix,]
  
  GS.sub <- list("KEGG" = GSEA.KEGG, "GOBP" = GSEA.GOBP,
                 'WP'=GSEA.WP,'LSC'=GSEA.LSC)
  return(GS.sub)
}

KEGG_gs_mtwt_ix<-grep('^(KEGG)',(GSEA.beat.tcga.mtwt.qlf01$GS_ID))
KEGG_gs_mtlwt_ix<-grep('^(KEGG)',(GSEA.beat.tcga.mtlwt.qlf01$GS_ID))
GOBP_gs_mtwt_ix<-grep('^(GOBP_)',(GSEA.beat.tcga.mtwt.qlf01$GS_ID))
GOBP_gs_mtlwt_ix<-grep('^(GOBP_)',(GSEA.beat.tcga.mtlwt.qlf01$GS_ID))
WP_gs_mtwt_ix<-grep('^(WP)',(GSEA.beat.tcga.mtwt.qlf01$GS_ID))
WP_gs_mtlwt_ix<-grep('^(WP)',(GSEA.beat.tcga.mtlwt.qlf01$GS_ID))
LSC_gs_mtwt_ix<-which((GSEA.beat.tcga.mtwt.qlf01$GS_ID)%in%LSC_list)
LSC_gs_mtlwt_ix<-which((GSEA.beat.tcga.mtlwt.qlf01$GS_ID)%in%LSC_list)


GSEA.beat.tcga.mtwt.qlf01.KEGG<-GSEA.beat.tcga.mtwt.qlf01[KEGG_gs_mtwt_ix,]
GSEA.beat.tcga.mtlwt.qlf01.KEGG<-GSEA.beat.tcga.mtlwt.qlf01[KEGG_gs_mtlwt_ix,]

GSEA.beat.tcga.mtwt.qlf01.GOBP<-GSEA.beat.tcga.mtwt.qlf01[GOBP_gs_mtwt_ix,]
GSEA.beat.tcga.mtlwt.qlf01.GOBP<-GSEA.beat.tcga.mtlwt.qlf01[GOBP_gs_mtlwt_ix,]

GSEA.beat.tcga.mtwt.qlf01.WP<-GSEA.beat.tcga.mtwt.qlf01[WP_gs_mtwt_ix,]
GSEA.beat.tcga.mtlwt.qlf01.WP<-GSEA.beat.tcga.mtlwt.qlf01[WP_gs_mtlwt_ix,]

GSEA.beat.tcga.mtwt.qlf01.LSC<-GSEA.beat.tcga.mtwt.qlf01[LSC_gs_mtwt_ix,]
GSEA.beat.tcga.mtlwt.qlf01.LSC<-GSEA.beat.tcga.mtlwt.qlf01[LSC_gs_mtlwt_ix,]

GS.sub <- list("KEGG" = GSEA.beat.tcga.mtwt.qlf01.KEGG, "GOBP" = GSEA.beat.tcga.mtwt.qlf01.GOBP,
               'WP'=GSEA.beat.tcga.mtwt.qlf01.WP,'LSC'=GSEA.beat.tcga.mtwt.qlf01.LSC)

require(openxlsx)
GS.shared.mtwt<-extract_subGS(GSEA.beat.tcga.mtwt.qlf01,LSC_list,shared=1)
GS.BEAT.sub.mtwt<-extract_subGS(GSEA.beat.mtwt.qlf,LSC_list,shared=0)
GS.TCGA.sub.mtwt<-extract_subGS(GSEA.tcga.mtwt.qlf,LSC_list,shared=0)

write.xlsx(GS.shared.mtwt, file = "Analysis/GEP/DEG_GSEA/GSEA/GSEA_BEAT_TCGA_shared_MTWT_cat.xlsx")
write.xlsx(GS.BEAT.sub.mtwt, file = "Analysis/GEP/DEG_GSEA/GSEA/GSEA_BEAT_MTWT_cat.xlsx")
write.xlsx(GS.TCGA.sub.mtwt, file = "Analysis/GEP/DEG_GSEA/GSEA/GSEA_TCGA_MTWT_cat.xlsx")

GS.shared.mtlwt <- list("KEGG" = GSEA.beat.tcga.mtlwt.qlf01.KEGG, "GOBP" = GSEA.beat.tcga.mtlwt.qlf01.GOBP,
                        'WP'=GSEA.beat.tcga.mtlwt.qlf01.WP,'LSC'=GSEA.beat.tcga.mtlwt.qlf01.LSC)
write.xlsx(GS.shared.mtlwt, file = "Analysis/GEP/DEG_GSEA/GSEA/GSEA_BEAT_TCGA_MTlikeWT_cat.xlsx")
#save.image('Analysis/GEP/DEG_GSEA/TP53GEP_DEG_GSEA02.Rdata')
load('Analysis/GEP/DEG_GSEA/TP53GEP_DEG_GSEA02.Rdata')

# 
# library(xlsx)
# write.xlsx(dataframe1, file="filename.xlsx", sheetName="sheet1", row.names=FALSE)
# write.xlsx(dataframe2, file="filename.xlsx", sheetName="sheet2", append=TRUE, row.names=FALSE)

#2022.08.18-->2022.11.09
#DEG shared between Mutlike vs. WT and Mut vs. WT
# mt_mtl_deg<-intersect(DEG.BEAT.TCGA.mtwt.shared$DEG,
#           DEG.BEAT.TCGA.mtlwt.shared$DEG)
mt_mtl_deg<-intersect(DEG.BEAT.TCGA.mtwt.shared.filt$DEG,
                      DEG.BEAT.TCGA.mtlwt.shared.filt$DEG)
# rownames(DEG.BEAT.TCGA.mtwt.shared)<-DEG.BEAT.TCGA.mtwt.shared$DEG
# rownames(DEG.BEAT.TCGA.mtlwt.shared)<-DEG.BEAT.TCGA.mtlwt.shared$DEG
rownames(DEG.BEAT.TCGA.mtwt.shared.filt)<-DEG.BEAT.TCGA.mtwt.shared.filt$DEG
rownames(DEG.BEAT.TCGA.mtlwt.shared.filt)<-DEG.BEAT.TCGA.mtlwt.shared.filt$DEG

View(DEG.BEAT.TCGA.mtwt.shared.filt[mt_mtl_deg,])
View(DEG.BEAT.TCGA.mtlwt.shared.filt[mt_mtl_deg,])
conc_ix<-c()
#for(i in 1:dim(mt_mtl_deg)){}
# pos_id<-DEG.BEAT.TCGA.mtwt.shared[mt_mtl_deg,]$DEG[(DEG.BEAT.TCGA.mtwt.shared[mt_mtl_deg,]$DEG01_FC>0 & DEG.BEAT.TCGA.mtlwt.shared[mt_mtl_deg,]$DEG01_FC>0)]
# neg_id<-DEG.BEAT.TCGA.mtwt.shared[mt_mtl_deg,]$DEG[(DEG.BEAT.TCGA.mtwt.shared[mt_mtl_deg,]$DEG01_FC<0 & DEG.BEAT.TCGA.mtlwt.shared[mt_mtl_deg,]$DEG01_FC<0)]
pos_id<-DEG.BEAT.TCGA.mtwt.shared.filt[mt_mtl_deg,]$DEG[(DEG.BEAT.TCGA.mtwt.shared.filt[mt_mtl_deg,]$DEG01_FC>0 & DEG.BEAT.TCGA.mtwt.shared.filt[mt_mtl_deg,]$DEG01_FC>0)]
neg_id<-DEG.BEAT.TCGA.mtwt.shared.filt[mt_mtl_deg,]$DEG[(DEG.BEAT.TCGA.mtwt.shared.filt[mt_mtl_deg,]$DEG01_FC<0 & DEG.BEAT.TCGA.mtwt.shared.filt[mt_mtl_deg,]$DEG01_FC<0)]

View(DEG.BEAT.TCGA.mtwt.shared.filt[mt_mtl_deg,][c(pos_id,neg_id),])
View(DEG.BEAT.TCGA.mtlwt.shared.filt[mt_mtl_deg,][c(pos_id,neg_id),])
#do both BEAT and TCGA but do BEAT first 

# y.beat.cpm<-as.data.frame(cpm(y.beat,log=T,normalized.lib.sizes = TRUE))
# y.tcga.cpm<-as.data.frame(cpm(y.tcga,log=T,normalized.lib.sizes = TRUE))
y.beat.cpm
y.tcga.cpm


mt_ix03<-which(colnames(y.beat.cpm) %in% BEAT_MUT_id)
mtl_ix03<-which(colnames(y.beat.cpm) %in% BEAT_MUT_like_id)
wt_ix03<-setdiff(1:dim(y.beat.cpm)[2],c(mt_ix03,mtl_ix03))
cls.beat<-rep('TP53_WT',dim(y.beat.cpm)[2]);cls.beat[mt_ix03]<-'TP53_MUT';cls.beat[mtl_ix03]<-'TP53MUT_like'
cls.beat.df<-data.frame('ID'=colnames(y.beat.cpm),'cls'= cls.beat ) %>% column_to_rownames(var='ID')
anno_beat_gep<-HeatmapAnnotation(df=cls.beat.df,
                                 col = list(cls=c('TP53_MUT'='red','TP53_WT'='blue1','TP53MUT_like'='purple')))
beat.mtwt.mtlwt.DEG<-Heatmap(t(scale(t(y.beat.cpm[c(pos_id,neg_id),c(BEAT_MUT_id,BEAT_MUT_like_id,BEAT_WT_id)]))), 
        show_row_names = F,show_column_names = F,top_annotation = anno_beat_gep[c(mt_ix03,mtl_ix03,wt_ix03)],
        column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 12),
        show_row_dend = F,cluster_rows = F,cluster_columns =F)
TCGA_MUT_id
colnames(y.tcga.cpm)<-gsub('.03A.*|.03B.*','',colnames(y.tcga.cpm))
# mt_ix03.tcga<-which(gsub('.03A.*|.03B.*','',colnames(y.tcga.cpm)) %in% TCGA_MUT_id)
# mtl_ix03.tcga<-which(gsub('.03A.*|.03B.*','',colnames(y.tcga.cpm)) %in% TCGA_MUT_like_id)
# wt_ix03.tcga<-setdiff(1:dim(y.tcga.cpm)[2],c(mt_ix03.tcga,mtl_ix03.tcga))
cls.tcga<-rep('TP53_WT',dim(y.tcga.cpm)[2]);cls.tcga[mt_ix03.tcga]<-'TP53_MUT';cls.tcga[mtl_ix03.tcga]<-'TP53MUT_like'
cls.tcga.df<-data.frame('ID'=colnames(y.tcga.cpm),'cls'= cls.tcga ) %>% column_to_rownames(var='ID')
anno_tcga_gep<-HeatmapAnnotation(df=cls.tcga.df,
                                 col = list(cls=c('TP53_MUT'='red','TP53_WT'='blue1','TP53MUT_like'='purple')))
tcga.mtwt.mtlwt.DEG<-Heatmap(scale_GEP(y.tcga.cpm[c(pos_id,neg_id),c(TCGA_MUT_id,TCGA_MUT_like_id,TCGA_WT_id)]), 
        show_row_names = F,show_column_names = F,top_annotation = anno_tcga_gep[c(mt_ix03.tcga,mtl_ix03.tcga,wt_ix03.tcga)],
        column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 12),
        show_row_dend = F,cluster_rows = F,cluster_columns =F)
jpeg('Analysis/GEP/DEG_GSEA/DEG/Heatmap/beat_mtwt_mtlwt_DEG_fdr005.jpeg',height = 7, width = 7, units = 'in', res=300)
print(beat.mtwt.mtlwt.DEG)
dev.off()
jpeg('Analysis/GEP/DEG_GSEA/DEG/Heatmap/tcga_mtwt_mtlwt_DEG_fdr005.jpeg',height = 7, width = 7, units = 'in', res=300)
print(tcga.mtwt.mtlwt.DEG)
dev.off()

####################################################################################
#2022.12.02
# Zohar and I annotated GSEA gene sets what to use in paper - 
#   shGSEA_BEAT_TCGA_mtwt_FDR0.05_ZS_Ylannot
# shGSEA_BEAT_TCGA_mtlwt_FDR010_YL_annot_Zsreview
# Use these to generate GSEA plots - BEAT AML as representative
#beat AML as representative

#retrieve original GSEA object
cprofiler_GSEA<-function(Seurat_DE_genes,msig_database_filt){
  library(dplyr)
  #calculate weight-ranks for fgsea calculation:
  
  Seurat_DE_genes$gs_weight<- (-log10(Seurat_DE_genes$PValue)*sign(Seurat_DE_genes$logFC))
  na_ix<-which(is.na(Seurat_DE_genes$gs_weight))
  if(any(na_ix)){
    Seurat_DE_genes<-Seurat_DE_genes[-na_ix,]
  }else{
    Seurat_DE_genes<-Seurat_DE_genes
  }
  
  #msigdb input outside of function
  #m_df<-msigdbr(species = 'Homo sapiens')# retrive all the database :C1~c7
  #I usually take: H(hallmark), C2(curated), C5(GO) and C6(oncogenic) + ZScurated lists+C4
  
  ##msig_database_filt
  
  Seurat_DE_genes_gs<-(Seurat_DE_genes)%>%
    dplyr::arrange(desc(gs_weight)) %>%
    tibble::rownames_to_column(var = 'gene_id') %>%
    dplyr::select(gene_id,gs_weight)
  
  
  
  
  dup_ix<-which(duplicated(Seurat_DE_genes_gs$gs_weight))
  print(length(dup_ix))
  #Seurat_DE_genes_gs_ori<-Seurat_DE_genes_gs
  tied_gs_weights<-Seurat_DE_genes_gs %>% group_by(gs_weight) %>% dplyr::summarize(n=n()) %>% filter(n>1) %>% arrange(desc(gs_weight))
  dup_ix_all<-c()
  
  for (i in 1:dim(tied_gs_weights)[1]){
    #print(as.numeric(tied_gs_weights$gs_weight[i]),digits = 20)
    ix<-which(tied_gs_weights$gs_weight[i] == Seurat_DE_genes_gs$gs_weight)
    #print(ix)
    dup_ix_all<-c(dup_ix_all,ix)}
  print(paste('duplicated genes:',length(dup_ix),sep=''))
  print(paste('all duplicated events:',length(dup_ix_all),sep=''))
  print(paste('percentages of duplicates:',(length(dup_ix_all)/dim(Seurat_DE_genes)[1])*100,'%',sep=''))
  
  pos_rank<-(Seurat_DE_genes_gs[Seurat_DE_genes_gs$gs_weight>=0,])%>% mutate(rank = rank(gs_weight,ties.method ='random')) %>%
    arrange(desc(rank))
  neg_rank<-(Seurat_DE_genes_gs[Seurat_DE_genes_gs$gs_weight<0,])%>% mutate(rank = rank(-1*gs_weight,ties.method ='random')) %>%
    arrange(rank) %>% mutate(rank = rank *(-1))
  Seurat_DE_genes_gs<-bind_rows(list(pos_rank,neg_rank))
  Seurat_DE_genes_gs<-Seurat_DE_genes_gs[,c('gene_id','rank')]
  
  Seurat_DE_genes_gs<-Seurat_DE_genes_gs %>% arrange(desc(rank))
  dup_ix_check<-which(duplicated(Seurat_DE_genes_gs$rank))
  print(paste('duplication check:',length(dup_ix_check),sep=''))
  
  Seurat_DE_genes_gs<-deframe(Seurat_DE_genes_gs)
  
  CP_GSEA<-clusterProfiler::GSEA(Seurat_DE_genes_gs, TERM2GENE =msig_database_filt,pvalueCutoff = 0.5,maxGSSize = 1000)#eps=1e-10
  #Res_out<-as.data.frame(fgseaRes[fgseaRes$padj<0.25,])
  Res_out<-CP_GSEA[order(CP_GSEA$p.adjust,decreasing = F), ]
  
  #return(Res_out)
  return(CP_GSEA)
}
msigdb_filt<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/GSEA/Compiled data/Gene set db/MSIGDB_ZS_gs_filtered_long_complete_gene_conv2022.rds')

GSEA.beat.mtwt.qlf.obj<-cprofiler_GSEA(DEG.beat.mtwt.qlf,msigdb_filt)
GSEA.beat.mtlwt.qlf.obj<-cprofiler_GSEA(DEG.beat.mtlwt.qlf,msigdb_filt)

# GSEA.tcga.mtwt.qlf.obj<-cprofiler_GSEA(DEG.tcga.mtwt.qlf,msigdb_filt)
# GSEA.tcga.mtlwt.qlf.obj<-cprofiler_GSEA(DEG.tcga.mtlwt.qlf,msigdb_filt)


GSEA.beat.mtwt.qlf.obj
GSEA.beat.mtlwt.qlf.obj
library(enrichplot)
library(ggplot2)
# GSEA_mtwt_BEAT_list<-c('TANG_SENESCENCE_TP53_TARGETS_DN','HALLMARK_TNFA_SIGNALING_VIA_NFKB',
# 'ZHOU_INFLAMMATORY_RESPONSE_LPS_UP','PRC2_EZH2_UP.V1_DN','ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_UP',
# 'BOQUEST_STEM_CELL_UP','HALLMARK_INFLAMMATORY_RESPONSE','PRC2_EZH2_UP.V1_UP','KONDO_EZH2_TARGETS',
# 'Top HSC','IVANOVA_HEMATOPOIESIS_STEM_CELL','IVANOVA_HEMATOPOIESIS_STEM_CELL_LONG_TERM','KEGG_OXIDATIVE_PHOSPHORYLATION',
# 'GOBP_MITOCHONDRIAL_CYTOCHROME_C_OXIDASE_ASSEMBLY','REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT',
# 'GOBP_MITOCHONDRIAL_ATP_SYNTHESIS_COUPLED_PROTON_TRANSPORT','GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY',
# 'WP_MITOCHONDRIAL_COMPLEX_I_ASSEMBLY_MODEL_OXPHOS_SYSTEM','REACTOME_MITOCHONDRIAL_TRANSLATION',
# 'HALLMARK_OXIDATIVE_PHOSPHORYLATION')
GSEA_mtwt_BEAT_list<-c('HALLMARK_TNFA_SIGNALING_VIA_NFKB',
                       'ZHOU_INFLAMMATORY_RESPONSE_LPS_UP',
                       'BOQUEST_STEM_CELL_UP','HALLMARK_INFLAMMATORY_RESPONSE','PRC2_EZH2_UP.V1_UP','KONDO_EZH2_TARGETS',
                       'Top HSC','IVANOVA_HEMATOPOIESIS_STEM_CELL_LONG_TERM','KEGG_OXIDATIVE_PHOSPHORYLATION',
                       'GOBP_MITOCHONDRIAL_CYTOCHROME_C_OXIDASE_ASSEMBLY','REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT',
                       'GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY',
                       'REACTOME_MITOCHONDRIAL_TRANSLATION',
                       'HALLMARK_OXIDATIVE_PHOSPHORYLATION')
GSEA.beat.mtwt.qlf.sub<-GSEA.beat.mtwt.qlf[GSEA_mtwt_BEAT_list,]
GSEA.beat.mtwt.qlf.sub$FDR_score<-(-log10(GSEA.beat.mtwt.qlf.sub$p.adjust)*sign(GSEA.beat.mtwt.qlf.sub$NES))


y<-mutate(GSEA.beat.mtwt.qlf.sub, ordering=abs(NES)) %>% arrange(desc(ordering)) 
mtwt.gsea.res<-ggplot(y, aes(NES, fct_reorder(Description, NES), fill = FDR_score), showCategory=(n*2)) + 
  geom_bar(stat='identity') + 
  scale_fill_gradient2(low='blue',mid='white',high='red')+ ylab(NULL)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"),
  legend.key=element_rect(fill="white"),
  axis.title = element_text(size=14),title = element_text(size=14),
  axis.text = element_text(size=14),
  legend.text=element_text(size=14))+  #theme(legend.text=element_text(size=18))
  guides(colour = guide_legend(override.aes = list(size=3)))

#2versions for MTlvsWT
GSEA_mtlwt_BEAT_list_v01<-c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_UP',
                            'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE',
                            'HALLMARK_INTERFERON_ALPHA_RESPONSE','HINATA_NFKB_TARGETS_FIBROBLAST_UP',
                            'HINATA_NFKB_TARGETS_KERATINOCYTE_UP','SANA_TNF_SIGNALING_UP',
                            'HALLMARK_IL6_JAK_STAT3_SIGNALING','P53_DN.V1_UP','Top HSC',
                            'GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR','REACTOME_INTERLEUKIN_6_SIGNALING',
                            'P53_DN.V2_UP','KONDO_EZH2_TARGETS','HSC II','IVANOVA_HEMATOPOIESIS_STEM_CELL',
                            'BOQUEST_STEM_CELL_UP','EPPERT_HSC_R','MOHANKUMAR_HOXA1_TARGETS_UP',
                            'HALLMARK_GLYCOLYSIS','BIOCARTA_MITOCHONDRIA_PATHWAY',
                            'REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_DNA_REPAIR_GENES',
                            'TP53_proteosome_205sig_dn','GOBP_FATTY_ACID_BETA_OXIDATION_USING_ACYL_COA_DEHYDROGENASE',
                            'LU_EZH2_TARGETS_UP','DANG_MYC_TARGETS_UP','GOBP_MITOCHONDRIAL_DNA_METABOLIC_PROCESS',
                            'GOBP_MITOCHONDRIAL_TRANSCRIPTION','KEGG_OXIDATIVE_PHOSPHORYLATION',
                            'GOBP_OXIDATIVE_PHOSPHORYLATION','HALLMARK_MYC_TARGETS_V1','HALLMARK_MYC_TARGETS_V2',
                            'GOBP_REGULATION_OF_MITOCHONDRIAL_GENE_EXPRESSION','HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                            'GOBP_MITOCHONDRIAL_TRANSLATION','GOBP_MITOCHONDRIAL_GENE_EXPRESSION',
                            'REACTOME_MITOCHONDRIAL_TRANSLATION')
#without MYC
GSEA_mtlwt_BEAT_list_v02<-c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_UP',
                            'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE',
                            'HALLMARK_INTERFERON_ALPHA_RESPONSE','HINATA_NFKB_TARGETS_FIBROBLAST_UP',
                            'HINATA_NFKB_TARGETS_KERATINOCYTE_UP','SANA_TNF_SIGNALING_UP',
                            'HALLMARK_IL6_JAK_STAT3_SIGNALING','P53_DN.V1_UP','Top HSC',
                            'GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR','REACTOME_INTERLEUKIN_6_SIGNALING',
                            'P53_DN.V2_UP','KONDO_EZH2_TARGETS','HSC II','IVANOVA_HEMATOPOIESIS_STEM_CELL',
                            'BOQUEST_STEM_CELL_UP','EPPERT_HSC_R','MOHANKUMAR_HOXA1_TARGETS_UP',
                            'HALLMARK_GLYCOLYSIS','BIOCARTA_MITOCHONDRIA_PATHWAY',
                            'REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_DNA_REPAIR_GENES',
                            'TP53_proteosome_205sig_dn','GOBP_FATTY_ACID_BETA_OXIDATION_USING_ACYL_COA_DEHYDROGENASE',
                            'LU_EZH2_TARGETS_UP','GOBP_MITOCHONDRIAL_DNA_METABOLIC_PROCESS',
                            'GOBP_MITOCHONDRIAL_TRANSCRIPTION','KEGG_OXIDATIVE_PHOSPHORYLATION',
                            'GOBP_OXIDATIVE_PHOSPHORYLATION',
                            'GOBP_REGULATION_OF_MITOCHONDRIAL_GENE_EXPRESSION','HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                            'GOBP_MITOCHONDRIAL_TRANSLATION','GOBP_MITOCHONDRIAL_GENE_EXPRESSION',
                            'REACTOME_MITOCHONDRIAL_TRANSLATION')
#without MYC and HOXA1 and Glycolysis (I needed to reduce list so I can fit to the paper)
GSEA_mtlwt_BEAT_list_v03<-c('HALLMARK_TNFA_SIGNALING_VIA_NFKB','ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_UP',
                            'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE',
                            'HALLMARK_INTERFERON_ALPHA_RESPONSE','HINATA_NFKB_TARGETS_FIBROBLAST_UP',
                            'HINATA_NFKB_TARGETS_KERATINOCYTE_UP','SANA_TNF_SIGNALING_UP',
                            'HALLMARK_IL6_JAK_STAT3_SIGNALING','P53_DN.V1_UP','Top HSC',
                            'GOBP_RESPONSE_TO_TUMOR_NECROSIS_FACTOR','REACTOME_INTERLEUKIN_6_SIGNALING',
                            'P53_DN.V2_UP','KONDO_EZH2_TARGETS','HSC II','IVANOVA_HEMATOPOIESIS_STEM_CELL',
                            'BOQUEST_STEM_CELL_UP','EPPERT_HSC_R','BIOCARTA_MITOCHONDRIA_PATHWAY',
                            'REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_DNA_REPAIR_GENES',
                            'TP53_proteosome_205sig_dn','GOBP_FATTY_ACID_BETA_OXIDATION_USING_ACYL_COA_DEHYDROGENASE',
                            'LU_EZH2_TARGETS_UP','GOBP_MITOCHONDRIAL_DNA_METABOLIC_PROCESS',
                            'GOBP_MITOCHONDRIAL_TRANSCRIPTION','KEGG_OXIDATIVE_PHOSPHORYLATION',
                            'GOBP_OXIDATIVE_PHOSPHORYLATION',
                            'GOBP_REGULATION_OF_MITOCHONDRIAL_GENE_EXPRESSION','HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                            'GOBP_MITOCHONDRIAL_TRANSLATION','GOBP_MITOCHONDRIAL_GENE_EXPRESSION',
                            'REACTOME_MITOCHONDRIAL_TRANSLATION')
#2023.01.28: subset the list for the visualization
GSEA_mtlwt_BEAT_list_v04<-c('HALLMARK_TNFA_SIGNALING_VIA_NFKB',
                            'HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_INTERFERON_GAMMA_RESPONSE',
                            'SANA_TNF_SIGNALING_UP',
                            'HALLMARK_IL6_JAK_STAT3_SIGNALING','Top HSC',
                            'KONDO_EZH2_TARGETS','HSC II','IVANOVA_HEMATOPOIESIS_STEM_CELL',
                            'BOQUEST_STEM_CELL_UP','EPPERT_HSC_R','BIOCARTA_MITOCHONDRIA_PATHWAY',
                            'GOBP_FATTY_ACID_BETA_OXIDATION_USING_ACYL_COA_DEHYDROGENASE',
                            'GOBP_MITOCHONDRIAL_TRANSCRIPTION','KEGG_OXIDATIVE_PHOSPHORYLATION',
                            'GOBP_REGULATION_OF_MITOCHONDRIAL_GENE_EXPRESSION','HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                            'GOBP_MITOCHONDRIAL_TRANSLATION','GOBP_MITOCHONDRIAL_GENE_EXPRESSION'
                            )
GSEA.beat.mtlwt.qlf.sub01<-GSEA.beat.mtlwt.qlf[GSEA_mtlwt_BEAT_list_v01,]
GSEA.beat.mtlwt.qlf.sub02<-GSEA.beat.mtlwt.qlf[GSEA_mtlwt_BEAT_list_v02,]
GSEA.beat.mtlwt.qlf.sub03<-GSEA.beat.mtlwt.qlf[GSEA_mtlwt_BEAT_list_v03,]
GSEA.beat.mtlwt.qlf.sub04<-GSEA.beat.mtlwt.qlf[GSEA_mtlwt_BEAT_list_v04,]


GSEA.beat.mtlwt.qlf.sub01$FDR_score<-(-log10(GSEA.beat.mtlwt.qlf.sub01$p.adjust)*sign(GSEA.beat.mtlwt.qlf.sub01$NES))
GSEA.beat.mtlwt.qlf.sub02$FDR_score<-(-log10(GSEA.beat.mtlwt.qlf.sub02$p.adjust)*sign(GSEA.beat.mtlwt.qlf.sub02$NES))
GSEA.beat.mtlwt.qlf.sub03$FDR_score<-(-log10(GSEA.beat.mtlwt.qlf.sub03$p.adjust)*sign(GSEA.beat.mtlwt.qlf.sub03$NES))
GSEA.beat.mtlwt.qlf.sub04$FDR_score<-(-log10(GSEA.beat.mtlwt.qlf.sub04$p.adjust)*sign(GSEA.beat.mtlwt.qlf.sub04$NES))

y01<-mutate(GSEA.beat.mtlwt.qlf.sub01, ordering=abs(NES)) %>% arrange(desc(ordering)) 
mtlwt.gsea.res01<-ggplot(y01, aes(NES, fct_reorder(Description, NES), fill = FDR_score)) + 
  geom_bar(stat='identity') + 
  #scale_fill_continuous(low='blue', high='red', guide=guide_colorbar(reverse=F)) + 
  scale_fill_gradient2(low='blue',mid='white',high='red')+ ylab(NULL)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),
        axis.title = element_text(size=14),title = element_text(size=14),
        axis.text = element_text(size=14),
        legend.text=element_text(size=14))+  #theme(legend.text=element_text(size=18))
  guides(colour = guide_legend(override.aes = list(size=3)))

y02<-mutate(GSEA.beat.mtlwt.qlf.sub02, ordering=abs(NES)) %>% arrange(desc(ordering)) 
mtlwt.gsea.res02<-ggplot(y02, aes(NES, fct_reorder(Description, NES), fill = FDR_score)) + 
  geom_bar(stat='identity') + 
  scale_fill_gradient2(low='blue',mid='white',high='red')+ylab(NULL)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),
        axis.title = element_text(size=14),title = element_text(size=14),
        axis.text = element_text(size=14),
        legend.text=element_text(size=14))+  #theme(legend.text=element_text(size=18))
  guides(colour = guide_legend(override.aes = list(size=3)))


y03<-mutate(GSEA.beat.mtlwt.qlf.sub03, ordering=abs(NES)) %>% arrange(desc(ordering)) 
mtlwt.gsea.res03<-ggplot(y03, aes(NES, fct_reorder(Description, NES), fill = FDR_score)) + 
  geom_bar(stat='identity') + 
  scale_fill_gradient2(low='blue',mid='white',high='red')+ylab(NULL)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),
        axis.title = element_text(size=14),title = element_text(size=14),
        axis.text = element_text(size=14),
        legend.text=element_text(size=14))+  #theme(legend.text=element_text(size=18))
  guides(colour = guide_legend(override.aes = list(size=3)))

y04<-mutate(GSEA.beat.mtlwt.qlf.sub04, ordering=abs(NES)) %>% arrange(desc(ordering)) 
mtlwt.gsea.res04<-ggplot(y04, aes(NES, fct_reorder(Description, NES), fill = FDR_score)) + 
  geom_bar(stat='identity') + 
  #scale_fill_continuous(low='blue', high='red', guide=guide_colorbar(reverse=F)) + 
  scale_fill_gradient2(low='blue',mid='white',high='red')+ ylab(NULL)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_rect(fill="white"),
        axis.title = element_text(size=14),title = element_text(size=14),
        axis.text = element_text(size=14),
        legend.text=element_text(size=14))+  #theme(legend.text=element_text(size=18))
  guides(colour = guide_legend(override.aes = list(size=3)))

# texts are too long - take y-axis label and gsea separaely
#+theme(axis.text.y=element_blank())
jpeg('Analysis/GEP/DEG_GSEA/GSEA/annotated ver/ZS_annotated_Dec022022/GSEA_figs/BEAT_MTWT_gsea_plot_20230222.jpg',height = 7, width = 10, units = 'in', res=300)
print(mtwt.gsea.res)
dev.off()
jpeg('Analysis/GEP/DEG_GSEA/GSEA/annotated ver/ZS_annotated_Dec022022/GSEA_figs/BEAT_MTWT_gsea_ytext_20230222.jpg',height = 7, width = 10, units = 'in', res=300)
print(mtwt.gsea.res+theme(axis.text.y=element_blank()))
dev.off()
jpeg('Analysis/GEP/DEG_GSEA/GSEA/annotated ver/ZS_annotated_Dec022022/GSEA_figs/BEAT_MTLWT_gsea01_plot.jpg',height = 7, width = 10, units = 'in', res=300)
print(mtlwt.gsea.res01)
dev.off()
jpeg('Analysis/GEP/DEG_GSEA/GSEA/annotated ver/ZS_annotated_Dec022022/GSEA_figs/BEAT_MTLWT_gsea01_ytext.jpg',height = 7, width = 10, units = 'in', res=300)
print(mtlwt.gsea.res01+theme(axis.text.y=element_blank()))
dev.off()
jpeg('Analysis/GEP/DEG_GSEA/GSEA/annotated ver/ZS_annotated_Dec022022/GSEA_figs/BEAT_MTLWT_gsea02_plot.jpg',height = 7, width = 10, units = 'in', res=300)
print(mtlwt.gsea.res02)
dev.off()
jpeg('Analysis/GEP/DEG_GSEA/GSEA/annotated ver/ZS_annotated_Dec022022/GSEA_figs/BEAT_MTLWT_gsea02_ytext.jpg',height = 7, width = 10, units = 'in', res=300)
print(mtlwt.gsea.res02+theme(axis.text.y=element_blank()))
dev.off()

jpeg('Analysis/GEP/DEG_GSEA/GSEA/annotated ver/ZS_annotated_Dec022022/GSEA_figs/BEAT_MTLWT_gsea03_plot.jpg',height = 7, width = 10, units = 'in', res=300)
print(mtlwt.gsea.res03)
dev.off()
jpeg('Analysis/GEP/DEG_GSEA/GSEA/annotated ver/ZS_annotated_Dec022022/GSEA_figs/BEAT_MTLWT_gsea03_ytext.jpg',height = 7, width = 10, units = 'in', res=300)
print(mtlwt.gsea.res03+theme(axis.text.y=element_blank()))
dev.off()


jpeg('Analysis/GEP/DEG_GSEA/GSEA/annotated ver/ZS_annotated_Dec022022/GSEA_figs/BEAT_MTLWT_gsea04_plot_20230222.jpg',height = 7, width = 10, units = 'in', res=300)
print(mtlwt.gsea.res04)
dev.off()
jpeg('Analysis/GEP/DEG_GSEA/GSEA/annotated ver/ZS_annotated_Dec022022/GSEA_figs/BEAT_MTLWT_gsea04_ytext_20230222.jpg',height = 7, width = 10, units = 'in', res=300)
print(mtlwt.gsea.res04+theme(axis.text.y=element_blank()))
dev.off()

#print out every individual GSEA plots.

GSEA.beat.mtwt.qlf.obj
GSEA.beat.mtlwt.qlf.obj

GSEA_mtwt_BEAT_list
GSEA_mtlwt_BEAT_list_v01
library(enrichplot)
GSEA_BEAT_mtwt_plot<-list()
GSEA_BEAT_mtlwt_plot<-list()

for(i in 1:length(GSEA_mtwt_BEAT_list)){
  GSEA_BEAT_mtwt_plot[[i]]<-gseaplot2(GSEA.beat.mtwt.qlf.obj,geneSetID =GSEA_mtwt_BEAT_list[i] ,
                                      # title=TP53mut_psebulk_GSEA$originalDescription['BIOCARTA_PROTEASOME_PATHWAY'],
                                      title=GSEA_mtwt_BEAT_list[i],base_size = 16,
                                      pvalue_table=F)
}
for(i in 1:length(GSEA_mtlwt_BEAT_list_v01)){
  GSEA_BEAT_mtlwt_plot[[i]]<-gseaplot2(GSEA.beat.mtlwt.qlf.obj,geneSetID =GSEA_mtlwt_BEAT_list_v01[i] ,
                                      # title=TP53mut_psebulk_GSEA$originalDescription['BIOCARTA_PROTEASOME_PATHWAY'],
                                      title=GSEA_mtlwt_BEAT_list_v01[i],base_size = 16,
                                      pvalue_table=F)
}
for(i in 1:length(GSEA_mtwt_BEAT_list)){
  
  jpeg(paste('Analysis/GEP/DEG_GSEA/GSEA/annotated ver/ZS_annotated_Dec022022/GSEA_figs/BEAT_MTWT/BEAT_',GSEA_mtwt_BEAT_list[i],'.jpg',sep=''),height = 7, width = 10, units = 'in', res=300)
  print(GSEA_BEAT_mtwt_plot[[i]])
  dev.off()
}
for(i in 1:length(GSEA_mtlwt_BEAT_list_v01)){
  
  jpeg(paste('Analysis/GEP/DEG_GSEA/GSEA/annotated ver/ZS_annotated_Dec022022/GSEA_figs/BEAT_MTLWT/BEAT_',GSEA_mtlwt_BEAT_list_v01[i],'.jpg',sep=''),height = 7, width = 10, units = 'in', res=300)
  print(GSEA_BEAT_mtlwt_plot[[i]])
  dev.off()
}

####################################################################################















#####################################################################################

#y.beat.mtwt.cpm<-as.data.frame(cpm(y.beat.mtwt,log=F,normalized.lib.sizes = TRUE))
#y.tcga.cpm<-as.data.frame(cpm(y.tcga,log=F,normalized.lib.sizes = TRUE))
#y.comb.cpm<-as.data.frame(log2(cpm(y.comb,log=F,normalized.lib.sizes = TRUE)+1))

#perform DEG analysis

#make group annotation
BEAT_cls.mtwt<-rep('WT',dim(y.beat.mtwt)[2]);BEAT_cls.mtwt[colnames(y.beat.mtwt)%in%BEAT_MUT_id]<-'TP53_MUT';BEAT_cls.mtwt[colnames(y.beat.mtwt)%in%BEAT_MUT_like_id]<-'TP53_MUTlike'
TCGA_cls.mtwt<-rep('WT',dim(y.tcga)[2]);TCGA_cls.mtwt[gsub('.03A.*|.03B.*','',colnames(y.tcga))%in%TCGA_MUT_id]<-'TP53_MUT';TCGA_cls.mtwt[gsub('.03A.*|.03B.*','',colnames(y.tcga))%in%TCGA_MUT_like_id]<-'TP53_MUTlike'
#data to use
y.beat.mtwt.n<-y.beat.mtwt[,c(mt_ix,wt_ix)]
y.beat.wt<-y.beat.mtwt[,c(mtl_ix,wt_ix)]
y.tcga.mtwt<-y.tcga[,c(mt_ix.tcga,wt_ix.tcga)]
y.tcga.wt<-y.tcga[,c(mtl_ix.tcga,wt_ix.tcga)]
y.comb.mtwt<-y.comb[,c(colnames(y.beat.mtwt),colnames(y.tcga))]

#all.equal(colnames(y.comb.mtwt),c(colnames(y.beat.mtwt),colnames(y.tcga)))

BEAT_cls.mtwt.n<-BEAT_cls.mtwt[c(mt_ix,wt_ix)]
BEAT_cls.mtlwt<-BEAT_cls.mtwt[c(mtl_ix,wt_ix)]
TCGA_cls.mtwt.n<-TCGA_cls.mtwt[c(mt_ix.tcga,wt_ix.tcga)]
TCGA_cls.wt<-TCGA_cls.mtwt[c(mtl_ix.tcga,wt_ix.tcga)]
mtl_ix.comb<-which(gsub('.03A.*|.03B.*','',colnames(y.comb.mtwt))%in% c(BEAT_MUT_like_id,TCGA_MUT_like_id))
wt_ix.comb<-which(!(gsub('.03A.*|.03B.*','',colnames(y.comb.mtwt))%in% c(BEAT_MUT_id,TCGA_MUT_id,BEAT_MUT_like_id,TCGA_MUT_like_id)))
Comb_cls<-rep('WT',dim(y.comb.mtwt)[2]);Comb_cls[gsub('.03A.*|.03B.*','',colnames(y.comb.mtwt))%in% c(BEAT_MUT_id,TCGA_MUT_id)]<-'TP53_MUT';Comb_cls[gsub('.03A.*|.03B.*','',colnames(y.comb.mtwt))%in% c(BEAT_MUT_like_id,TCGA_MUT_like_id)]<-'TP53_MUTlike'
Comb_cls.wt<-Comb_cls[c(mtl_ix.comb,wt_ix.comb)]
y.comb.wt<-y.comb.mtwt[,c(mtl_ix.comb,wt_ix.comb)]

#y.comb
beat.mtwt.group<-factor(ifelse(BEAT_cls.mtwt.n=='TP53_MUT',1,0))
beat.mtlwt.group<-factor(ifelse(BEAT_cls.mtlwt=='TP53_MUTlike',1,0))

tcga.mtwt.group<-factor(ifelse(TCGA_cls.mtwt.n=='TP53_MUT',1,0))
tcga.mtlwt.group<-factor(ifelse(TCGA_cls.wt=='TP53_MUTlike',1,0))


design.beat.mtwt<-model.matrix(~beat.mtwt.group)
design.beat.mtlwt<-model.matrix(~beat.mtlwt.group)
design.tcga.mtwt<-model.matrix(~tcga.mtwt.group)
design.tcga.mtlwt<-model.matrix(~tcga.mtlwt.group)

y.beat.mtwt.n<-estimateDisp(y.beat.mtwt.n,design = design.beat.mtwt)
y.beat.mtlwt<-estimateDisp(y.beat.wt,design = design.beat.mtlwt)
y.tcga.mtwt<-estimateDisp(y.tcga.mtwt,design = design.tcga.mtwt)
y.tcga.mtlwt<-estimateDisp(y.tcga.wt,design = design.tcga.mtlwt)

# plotBCV(y.beat.mtwt.n)
# plotBCV(y.beat.mtlwt)
# plotBCV(y.tcga.mtwt)
# plotBCV(y.tcga.mtlwt)
# plotBCV(y.comb.mtwt)
# plotBCV(y.comb.mtlwt)


fit.beat.mtwt<-glmQLFit(y.beat.mtwt.n,design = design.beat.mtwt,robust=T)
fit.beat.mtlwt<-glmQLFit(y.beat.mtlwt,design = design.beat.mtlwt,robust=T)
fit.tcga.mtwt<-glmQLFit(y.tcga.mtwt,design = design.tcga.mtwt,robust=T)
fit.tcga.mtlwt<-glmQLFit(y.tcga.mtlwt,design = design.tcga.mtlwt,robust=T)

qlf.beat.mtwt.qlf <- glmQLFTest(fit.beat.mtwt,coef=2)
qlf.beat.mtwt.treat <- glmTreat(fit.beat.mtwt,coef=2,lfc=1)
qlf.beat.mtlwt.qlf <- glmQLFTest(fit.beat.mtlwt,coef=2)
qlf.beat.mtlwt.treat <- glmTreat(fit.beat.mtlwt,coef=2,lfc=1)

qlf.tcga.mtwt.qlf <- glmQLFTest(fit.tcga.mtwt,coef=2)
qlf.tcga.mtwt.treat <- glmTreat(fit.tcga.mtwt,coef=2,lfc=1)
qlf.tcga.mtlwt.qlf <- glmQLFTest(fit.tcga.mtlwt,coef=2)
qlf.tcga.mtlwt.treat <- glmTreat(fit.tcga.mtlwt,coef=2,lfc=1)

summary(decideTestsDGE(qlf.beat.mtwt.qlf))
summary(decideTestsDGE(qlf.beat.mtwt.treat))
summary(decideTestsDGE(qlf.beat.mtlwt.qlf))
summary(decideTestsDGE(qlf.beat.mtlwt.treat))

summary(decideTestsDGE(qlf.tcga.mtwt.qlf))
summary(decideTestsDGE(qlf.tcga.mtwt.treat))
summary(decideTestsDGE(qlf.tcga.mtlwt.qlf))
summary(decideTestsDGE(qlf.tcga.mtlwt.treat))

#combined DGE

comb.mtlwt.group<-factor(ifelse(Comb_cls.wt=='TP53_MUTlike',1,0))
design.comb.mtlwt<-model.matrix(~comb.mtlwt.group)
y.comb.mtlwt<-estimateDisp(y.comb.wt,design = design.comb.mtlwt)
plotBCV(y.comb.mtlwt)
fit.comb.mtlwt<-glmQLFit(y.comb.mtlwt,design = design.comb.mtlwt,robust=T)
qlf.comb.mtlwt.qlf <- glmQLFTest(fit.comb.mtlwt,coef=2)
qlf.comb.mtlwt.treat <- glmTreat(fit.comb.mtlwt,coef=2,lfc=1)

# View(topTags(qlf.beat.mtwt.qlf,Inf)$table)

DEG.beat.mtwt.qlf<-(topTags(qlf.beat.mtwt.qlf,Inf)$table)
DEG.beat.mtwt.treat<-(topTags(qlf.beat.mtwt.treat,Inf)$table)
DEG.beat.mtlwt.qlf<-(topTags(qlf.beat.mtlwt.qlf,Inf)$table)
DEG.beat.mtlwt.treat<-(topTags(qlf.beat.mtlwt.treat,Inf)$table)

DEG.tcga.mtwt.qlf<-(topTags(qlf.tcga.mtwt.qlf,Inf)$table)
DEG.tcga.mtwt.treat<-(topTags(qlf.tcga.mtwt.treat,Inf)$table)
DEG.tcga.mtlwt.qlf<-(topTags(qlf.tcga.mtlwt.qlf,Inf)$table)
DEG.tcga.mtlwt.treat<-(topTags(qlf.tcga.mtlwt.treat,Inf)$table)

DEG.comb.mtlwt.qlf<-(topTags(qlf.comb.mtlwt.qlf,Inf)$table)
DEG.comb.mtlwt.treat<-(topTags(qlf.comb.mtlwt.treat,Inf)$table)
# 
# write.csv(DEG.beat.mtwt.qlf,'DEG_GSEA/DEG/Reanalysis_aug02/BEAT_MTWT_DEG_qlf_re.csv')
# write.csv(DEG.beat.mtwt.treat,'DEG_GSEA/DEG/Reanalysis_aug02/BEAT_MTWT_DEG_treat_re.csv')
# write.csv(DEG.beat.mtlwt.qlf,'DEG_GSEA/DEG/Reanalysis_aug02/BEAT_MTlikeWT_DEG_qlf_re.csv')
# write.csv(DEG.beat.mtlwt.treat,'DEG_GSEA/DEG/Reanalysis_aug02/BEAT_MTlikeWT_DEG_treat_re.csv')
# 
# write.csv(DEG.tcga.mtwt.qlf,'DEG_GSEA/DEG/Reanalysis_aug02/TCGA_MTWT_DEG_qlf_re.csv')
# write.csv(DEG.tcga.mtwt.treat,'DEG_GSEA/DEG/Reanalysis_aug02/TCGA_MTWT_DEG_treat_re.csv')
# write.csv(DEG.tcga.mtlwt.qlf,'DEG_GSEA/DEG/Reanalysis_aug02/TCGA_MTlikeWT_DEG_qlf_re.csv')
# write.csv(DEG.tcga.mtlwt.treat,'DEG_GSEA/DEG/Reanalysis_aug02/TCGA_MTlikeWT_DEG_treat_re.csv')
# 
# write.csv(DEG.comb.mtlwt.qlf,'DEG_GSEA/DEG/Reanalysis_aug02/Comb_MTlikeWT_DEG_qlf_re.csv')
# write.csv(DEG.comb.mtlwt.treat,'DEG_GSEA/DEG/Reanalysis_aug02/Comb_MTlikeWT_DEG_treat_re.csv')


shared_DEG_updn_sig<-function(DEG_01.df,DEG_02.df,pval_cutoff=0.05){
  #pval_cutoff
  DEG01<-DEG_01.df
  DEG02<-DEG_02.df
  OL_DEG<-intersect(rownames(DEG01)[DEG01$FDR<pval_cutoff],
                    rownames(DEG02)[DEG02$FDR<pval_cutoff])
  shared.all<-data.frame('DEG'=OL_DEG,'DEG01_FC'=rep(NA,length(OL_DEG)),'DEG02_FC'=rep(NA,length(OL_DEG)),
                         'DEG01_FDR'=rep(NA,length(OL_DEG)),'DEG02_FDR'=rep(NA,length(OL_DEG)))
  for(i in 1:length(shared.all$DEG)){
    DEG01_ix<-which(rownames(DEG01) %in% shared.all$DEG[i])
    DEG02_ix<-which(rownames(DEG02) %in% shared.all$DEG[i])
    shared.all$DEG01_FC[i]<-DEG01$logFC[DEG01_ix]
    shared.all$DEG02_FC[i]<-DEG02$logFC[DEG02_ix]
    shared.all$DEG01_FDR[i]<-DEG01$FDR[DEG01_ix]
    shared.all$DEG02_FDR[i]<-DEG02$FDR[DEG02_ix]
  }
  shared.all.clean<-shared.all[c(which(shared.all$DEG01_FC>0 & shared.all$DEG02_FC>0),
                                 which(shared.all$DEG01_FC<0 & shared.all$DEG02_FC<0)),]
  return(shared.all.clean)
}



DEG.BEAT.TCGA.mtwt.shared<-shared_DEG_updn_sig(DEG.beat.mtwt.qlf,DEG.tcga.mtwt.qlf,0.05)
DEG.BEAT.TCGA.mtlwt.shared<-shared_DEG_updn_sig(DEG.beat.mtlwt.qlf,DEG.tcga.mtlwt.qlf,0.05)
# write.csv(DEG.BEAT.TCGA.mtwt.shared,'DEG_GSEA/DEG/Reanalysis_aug02/BEAT_TCGA_OL_MTWT_DEG_qlf_re.csv')
# write.csv(DEG.BEAT.TCGA.mtlwt.shared,'DEG_GSEA/DEG/Reanalysis_aug02/BEAT_TCGA_OL_MTWT_DEG_qlf_re.csv')

#load('DEG_GSEA/TP53GEP_DEG_GSEA_reanalysis.Rdata')
#rm(y.beat.cpm,y.tcga.cpm)
#genes that are only in beat TP53mut-like and not in TP53WT

BEAT_MUT_id<-BEAT_clinical_final$LabId[BEAT_clinical_final$cls.GLM=='TP53_MUT']
BEAT_MUT_like_id<-BEAT_clinical_final$LabId[BEAT_clinical_final$cls.GLM=='TP53_WTunk_MT_like']#actually just WT_MTlike, no unknown added
TCGA_MUT_id<-TCGA_clinical_final$patient_RNA_id[TCGA_clinical_final$cls.GLM=='TP53_MUT']
TCGA_MUT_like_id<-TCGA_clinical_final$patient_RNA_id[TCGA_clinical_final$cls.GLM=='TP53_unk_MT_like']

#need to have 2 diferent data: MT vs. WT and MTlike vs. WTlike
#data to use
####
# y.beat.mtwt<-y.beat[,TP53MT_WT_ix]
y.beat.mtwt.ordered<-y.beat[,c(mt_ix,mtl_ix,wt_ix)]
# y.beat.mtwt.n<-y.beat.mtwt[,c(mt_ix,wt_ix)]
# y.beat.wt<-y.beat.mtwt[,c(mtl_ix,wt_ix)]
# y.tcga.mtwt<-y.tcga[,c(mt_ix.tcga,wt_ix.tcga)]
# y.tcga.wt<-y.tcga[,c(mtl_ix.tcga,wt_ix.tcga)]

#y.comb.mtwt<-y.comb[,c(colnames(y.beat.mtwt),colnames(y.tcga))]

####
# get index of each group (MUT, MUtlike and WT)
#beat index
mt_ix01<-which(colnames(y.beat.mtwt) %in% BEAT_MUT_id)
mtl_ix01<-which(colnames(y.beat.mtwt) %in% BEAT_MUT_like_id)
wt_ix01<-setdiff(1:dim(y.beat.mtwt)[2],c(mt_ix01,mtl_ix01))

#tcga index
mt_ix.tcga01<-which(gsub('.03A.*|.03B.*','',colnames(y.tcga)) %in% TCGA_MUT_id)
mtl_ix.tcga01<-which(gsub('.03A.*|.03B.*','',colnames(y.tcga)) %in% TCGA_MUT_like_id)
wt_ix.tcga01<-setdiff(1:dim(y.tcga)[2],c(mt_ix.tcga01,mtl_ix.tcga01))

#make a vector for colnames of comb data
comb_colid<-colnames(y.comb.mtwt)
comb_colid[414:591]<-gsub('.03A.*|.03B.*','',comb_colid[414:591])
#comb data index
mt.ix.comb01<-which((comb_colid)
                    %in%c(BEAT_MUT_id,TCGA_MUT_id))
mtl_ix.comb01<-which((comb_colid)
                     %in%c(BEAT_MUT_like_id,TCGA_MUT_like_id))
wt_ix.comb01<-setdiff(1:length(comb_colid),c(mt.ix.comb01,mtl_ix.comb01))



mtl_gene_ct<-rowSums(y.beat.mtwt[,mtl_ix]$counts>=1)
wt_gene_ct<-rowSums(y.beat.mtwt[,wt_ix]$counts<1) #number of patients that have 0 counts (counting TRUE)
beat.cutoff<-0.65
final_ix<-which(mtl_gene_ct>(length(mtl_ix)*beat.cutoff) & wt_gene_ct>(length(wt_ix)*beat.cutoff))
View(rownames_to_column(as.data.frame(cbind(mtl_gene_ct[final_ix],wt_gene_ct[final_ix]))))
cand_genes<-rownames_to_column(as.data.frame(cbind(mtl_gene_ct[final_ix],wt_gene_ct[final_ix])))$rowname
#View(rownames_to_column(DEG.beat.mtlwt.treat,var='Symbol'))
View(rownames_to_column(DEG.beat.mtlwt.qlf[cand_genes,]))
#cand_genes.beat.dge<-rownames_to_column(DEG.beat.mtlwt.treat[cand_genes,]) %>% filter(FDR<0.05)
cand_genes.beat.dge<-rownames_to_column(DEG.beat.mtlwt.qlf[cand_genes,]) %>% filter(FDR<0.05)

#are these genes shared and statistically significant in both BEAT and TCGA?
rownames(DEG.BEAT.TCGA.mtlwt.shared)<-DEG.BEAT.TCGA.mtlwt.shared$DEG
cand_genes.beat.dge.shared<-DEG.BEAT.TCGA.mtlwt.shared[cand_genes,] %>% filter(DEG01_FDR<0.05 &DEG02_FDR<0.05)

write.csv(cand_genes.beat.dge,'DEG_GSEA/DEG/Marker genes/BEAT_DEG_mtl_wt_065.csv')
write.csv(cand_genes.beat.dge.shared,'DEG_GSEA/DEG/Marker genes/BEAT_065_TCGA_shared_DEG_mtl_wt.csv')



#TCGA
#TCGA_MUT_id<-TCGA_clinical_final$patient_RNA_id[TCGA_clinical_final$cls.GLM=='TP53_MUT']
#TCGA_MUT_like_id<-TCGA_clinical_final$patient_RNA_id[TCGA_clinical_final$cls.GLM=='TP53_unk_MT_like']

mtl_gene_ct.tcga<-rowSums(y.tcga[,mtl_ix.tcga01]$counts>=1)
wt_gene_ct.tcga<-rowSums(y.tcga[,wt_ix.tcga01]$counts<1)

final_ix.tcga<-which(mtl_gene_ct.tcga>(length(mtl_ix.tcga01)*0.65) & wt_gene_ct.tcga>(length(wt_ix.tcga01)*0.65))
View(rownames_to_column(as.data.frame(cbind(mtl_gene_ct.tcga[final_ix.tcga],wt_gene_ct.tcga[final_ix.tcga]))))
cand_genes.tcga<-rownames_to_column(as.data.frame(cbind(mtl_gene_ct.tcga[final_ix.tcga],wt_gene_ct.tcga[final_ix.tcga])))$rowname
#View(rownames_to_column(DEG.beat.mtlwt.treat,var='Symbol'))
View(rownames_to_column(DEG.tcga.mtlwt.qlf[cand_genes.tcga,]))

intersect(rownames(DEG.tcga.mtlwt.qlf[cand_genes.tcga,]),rownames(DEG.beat.mtlwt.qlf[cand_genes,]))

#try for combined data
#dim(y.comb)
dim(y.comb.mtwt)
# 
# mt_ix.comb<-c(which(colnames(y.comb.mtwt) %in% BEAT_MUT_id),which(gsub('.03A.*|.03B.*','',colnames(y.comb.mtwt)) %in% TCGA_MUT_id))
# mtl_ix.comb<-c(which(colnames(y.comb.mtwt) %in% BEAT_MUT_like_id),which(gsub('.03A.*|.03B.*','',colnames(y.comb.mtwt)) %in% TCGA_MUT_like_id))
# wt_ix.comb<-setdiff(1:dim(y.comb.mtwt)[2],c(mt_ix.comb,mtl_ix.comb))

mtl_gene_ct.comb<-rowSums(y.comb.mtwt[,mtl_ix.comb01]$counts>=1)
wt_gene_ct.comb<-rowSums(y.comb.mtwt[,wt_ix.comb01]$counts<1)

cutoff_comb<-0.60
final_ix.comb<-which(mtl_gene_ct.comb>(length(mtl_ix.comb01)*cutoff_comb) & wt_gene_ct.comb>(length(wt_ix.comb01)*cutoff_comb))
View(rownames_to_column(as.data.frame(cbind(mtl_gene_ct.comb[final_ix.comb],wt_gene_ct.comb[final_ix.comb]))))
cand_genes.comb<-rownames_to_column(as.data.frame(cbind(mtl_gene_ct.comb[final_ix.comb],wt_gene_ct.comb[final_ix.comb])))$rowname

# View(rownames_to_column(DEG.beat.mtlwt.qlf[cand_genes.comb,]))
# View(rownames_to_column(DEG.tcga.mtlwt.qlf[cand_genes.comb,]))
# View(DEG.BEAT.TCGA.mtlwt.shared[cand_genes.comb,])

cand_genes.comb.dge<-(DEG.BEAT.TCGA.mtlwt.shared[cand_genes.comb,]) %>% filter(DEG01_FDR<0.05 & DEG02_FDR<0.05)
write.csv(cand_genes.comb.dge,'DEG_GSEA/DEG/Marker genes/Comb_shared_DEG_marker_mtl_wt_060.csv')
#save.image('DEG_GSEA/DEG/TP53Mutlike_Marker_genes.RData')
load('DEG_GSEA/DEG/TP53Mutlike_Marker_genes.RData')
gc()
#generate heatmap for the marker genes
###
mt_ix.comb<-c(which(colnames(y.comb.mtwt) %in% BEAT_MUT_id),which(gsub('.03A.*|.03B.*','',colnames(y.comb.mtwt)) %in% TCGA_MUT_id))
mtl_ix.comb<-c(which(colnames(y.comb.mtwt) %in% BEAT_MUT_like_id),which(gsub('.03A.*|.03B.*','',colnames(y.comb.mtwt)) %in% TCGA_MUT_like_id))
wt_ix.comb<-setdiff(1:dim(y.comb.mtwt)[2],c(mt_ix.comb,mtl_ix.comb))

y.beat.mtwt.ordered<-y.beat.mtwt[,c(mt_ix,mtl_ix,wt_ix)]
y.tcga.ordered<-y.tcga[,c(mt_ix.tcga,mtl_ix.tcga,wt_ix.tcga)]
y.comb.mtwt.ordered<-y.comb.mtwt[,c(mt_ix.comb,mtl_ix.comb,wt_ix.comb)]

###

# DEG.beat.mtwt.qlf
# DEG.beat.mtlwt.qlf
# DEG.tcga.mtwt.qlf
# DEG.tcga.mtlwt.qlf
# 
# DEG.comb.mtwt.qlf
# DEG.comb.mtlwt.qlf

cand_genes.comb.dge
y.beat.mtwt.ord.cpm<-as.data.frame(cpm(y.beat.mtwt.ordered,log=T,normalized.lib.sizes = TRUE))
#y.beat.mtwt.cpm<-as.data.frame(cpm(y.beat.mtwt,log=T,normalized.lib.sizes = TRUE))
y.tcga.ordered.cpm<-as.data.frame(cpm(y.tcga.ordered,log=T,normalized.lib.sizes = TRUE))

#c(mt_ix,mtl_ix,wt_ix)
#y.beat.mtwt.cpm.dge<-(y.beat.mtwt.cpm[cand_genes.comb.dge$DEG,])#comb DEG
# y.beat.mtwt.ord.cpm.dge<-(y.beat.mtwt.ord.cpm[cand_genes.comb.dge$DEG,])#comb DEG
# y.tcga.ord.cpm.dge<-(y.tcga.ordered.cpm[cand_genes.comb.dge$DEG,])#comb DEG


BEAT_clinical_final.sub<-BEAT_clinical_final[!BEAT_clinical_final$cls.GLM=='TP53_unk_WT_like',]
BEAT_clinical_final.sub<-BEAT_clinical_final.sub[,c('LabId','cls.GLM')]
anno.df<-data.frame('LabId'=colnames(y.beat.mtwt.ord.cpm.dge))
anno.df<-right_join(anno.df,BEAT_clinical_final.sub,by='LabId') #%>%t()
anno.df<-anno.df %>% column_to_rownames(var='LabId')
anno.df$cls.GLM[which(anno.df$cls.GLM=='TP53_WTunk_MT_like')]<-'TP53_MUT_like'

anno_beat_dge<-HeatmapAnnotation(df=anno.df,col = list(cls.GLM=c('TP53_MUT'='red','TP53_WT'='blue','TP53_MUT_like'='purple')))

beat_hm_marker<-Heatmap(scale_GEP(y.beat.mtwt.ord.cpm.dge), 
                        show_row_names = T,show_column_names = F,top_annotation = anno_beat_dge,
                        column_names_max_height = unit(4, "cm"),
                        column_names_gp = gpar(fontsize = 6),
                        row_names_max_width = unit(6, "cm"),
                        row_names_gp = gpar(fontsize = 12),
                        show_row_dend = F,cluster_rows = F,cluster_columns =F)



TCGA_clinical_final.sub<-TCGA_clinical_final
TCGA_clinical_final.sub<-TCGA_clinical_final.sub[,c('patient_RNA_id','cls.GLM')]
anno.tcga.df<-data.frame('patient_RNA_id'=gsub('.03A.*|.03B.*','',colnames(y.tcga.ord.cpm.dge)))
anno.tcga.df<-right_join(anno.tcga.df,TCGA_clinical_final.sub,by='patient_RNA_id') #%>%t()
anno.tcga.df<-anno.tcga.df %>% column_to_rownames(var='patient_RNA_id')
anno.tcga.df$cls.GLM[which(anno.tcga.df$cls.GLM=='TP53_unk_MT_like')]<-'TP53_MUT_like'

anno_tcga_dge<-HeatmapAnnotation(df=anno.tcga.df,col = list(cls.GLM=c('TP53_MUT'='red','TP53_WT'='blue','TP53_MUT_like'='purple')))

tcga_hm_marker<-Heatmap(scale_GEP(y.tcga.ord.cpm.dge), 
                        show_row_names = T,show_column_names = F,top_annotation = anno_tcga_dge,
                        column_names_max_height = unit(4, "cm"),
                        column_names_gp = gpar(fontsize = 6),
                        row_names_max_width = unit(6, "cm"),
                        row_names_gp = gpar(fontsize = 12),
                        show_row_dend = F,cluster_rows = F,cluster_columns =F)
jpeg('DEG_GSEA/DEG/Marker genes/BEAT_marker_10_logcpmz.jpg',height = 7, width = 10, units = 'in', res=300)
print(beat_hm_marker)
dev.off()
jpeg('DEG_GSEA/DEG/Marker genes/TCGA_marker_10_logcpmz.jpg',height = 7, width = 10, units = 'in', res=300)
print(tcga_hm_marker)
dev.off()

#do binary version: check the Zohar's idea on Box (All Files> YoonKyu> p53 mutant Like Project> Mut Like Gene List Sample Data presentation)
#looking at broader gene now
#ZS: Would like to see if we can find a group of genes where 1 (or more) of them is expressed 
#in all/most mut-like samples and none/few wt samples (meaning, they are not co-expressed 
#but all/most samples express at least one of them). See spreadsheet for ideas on expressing 
#this data. Use only 1 or 0 to express expression versus not expressed. Will look at magnitude later.
BEAT.mtwt.bin<-as.data.frame(ifelse(y.beat.mtwt.ordered$counts>0,1,0))
TCGA.mtwt.bin<-as.data.frame(ifelse(y.tcga.ordered$counts>0,1,0))
Comb.mtwt.bin<-as.data.frame(ifelse(y.comb.mtwt.ordered$counts>0,1,0))

#beat index
mt_ix02<-which(colnames(BEAT.mtwt.bin) %in% BEAT_MUT_id)
mtl_ix02<-which(colnames(BEAT.mtwt.bin) %in% BEAT_MUT_like_id)
wt_ix02<-setdiff(1:dim(BEAT.mtwt.bin)[2],c(mt_ix02,mtl_ix02))

#tcga index
mt_ix.tcga02<-which(gsub('.03A.*|.03B.*','',colnames(TCGA.mtwt.bin)) %in% TCGA_MUT_id)
mtl_ix.tcga02<-which(gsub('.03A.*|.03B.*','',colnames(TCGA.mtwt.bin)) %in% TCGA_MUT_like_id)
wt_ix.tcga02<-setdiff(1:dim(TCGA.mtwt.bin)[2],c(mt_ix.tcga02,mtl_ix.tcga02))

mt_ix.comb.bin<-c(which(colnames(Comb.mtwt.bin) %in% BEAT_MUT_id),which(gsub('.03A.*|.03B.*','',colnames(Comb.mtwt.bin)) %in% TCGA_MUT_id))
mtl_ix.comb.bin<-c(which(colnames(Comb.mtwt.bin) %in% BEAT_MUT_like_id),which(gsub('.03A.*|.03B.*','',colnames(Comb.mtwt.bin)) %in% TCGA_MUT_like_id))
wt_ix.comb.bin<-setdiff(1:dim(Comb.mtwt.bin)[2],c(mt_ix.comb.bin,mtl_ix.comb.bin))

#beat aml count
mtl_gene_ct.beat.bin<-rowSums(BEAT.mtwt.bin[,mtl_ix02]>=1)
wt_gene_ct.beat.bin<-rowSums(BEAT.mtwt.bin[,wt_ix02]<1)
sample_colsums.beat<-colSums(BEAT.mtwt.bin)#genes simple filtered: 43494
write.csv(sample_colsums.beat,'beat_colsums.csv')
#tcga aml count
mtl_gene_ct.tcga.bin<-rowSums(TCGA.mtwt.bin[,mtl_ix.tcga02]>=1)
wt_gene_ct.tcga.bin<-rowSums(TCGA.mtwt.bin[,wt_ix.tcga02]<1)
sample_colsums.tcga<-colSums(TCGA.mtwt.bin)#genes simple filtered: 20K~
write.csv(sample_colsums.tcga,'tcga_colsums.csv')

#comb aml count
mtl_gene_ct.comb.bin<-rowSums(Comb.mtwt.bin[,mtl_ix02]>=1)
wt_gene_ct.comb.bin<-rowSums(Comb.mtwt.bin[,wt_ix02]<1)
sample_colsums.comb<-colSums(Comb.mtwt.bin)#genes simple filtered: 20K~

#
cutoff_beat.bin.mtl<-0.52
cutoff_beat.bin.wt<-0.80
beat_mtl_mgene<-which(mtl_gene_ct.beat.bin>=(length(mtl_ix02)*cutoff_beat.bin.mtl))
#mtl genes need another condition that it   
beat_wt_mgene<-which(wt_gene_ct.beat.bin>=(length(wt_ix02)*cutoff_beat.bin.wt))
beat_final_ix.bin<-intersect(beat_mtl_mgene,beat_wt_mgene)
target_genes<-rownames(BEAT.mtwt.bin)[beat_final_ix.bin]
Heatmap((BEAT.mtwt.bin[beat_final_ix.bin,]), 
        show_row_names = T,show_column_names = F,top_annotation = anno_beat_dge,
        column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 12),
        show_row_dend = F,cluster_rows = F,cluster_columns =F)

#check in DEG
a<-(DEG.beat.mtlwt.qlf[target_genes,]%>%filter(FDR<0.05))
write.csv(a,'DEG_GSEA/DEG/BEAT_detailed_markerDEG.csv')
#tcga count
mtl_gene_ct.tcga.bin<-rowSums(TCGA.mtwt.bin[,mtl_ix.tcga02]>=1)
wt_gene_ct.tcga.bin<-rowSums(TCGA.mtwt.bin[,wt_ix.tcga02]<1)
sample_colsums.tcga<-colSums(TCGA.mtwt.bin)#OL genes simple filtered: 22651

cutoff_tcga.bin.mtl<-0.52
cutoff_tcga.bin.wt<-0.80
tcga_mtl_mgene<-which(mtl_gene_ct.tcga.bin>=(length(mtl_ix.tcga02)*cutoff_tcga.bin.mtl))
#mtl genes need another condition that it   
tcga_wt_mgene<-which(wt_gene_ct.tcga.bin>=(length(wt_ix.tcga02)*cutoff_tcga.bin.wt))
tcga_final_ix.bin<-intersect(tcga_mtl_mgene,tcga_wt_mgene)
target_genes.tcga<-rownames(TCGA.mtwt.bin)[tcga_final_ix.bin]
Heatmap((TCGA.mtwt.bin[tcga_final_ix.bin,]), 
        show_row_names = T,show_column_names = F,top_annotation = anno_tcga_dge,
        column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 12),
        show_row_dend = F,cluster_rows = F,cluster_columns =F)

#check in DEG
a<-(DEG.tcga.mtlwt.qlf[target_genes.tcga,]%>%filter(FDR<0.05))
dim(a)
write.csv(a,'DEG_GSEA/DEG/TCGA_detailed_markerDEG.csv')


#combined list count
mtl_gene_ct.comb.bin<-rowSums(Comb.mtwt.bin[,mtl_ix.comb.bin]>=1)
wt_gene_ct.comb.bin<-rowSums(Comb.mtwt.bin[,wt_ix.comb.bin]<1)
sample_colsums.comb<-colSums(Comb.mtwt.bin)#OL genes simple filtered: 21289

cutoff_comb.bin.mtl<-0.5
cutoff_comb.bin.wt<-0.8
#
comb_mtl_mgene<-which(mtl_gene_ct.comb.bin>=(length(mtl_ix.comb.bin)*cutoff_comb.bin.mtl))
#mtl genes need another condition that it   
comb_wt_mgene<-which(wt_gene_ct.comb.bin>=(length(wt_ix.comb.bin)*cutoff_comb.bin.wt))
comb_final_ix.bin<-intersect(comb_mtl_mgene,comb_wt_mgene)
target_genes.comb<-rownames(Comb.mtwt.bin)[comb_final_ix.bin]

comb.cls<-rep(NA,dim(Comb.mtwt.bin)[2]);comb.cls[mt_ix.comb.bin]<-'TP53_MUT';comb.cls[mtl_ix.comb.bin]<-'TP53_MUT_like';comb.cls[wt_ix.comb.bin]<-'TP53_WT'
anno_comb_dge<-data.frame('LabID'=colnames(Comb.mtwt.bin),'cls.GLM'=comb.cls) %>% column_to_rownames(var='LabID')
anno_comb_dge<-HeatmapAnnotation(df=anno_comb_dge,col = list(cls.GLM=c('TP53_MUT'='red','TP53_WT'='blue','TP53_MUT_like'='purple')))

Heatmap((Comb.mtwt.bin[comb_final_ix.bin,]), 
        show_row_names = T,show_column_names = F,top_annotation = anno_comb_dge,
        column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 12),
        show_row_dend = F,cluster_rows = F,cluster_columns =F)

#check in DEG
a<-(DEG.BEAT.TCGA.mtlwt.shared[target_genes.comb,]%>%filter(DEG01_FDR<0.05 &DEG02_FDR<0.05))
write.csv(a,'DEG_GSEA/DEG/COMB_detailed_markerDEG.csv')

#check by rank high in WT first()



#comb data: using heat map TCGA has higher expressions in general

aa<-(cbind(mtl_gene_ct.beat.bin,wt_gene_ct.beat.bin)[which(mtl_gene_ct.beat.bin>=(length(mtl_ix02)*0.3) & wt_gene_ct.beat.bin>=length(wt_ix02)*0.85),])
aaa<-(cbind(aa,DEG.beat.mtlwt.qlf[rownames(aa),c(4,5)]))
#write.csv(aaa,'DEG_GSEA/DEG/BEAT_detailed_markerDEG_rank.csv')
beat_dge_hm_03_085<-Heatmap(as.matrix(BEAT.mtwt.bin[rownames(aaa),]), 
                            show_row_names = T,show_column_names = F,top_annotation = anno_beat_dge,
                            column_names_max_height = unit(3, "cm"),
                            column_names_gp = gpar(fontsize = 6),
                            row_names_max_width = unit(2, "cm"),
                            row_names_gp = gpar(fontsize = 6),
                            show_row_dend = F,cluster_rows = F,cluster_columns =F)

bb<-(cbind(mtl_gene_ct.tcga.bin,wt_gene_ct.tcga.bin)[which(mtl_gene_ct.tcga.bin>=(length(mtl_ix.tcga02)*0.3) & wt_gene_ct.tcga.bin>=length(wt_ix.tcga02)*0.85),])
bbb<-(cbind(bb,DEG.tcga.mtlwt.qlf[rownames(bb),c(4,5)]))
#write.csv(aaa,'DEG_GSEA/DEG/TCGA_detailed_markerDEG_rank.csv')
tcga_dge_hm_03_085<-Heatmap((TCGA.mtwt.bin[rownames(bbb),]), 
                            show_row_names = T,show_column_names = F,top_annotation = anno_tcga_dge,
                            column_names_max_height = unit(3, "cm"),
                            column_names_gp = gpar(fontsize = 6),
                            row_names_max_width = unit(3, "cm"),
                            row_names_gp = gpar(fontsize = 8),
                            show_row_dend = F,cluster_rows = F,cluster_columns =F)

intersect(rownames(aaa),rownames(bbb))

cc<-(cbind(mtl_gene_ct.comb.bin,wt_gene_ct.comb.bin)[which(mtl_gene_ct.comb.bin>=length(mtl_ix.comb.bin)*0.3 & wt_gene_ct.comb.bin>=length(wt_ix.comb.bin)*0.85),])
#View(aa)
#write.csv(aa,'DEG_GSEA/DEG/COMB_detailed_markerDEG_rank.csv')
ccc<-DEG.BEAT.TCGA.mtlwt.shared[rownames(cc),]%>%filter(DEG01_FDR<0.05&DEG02_FDR<0.05)
comb_dge_hm_03_085<-Heatmap((Comb.mtwt.bin[rownames(ccc),]), 
                            show_row_names = T,show_column_names = F,top_annotation = anno_comb_dge,
                            column_names_max_height = unit(4, "cm"),
                            column_names_gp = gpar(fontsize = 6),
                            row_names_max_width = unit(3, "cm"),
                            row_names_gp = gpar(fontsize = 6),
                            show_row_dend = F,cluster_rows = F,cluster_columns =F)
View(DEG.BEAT.TCGA.mtlwt.shared[rownames(aa),])
View(cbind(aa,DEG.BEAT.TCGA.mtlwt.shared[rownames(aa),4]))

jpeg('DEG_GSEA/DEG/Marker genes/beat_dge_hm_030_085.jpg',height = 7, width = 10, units = 'in', res=300)
print(beat_dge_hm_03_085)
dev.off()
jpeg('DEG_GSEA/DEG/Marker genes/tcga_dge_hm_030_085.jpg',height = 7, width = 10, units = 'in', res=300)
print(tcga_dge_hm_03_085)
dev.off()
jpeg('DEG_GSEA/DEG/Marker genes/comb_dge_hm_030_085.jpg',height = 7, width = 10, units = 'in', res=300)
print(comb_dge_hm_03_085)
dev.off()



###################################################
##########################################################################\######################

#do the auc figure

library(readxl)
library(dplyr)
library(reshape2)
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP")

BEAT_cyto_clean<-readRDS('class ids/BEAT_cyto_clean.rds')
BEAT_wt_id<-BEAT_cyto_clean$LabId
BEAT_Drug_AUC<-read_excel("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/BEAT_AML_drug_sensitivity_AUC.xlsx")
#high auc resist and low auc sensitive; so transform the score to high=sensitive
BEAT_Drug_AUC.scale<-BEAT_Drug_AUC %>% group_by(inhibitor)%>%mutate(auc_scale=scale(auc)*(-1))


patiend_id<-unique(BEAT_Drug_AUC.scale$lab_id)
BEAT_Drug_AUC.rev<-BEAT_Drug_AUC.scale %>% filter(lab_id %in% BEAT_wt_id)

##intersect(BEAT_wt_id, BEAT_Drug_AUC.rev$lab_id)
BEAT_Drug_AUC.rev$cls.GLM<-rep(0,dim(BEAT_Drug_AUC.rev)[1])
BEAT_Drug_AUC.rev$cls.GLM[which(BEAT_Drug_AUC.rev$lab_id %in% BEAT_cyto_clean$LabId[BEAT_cyto_clean$cls.GLM.x=='TP53_MUT'])]<-'TP53_MUT'
BEAT_Drug_AUC.rev$cls.GLM[which(BEAT_Drug_AUC.rev$lab_id %in% BEAT_cyto_clean$LabId[BEAT_cyto_clean$cls.GLM.x=='TP53_WT'])]<-'TP53_WT'
BEAT_Drug_AUC.rev$cls.GLM[which(BEAT_Drug_AUC.rev$lab_id %in% BEAT_cyto_clean$LabId[BEAT_cyto_clean$cls.GLM.x=='TP53MUT_like'])]<-'TP53MUT_like'

venetoclax.df<-BEAT_Drug_AUC.rev[BEAT_Drug_AUC.rev$inhibitor=='Venetoclax',]

veneto.sum<-data_summary(venetoclax.df,varname='auc_scale',groupnames=c('cls.GLM'))
# veneto.p<-ggplot(veneto.sum,aes(x=cls.GLM,y=auc_scale,group=cls.GLM)) +
#   geom_bar(stat='identity',fill='steelblue',width=0.5)+theme_minimal()+
#   geom_errorbar( aes(x=cls.GLM, ymin=auc_scale-se, ymax=auc_scale+se), width=0.1, colour="orange", alpha=0.9, size=0.3)
veneto.p<-ggplot(veneto.sum,aes(x=cls.GLM,y=auc_scale,group=cls.GLM)) +
  geom_bar(stat='identity',fill='steelblue',width=0.5)+theme_minimal()+
  geom_errorbar( aes(x=cls.GLM, ymin=auc_scale-se, ymax=auc_scale+se), width=0.1, colour="orange", alpha=0.9, size=0.3)

veneto.p.scale<-ggdotplot(venetoclax.df,x='cls.GLM',y='auc_scale',add = 'mean_se',color='black',fill='black',size=0.1,
                          position=position_jitter(width=0.125,height=0.4),binwidth = 1,alpha=0.5,order=c('TP53_MUT','TP53MUT_like','TP53_WT'))+
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=15)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53MUT_like','TP53_WT')))+
  font("xlab", size = 18)+
  font("ylab", size = 18)




veneto.p<-ggdotplot(venetoclax.df,x='cls.GLM',y='auc',add = 'mean_se',color='black',fill='black',size=3,
                    position=position_jitter(0.3),binwidth = 1.5,alpha=0.5,order=c('TP53_MUT','TP53MUT_like','TP53_WT'))+
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=15)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53MUT_like','TP53_WT')))+
  font("xlab", size = 18)+
  font("ylab", size = 18)

stat.veneto <- venetoclax.df %>%
  #group_by(cls.GLM) %>%
  wilcox_test(auc_scale~cls.GLM) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
p.adjust(stat.veneto[c(1,3),]$p,'BH')

jpeg('auc_veneto_scale.jpeg',width=1024,height=640)
veneto.p.scale
dev.off()

jpeg('auc_veneto_auc.jpeg',width=1024,height=640)
veneto.p
dev.off()



