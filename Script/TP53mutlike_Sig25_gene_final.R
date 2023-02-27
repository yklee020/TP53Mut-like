#2022.12.13
#TP53mutlike_Sig25_gene_final.r
#Yoonkyu Lee
#Sachs Lab

##################################
#create visualization for sig_gene_25

source("C:/Users/yklee/Desktop/Sachs Lab/R_project/TP53/TP53GEP_functions.r")
packages<-c('edgeR','limma','dplyr','tidyr','sva','tidyverse',
            'ggplot2','glmnet','caret','pROC','survival','survminer','HGNChelper')

for(p in packages){
  library(p,character.only = T)
}

ridge_model.fin<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/mtl_sig_gene/Sig_gene_ridge_all_final_model.rds")
Sig_gene<-predict(ridge_model.fin,type='coef')@Dimnames[[1]][-1]

BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')

#below GEP data has genes were transformed into same reference
BEAT_RNA_counts<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_counts_mtmtlwt.rds")
TCGA_AML_counts<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/TCGA_AML_counts_mtmtlwt.rds")
#BEAT_TCGA_merge.corrected<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_TCGA_OL_data/BEAT461_TCGA_combat_count.rds")


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
y.beat.cpm<-as.data.frame(cpm(y.beat,log=T,normalized.lib.sizes = TRUE))
y.tcga.cpm<-as.data.frame(cpm(y.tcga,log=T,normalized.lib.sizes = TRUE))

all.equal(colnames(y.beat),BEAT_clinical_final$LabId)#TRUE
all.equal(colnames(y.tcga),TCGA_clinical_final$patient_RNA_id)
# TP53WT_ix<-which(BEAT_clinical_final$TP53_mut_stat != 'unknown' & BEAT_clinical_final$TP53_mut_stat != 'TP53_MUT')
# TP53WT_ix.tcga<-which(TCGA_clinical_final$TP53_mut_stat != 'TP53_MUT')
TP53WT_id<-which(BEAT_clinical_final$TP53_mut_stat != 'unknown' & BEAT_clinical_final$TP53_mut_stat != 'TP53_MUT')


#need to have 2 diferent data: MT vs. WT and MTlike vs. WTlike
#y.beat.mtwt<-y.beat[,c(BEAT_MUT_id,BEAT_WT_id)]
y.beat.cpm.mtlwt<-y.beat.cpm[,c(BEAT_MUT_like_id,BEAT_WT_id)]

#y.tcga
#y.tcga.mtwt<-y.tcga[,c(TCGA_MUT_id,TCGA_WT_id)]
y.tcga.cpm.mtlwt<-y.tcga.cpm[,c(TCGA_MUT_like_id,TCGA_WT_id)]

coef(ridge_model.fin)
#BEAT_AML and TCGA shows CREBBP,REC8, UGCG, OASL, THBD and PEAR1 as up regulated Put them first up and down
Sig_gene
library(ComplexHeatmap)
mt_ix03<-which(colnames(y.beat.cpm) %in% BEAT_MUT_id)
mtl_ix03<-which(colnames(y.beat.cpm) %in% BEAT_MUT_like_id)
wt_ix03<-setdiff(1:dim(y.beat.cpm)[2],c(mt_ix03,mtl_ix03))
cls.beat<-rep('TP53_WT',dim(y.beat.cpm)[2]);cls.beat[mt_ix03]<-'TP53_MUT';cls.beat[mtl_ix03]<-'TP53MUT_like'
cls.beat.df<-data.frame('ID'=colnames(y.beat.cpm),'cls'= cls.beat ) %>% column_to_rownames(var='ID')
anno_beat_gep<-HeatmapAnnotation(df=cls.beat.df,
                                 col = list(cls=c('TP53MUT_like'='purple','TP53_WT'='blue1')),simple_anno_size = unit(0.3, "cm"))
beat.sig25<-Heatmap(t(scale(t(y.beat.cpm[Sig_gene,c(BEAT_MUT_like_id,BEAT_WT_id)]))), 
                    show_row_names = T,show_column_names = F,top_annotation = anno_beat_gep[c(mtl_ix03,wt_ix03)],
                    column_names_max_height = unit(4, "cm"),
                    column_names_gp = gpar(fontsize = 6),
                    row_names_max_width = unit(6, "cm"),
                    row_names_gp = gpar(fontsize = 12),
                    show_row_dend = F,cluster_rows = T,cluster_columns =F)
TCGA_MUT_id
mt_ix03.tcga<-which(gsub('.03A.*|.03B.*','',colnames(y.tcga.cpm)) %in% TCGA_MUT_id)
mtl_ix03.tcga<-which(gsub('.03A.*|.03B.*','',colnames(y.tcga.cpm)) %in% TCGA_MUT_like_id)
wt_ix03.tcga<-setdiff(1:dim(y.tcga.cpm)[2],c(mt_ix03.tcga,mtl_ix03.tcga))
cls.tcga<-rep('TP53_WT',dim(y.tcga.cpm)[2]);cls.tcga[mt_ix03.tcga]<-'TP53_MUT';cls.tcga[mtl_ix03.tcga]<-'TP53MUT_like'
cls.tcga.df<-data.frame('ID'=colnames(y.tcga.cpm),'cls'= cls.tcga ) %>% column_to_rownames(var='ID')
anno_tcga_gep<-HeatmapAnnotation(df=cls.tcga.df,
                                 col = list(cls=c('TP53MUT_like'='purple','TP53_WT'='blue1')),simple_anno_size = unit(0.3, "cm"))
tcga.sig25<-Heatmap(scale_GEP(y.tcga.cpm[Sig_gene,c(TCGA_MUT_like_id,TCGA_WT_id)]), 
                    show_row_names = T,show_column_names = F,top_annotation = anno_tcga_gep[c(mtl_ix03.tcga,wt_ix03.tcga)],
                    column_names_max_height = unit(4, "cm"),
                    column_names_gp = gpar(fontsize = 6),
                    row_names_max_width = unit(6, "cm"),
                    row_names_gp = gpar(fontsize = 12),
                    show_row_dend = F,cluster_rows = T,cluster_columns =F)

anno_beat_gep<-HeatmapAnnotation(df=cls.beat.df,
                                 col = list(cls=c('TP53MUT_like'='purple','TP53_WT'='blue1')),
                                 show_legend = c(F,F),show_annotation_name = F,simple_anno_size = unit(0.3, "cm"))
beat.sig25.nolab<-Heatmap(t(scale(t(y.beat.cpm[Sig_gene,c(BEAT_MUT_like_id,BEAT_WT_id)]))), 
                    show_row_names = F,show_column_names = F,top_annotation = anno_beat_gep[c(mtl_ix03,wt_ix03)],
                    column_names_max_height = unit(4, "cm"),
                    column_names_gp = gpar(fontsize = 6),
                    row_names_max_width = unit(6, "cm"),
                    row_names_gp = gpar(fontsize = 12),
                    show_row_dend = F,cluster_rows = T,cluster_columns =F)

anno_tcga_gep<-HeatmapAnnotation(df=cls.tcga.df,
                                 col = list(cls=c('TP53MUT_like'='purple','TP53_WT'='blue1')),
                                 show_legend = c(F,F),show_annotation_name = F,simple_anno_size = unit(0.3, "cm"))
tcga.sig25.nolab<-Heatmap(scale_GEP(y.tcga.cpm[Sig_gene,c(TCGA_MUT_like_id,TCGA_WT_id)]), 
                    show_row_names = F,show_column_names = F,top_annotation = anno_tcga_gep[c(mtl_ix03.tcga,wt_ix03.tcga)],
                    column_names_max_height = unit(4, "cm"),
                    column_names_gp = gpar(fontsize = 6),
                    row_names_max_width = unit(6, "cm"),
                    row_names_gp = gpar(fontsize = 12),
                    show_row_dend = F,cluster_rows = T,cluster_columns =F)


pdf("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Sig_gene_25/Sig_gene_25_heatmap.pdf")
beat.sig25
tcga.sig25
dev.off()

jpeg('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Sig_gene_25/BEAT_Sig_gene25.jpeg',width=6.5,height=5,units='in',res=300)
beat.sig25
dev.off()
jpeg('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Sig_gene_25/TCGA_Sig_gene25.jpeg',width=6.5,height=5,units='in',res=300)
tcga.sig25
dev.off()
jpeg('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Sig_gene_25/BEAT_Sig_gene25_nolabel.jpeg',width=3.5,height=3.75,units='in',res=300)
beat.sig25.nolab
dev.off()
jpeg('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Sig_gene_25/TCGA_Sig_gene25_nolabel.jpeg',width=3.5,height=3.75,units='in',res=300)
tcga.sig25.nolab
dev.off()