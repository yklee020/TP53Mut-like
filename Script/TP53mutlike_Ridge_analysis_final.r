#2022.12.1
#TP53mutlike_Ridge_analysis_final.r

#2022.12.12
#run regression models using just Overlap GEP in BEAT and TCGA

source("C:/Users/yklee/Desktop/Sachs Lab/R_project/TP53/TP53GEP_functions.r")
packages<-c('edgeR','limma','dplyr','tidyr','sva','tidyverse',
            'ggplot2','glmnet','caret','pROC','survival','survminer','reshape2')

for(p in packages){
  library(p,character.only = T)
}

################################################################################
BEAT_RNA_counts.461.clinic<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_AML/BEAT_AML_preprocessed/BEAT_RNA_counts_461_clinic.rds")
TCGA_AML_counts<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/TCGA/TCGA_preprocessed/TCGA_RNA_counts.rds")

#BEAT AML clinical data
#BEAT_ct_clinc.RNA_avail.all.01<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_AML/BEAT_AML_preprocessed/BEAT_ct_clinc_RNA_avail_all_01.rds')
#BEAT_ct_clinc.RNA_avail.all.01.WES<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_AML/BEAT_AML_preprocessed/BEAT_ct_clinc_RNA_avail_all_WES_01.rds')
BEAT_ct_clinc.RNA_avail.all.01.WES<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_AML/BEAT_AML_preprocessed/BEAT_ct_RNA_01.WES.allmut.IP.rds')
TCGA_clinical.rna.GDC.178<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/TCGA/TCGA_preprocessed/TCGA_clinical_GDC_rna_178_final.rds')

#overlap genes
BEAT_genes<-rownames(BEAT_RNA_counts.461.clinic)
TCGA_genes<-rownames(TCGA_AML_counts)
OL_genes<-intersect(BEAT_genes,TCGA_genes)
BEAT_TCGA_count<-cbind(BEAT_RNA_counts.461.clinic[OL_genes, ],
                       TCGA_AML_counts[OL_genes, ])

y.beat<-DGEList(counts=BEAT_RNA_counts.461.clinic)
y.beat<-calcNormFactors(y.beat)
y.beat.cpm<-as.data.frame(cpm(y.beat,log=FALSE,normalized.lib.sizes = TRUE))
y.tcga<-DGEList(counts=TCGA_AML_counts)
y.tcga<-calcNormFactors(y.tcga)
y.tcga.cpm<-as.data.frame(cpm(y.tcga,log=FALSE,normalized.lib.sizes = TRUE))


BEAT_AML_cpm.OL<-y.beat.cpm[OL_genes, ]
TCGA_AML_cpm.OL<-y.tcga.cpm[OL_genes, ]
BEAT_TCGA_cpm.OL<-cbind(BEAT_AML_cpm.OL,TCGA_AML_cpm.OL)
keep<-rowSums(BEAT_TCGA_cpm.OL>0)>=5
BEAT_TCGA_cpm.OL<-BEAT_TCGA_cpm.OL[keep,]
beat_tcga_group<-c(rep('BEAT',dim(BEAT_AML_cpm.OL)[2]),
                   rep('TCGA',dim(TCGA_AML_cpm.OL)[2]))
rm(y.beat,y.beat.cpm,y.tcga,y.tcga.cpm,BEAT_AML_cpm.OL,TCGA_AML_cpm.OL)
gc()

#try combat-seq to correct batch effects
#BEAT_TCGA_merge.corrected<-ComBat_seq(counts = as.matrix(BEAT_TCGA_count), batch = beat_tcga_group)
#saveRDS(BEAT_TCGA_merge.corrected,"C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_TCGA_OL_data/BEAT461_TCGA_combat_count.rds")
BEAT_TCGA_merge.corrected<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_TCGA_OL_data/BEAT461_TCGA_combat_count.rds")
y.comb<-DGEList(counts=BEAT_TCGA_merge.corrected)
y.comb<-calcNormFactors(y.comb)
#y.comb.cpm<-as.data.frame(cpm(y.comb,log=FALSE,normalized.lib.sizes = TRUE))
#y.comb.logcpm<-as.data.frame(cpm(y.comb,log=T,normalized.lib.sizes = TRUE))
y.comb.logcpm<-as.data.frame(log2(cpm(y.comb,log=FALSE,normalized.lib.sizes = TRUE)+1))
#rm(BEAT_TCGA_merge.corrected)
#saveRDS(y.comb.cpm,"C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_TCGA_OL_data/BEAT_TCGA_combat_cpm.rds")
beat_tcga_group<-c(rep('BEAT',461),
                   rep('TCGA',178))

# BEAT_ct_clinc.RNA_avail.all.01
# BEAT_ct_clinc.RNA_avail.all.01.WES
# TCGA_clinical.rna.GDC.178
#make a group annotaions for mutational subtypes


#RNAseq data name order and clinical id order are different - match them for the mutational status using ids
beat_mut_id01<-BEAT_ct_clinc.RNA_avail.all.01.WES$LabId[BEAT_ct_clinc.RNA_avail.all.01.WES$TP53_WES_class=='TP53_MUT']
beat_wt_id01<-BEAT_ct_clinc.RNA_avail.all.01.WES$LabId[BEAT_ct_clinc.RNA_avail.all.01.WES$TP53_WES_class=='WT']
beat_unk_id01<-BEAT_ct_clinc.RNA_avail.all.01.WES$LabId[BEAT_ct_clinc.RNA_avail.all.01.WES$TP53_WES_class=='unknown']

#tcga_mut_id<-TCGA_clinical.rna.prep.178$`Patient Identifier`[which(TCGA_clinical.rna.prep.178$TP53=='TP53_MUT')]
tcga_mut_id<-TCGA_clinical.rna.GDC.178$patient_RNA_id[which(TCGA_clinical.rna.GDC.178$TP53MUT_status==1)]
tcga_wt_id<-TCGA_clinical.rna.GDC.178$patient_RNA_id[which(TCGA_clinical.rna.GDC.178$TP53MUT_status==0)]

#integrate BEAT AML and TCGA mutational status
BEAT_TCGA_mut_stat<-rep(0,dim(BEAT_TCGA_merge.corrected)[2])

BEAT_tp53mutix<-which(colnames(BEAT_TCGA_merge.corrected)%in% beat_mut_id01)
BEAT_tp53wtix<-which(colnames(BEAT_TCGA_merge.corrected)%in% beat_wt_id01)
BEAT_tp53unkix<-which(colnames(BEAT_TCGA_merge.corrected)%in% beat_unk_id01)

TCGA_tp53mutix<-which(gsub('.03A.*|.03B.*','',colnames(BEAT_TCGA_merge.corrected))
                      %in% gsub('-','.',tcga_mut_id))
TCGA_tp53wtix<-which(gsub('.03A.*|.03B.*','',colnames(BEAT_TCGA_merge.corrected))
                     %in% gsub('-','.',tcga_wt_id))

BEAT_TCGA_mut_stat[BEAT_tp53mutix]<-'TP53_MUT'
BEAT_TCGA_mut_stat[BEAT_tp53wtix]<-'TP53_WT'
BEAT_TCGA_mut_stat[BEAT_tp53unkix]<-'unknown'

BEAT_TCGA_mut_stat[TCGA_tp53mutix]<-'TP53_MUT'
BEAT_TCGA_mut_stat[TCGA_tp53wtix]<-'TP53_WT'

setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis/LASSO/GEP")

#prep for the LASSO data
GEP.Data.corrected<-as.data.frame(t(y.comb.logcpm))
GEP.Data.corrected$sampleID <- rownames(GEP.Data.corrected)
GEP.Data.corrected$group<-beat_tcga_group
GEP.Data.corrected$TP53_mut_stat<-BEAT_TCGA_mut_stat

#Things to do:
#1. Train/test on BEAT AML and test on TCGA data.
##1.1. separate BEAT AML and TCGA from the BEAT_TCGA_comb.corrected PCA data for LASSO
##1.2. as a sanity check, run PCA for BEAT AML and project TCGA on the same PC space (with corrected data)
##1.3. find BEAT AML PCs that are different bw TP53_MUT vs. WT
##1.4. run LASSO train/test  


##1.1 separate BEAT AML and TCGA from the BEAT_TCGA_comb.corrected GEP data
#BEAT AML(n=461) and TCGA AML (n=178)
GEP.Data.corrected#1:642 are PCs
BEAT_GEP.Data.corrected<-GEP.Data.corrected[1:461,]
TCGA_GEP.Data.corrected<-GEP.Data.corrected[462:dim(GEP.Data.corrected)[1],]

TP53MT_WT_ix<-which(BEAT_GEP.Data.corrected$TP53_mut_stat != 'unknown')
TP53_MT_nMT<-rep(0,dim(BEAT_GEP.Data.corrected)[1])
TP53_MT_nMT[which(BEAT_GEP.Data.corrected$TP53_mut_stat=='TP53_MUT')]<-1
TP53_MT_WT<-TP53_MT_nMT[TP53MT_WT_ix]

#do 2ways - 1) all the OL genes and DE genes among OL genes 
#do DE analysis in BEAT data - use DE genes for the regression
#do 2 DE for BEAT AML TP53MUT vs. WT and TP53MUT vs. all the others

#y.comb
y.count.beat.mtwt<-y.comb[,TP53MT_WT_ix]
y.count.beat.mtnmt<-y.comb[,1:461]
y.count.tcga<-y.comb[,462:dim(y.comb)[2]]

beat.tp53_mtwt.group<-factor(TP53_MT_WT)
beat.tp53_mtnmt.group<-factor(ifelse(BEAT_GEP.Data.corrected$TP53_mut_stat=='TP53_MUT',1,0))
tcga.tp53_mtnmt.group<-factor(ifelse(TCGA_GEP.Data.corrected$TP53_mut_stat=='TP53_MUT',1,0))



##1.4. run LASSO train/test  

TP53_group_sub<-as.factor(BEAT_GEP.Data.corrected$TP53_mut_stat[TP53MT_WT_ix])
TP53_group_sub.y<-(ifelse(TP53_group_sub=='TP53_MUT',1,0))

#subset data for MT_WT and unknown and with all vs. DEGs
#1.all genes
BEAT_GEP_TP53MT_WT.all<-data.matrix(BEAT_GEP.Data.corrected[TP53MT_WT_ix,-(21619:21621)])
BEAT_GEP_TP53nMT.all<-data.matrix(BEAT_GEP.Data.corrected[-TP53MT_WT_ix,-(21619:21621)])

setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP")
# save.image('WES_save_001.Rdata')
# load('WES_save_001.Rdata')

set.seed(9011230)
#set.seed(90112312)
#set.seed(9011233)
trn.idx<-TP53_group_sub.y %>% createDataPartition(p=0.6, list = FALSE)
trn.data<-data.matrix(BEAT_GEP_TP53MT_WT.all[trn.idx,])
trn.y<-(TP53_group_sub.y[trn.idx])
test.data<-data.matrix(BEAT_GEP_TP53MT_WT.all[-trn.idx,])
test.y<-(TP53_group_sub.y[-trn.idx])

#set 10-fold cv for the model comparison
set.seed(284)
flds <- createFolds(trn.y, k = 10, list = TRUE, returnTrain = FALSE)
foldids = rep(1,length(trn.y))
foldids[flds$Fold02] = 2
foldids[flds$Fold03] = 3
foldids[flds$Fold04] = 4
foldids[flds$Fold05] = 5
foldids[flds$Fold06] = 6
foldids[flds$Fold07] = 7
foldids[flds$Fold08] = 8
foldids[flds$Fold09] = 9
foldids[flds$Fold10] = 10
#try ridge regression and elastic net 
cv.lasso<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=1,keep = TRUE,
                    type.measure = 'class',foldid=foldids)
cv.ridge<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=0,keep = TRUE,
                    type.measure = 'class',foldid=foldids)
cv.elastic<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=0.5,keep = TRUE,
                      type.measure = 'class',foldid=foldids)
# cv.lasso<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=1,keep = TRUE,
#                     type.measure = 'class')
# cv.ridge<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=0,keep = TRUE,
#                     type.measure = 'class')
# cv.elastic<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=0.5,keep = TRUE,
#                       type.measure = 'class')

plot(cv.lasso)
plot(cv.ridge)
plot(cv.elastic)
# model.lasso <- glmnet(trn.data, trn.y, alpha = 1, family = "binomial",
#                       lambda = cv.lasso$lambda.1se)
# model.ridge <- glmnet(trn.data, trn.y, alpha = 0, family = "binomial",
#                       lambda = cv.ridge$lambda.1se)
# model.elastic <- glmnet(trn.data, trn.y, alpha = 0.5, family = "binomial",
#                         lambda = cv.elastic$lambda.1se)

# 
model.lasso <- glmnet(trn.data, trn.y, alpha = 1, family = "binomial",
                      lambda = cv.lasso$lambda.min)
model.ridge <- glmnet(trn.data, trn.y, alpha = 0, family = "binomial",
                      lambda = cv.ridge$lambda.min)
model.elastic <- glmnet(trn.data, trn.y, alpha = 0.5, family = "binomial",
                        lambda = cv.elastic$lambda.min)

prob.lasso.res <- model.lasso %>% predict(newx = trn.data,s=cv.lasso$lambda.1se,type='response')
prob.ridge.res <- model.ridge %>% predict(newx = trn.data,s=cv.ridge$lambda.1se,type='response')
prob.elastic.res <- model.elastic %>% predict(newx = trn.data,s=cv.elastic$lambda.1se,type='response')


roc_trn.lasso<-pROC_summary(trn.y,prob.lasso.res)
PR_trn.lasso<-Prec_Recall_summary(trn.y,prob.lasso.res)

roc_trn.ridge<-pROC_summary(trn.y,prob.ridge.res)
PR_trn.ridge<-Prec_Recall_summary(trn.y,prob.ridge.res)

roc_trn.elastic<-pROC_summary(trn.y,prob.elastic.res)
PR_trn.elastic<-Prec_Recall_summary(trn.y,prob.elastic.res)


prob.lasso.res.test <- model.lasso %>% predict(newx = test.data,type='response')
prob.ridge.res.test <- model.ridge %>% predict(newx = test.data,type='response')
prob.elastic.res.test <- model.elastic %>% predict(newx = test.data,type='response')

hist(prob.lasso.res.test)
hist(prob.ridge.res.test)
hist(prob.elastic.res.test)


roc_test.lasso<-pROC_summary(test.y,prob.lasso.res.test)
PR_test.lasso<-Prec_Recall_summary(test.y,prob.lasso.res.test)

roc_test.ridge<-pROC_summary(test.y,prob.ridge.res.test)
PR_test.ridge<-Prec_Recall_summary(test.y,prob.ridge.res.test)

roc_test.elastic<-pROC_summary(test.y,prob.elastic.res.test)
PR_test.elastic<-Prec_Recall_summary(test.y,prob.elastic.res.test)





#get plot and manipulate it
roc_trn.ridge.0<-pROC_summary.0(trn.y,prob.ridge.res)
PR_trn.ridge.0<-Prec_Recall_summary.0(trn.y,prob.ridge.res)

roc_test.ridge.0<-pROC_summary.0(test.y,prob.ridge.res.test)
PR_test.ridge.0<-Prec_Recall_summary.0(test.y,prob.ridge.res.test)



#2022.12.12
#generate figures for the publication
ridge.trn.roc.obj<-roc(trn.y, as.vector(prob.ridge.res), plot=TRUE, 
                       legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)
ridge.trn.roc.auc<-round(auc(trn.y,as.numeric(prob.ridge.res)),3)

ridge.test.roc.obj<-roc(test.y, as.vector(prob.ridge.res.test), plot=TRUE, 
                      legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)
ridge.test.roc.auc<-round(auc(test.y,as.numeric(prob.ridge.res.test)),3)


ridge.trn.auroc<-ggroc(ridge.trn.roc.obj, colour = 'steelblue', size = 2,legacy.axes = T) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', ridge.trn.roc.auc, ')')) +
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+geom_abline(intercept =0 , slope = 1,linetype=2)+coord_equal()

ridge.trn.auprc<-ggplot(data.frame(PR_trn.ridge.0$curve),aes(x=X1,y=X2)) + 
  geom_path(colour = 'steelblue',size=2) +
  labs(x="Recall",y="Precision",
       title=format(PR_trn.ridge.0$auc.integral,digits=3),
       colour="Threshold") + 
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+coord_equal()

ridge.test.auroc<-ggroc(ridge.test.roc.obj, colour = 'steelblue', size = 2,legacy.axes = T) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', ridge.test.roc.auc, ')')) +
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+geom_abline(intercept =0 , slope = 1,linetype=2)+coord_equal()

ridge.test.auprc<-ggplot(data.frame(PR_test.ridge.0$curve),aes(x=X1,y=X2)) + 
  geom_path(colour = 'steelblue',size=2) +
  labs(x="Recall",y="Precision",
          title=format(PR_test.ridge.0$auc.integral,digits=3),
          colour="Threshold") + 
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+coord_equal()

setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/")

jpeg('BEAT_ridge_trn_auroc.jpeg',width=4,height=4,units='in',res=300)
ridge.trn.auroc
dev.off()
jpeg('BEAT_ridge_trn_auprc.jpeg',width=4,height=4,units='in',res=300)
ridge.trn.auprc
dev.off()
jpeg('BEAT_ridge_test_auroc.jpeg',width=4,height=4,units='in',res=300)
ridge.test.auroc
dev.off()
jpeg('BEAT_ridge_test_auprc.jpeg',width=4,height=4,units='in',res=300)
ridge.test.auprc
dev.off()

#old version
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP")
load('GLM_CV_GEP001.Rdata')

#save.image('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Final_GLM_CV_GEP001.Rdata')#case where ridge has the highest performance(lambda.min)
load('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Final_GLM_CV_GEP001.Rdata')
#apply to unknown and aggregate the clinical data

BEAT_GEP_TP53nMT.all
#group_TP53nMT<-pcaData.beat.correct[-TP53MT_WT_ix,462:464]
# 
# group_TP53nMT<-(BEAT_GEP.Data.corrected[-TP53MT_WT_ix,21619:21621])
# New.data<-BEAT_GEP_TP53nMT.all
# 
# prob.lasso.res.newdat <- model.lasso %>% predict(newx = New.data,s=cv.lasso$lambda.min,type='response')
# prob.ridge.res.newdat <- model.ridge %>% predict(newx = New.data,s=cv.ridge$lambda.min,type='response')
# prob.elastic.res.newdat <- model.elastic %>% predict(newx = New.data,s=cv.elastic$lambda.min,type='response')
# 
# hist(prob.lasso.res.newdat)
# hist(prob.ridge.res.newdat)
# hist(prob.elastic.res.newdat)
# 
# # predicted.classes <- ifelse(probabilities.res.newdat > 0.213, 'MT', 'WT')
# # TP53nMT_pred_cyto.df<-cbind(predicted.classes,group_TP53nMT)
# TP53nMT_pred_cyto.df<-as.data.frame(cbind(prob.lasso.res.newdat,
#                                           prob.ridge.res.newdat,
#                                           prob.elastic.res.newdat,
#                                           group_TP53nMT))
# colnames(TP53nMT_pred_cyto.df)[1:3]<-c('lasso.s1','ridge.s1','elastic.s1')

#########################################################################
#see if the models work in TCGA data first
tcga.y<-ifelse(TCGA_GEP.Data.corrected$TP53_mut_stat=='TP53_MUT',1,0)
test.tcga<-data.matrix(TCGA_GEP.Data.corrected[,-(21619:21621)])

prob.lasso.res.tcga <- model.lasso %>% predict(newx = test.tcga,s=cv.lasso$lambda.1se,type='response')
prob.ridge.res.tcga <- model.ridge %>% predict(newx = test.tcga,s=cv.ridge$lambda.1se,type='response')
prob.elastic.res.tcga <- model.elastic %>% predict(newx = test.tcga,s=cv.elastic$lambda.1se,type='response')

# hist(prob.lasso.res.tcga)
# hist(prob.ridge.res.tcga)
# hist(prob.elastic.res.tcga)

# roc_tcga.lasso<-pROC_summary(tcga.y,prob.lasso.res.tcga)
# PR_tcga.lasso<-Prec_Recall_summary(tcga.y,prob.lasso.res.tcga)

roc_tcga.ridge<-pROC_summary(tcga.y,prob.ridge.res.tcga)
PR_tcga.ridge<-Prec_Recall_summary(tcga.y,prob.ridge.res.tcga)

# roc_tcga.elastic<-pROC_summary(tcga.y,prob.elastic.res.tcga)
# PR_tcga.elastic<-Prec_Recall_summary(tcga.y,prob.elastic.res.tcga)

#get TCGA plot:

#get plot and manipulate it
roc_tcga.ridge.0<-pROC_summary.0(tcga.y,prob.ridge.res.tcga)
PR_tcga.ridge.0<-Prec_Recall_summary.0(tcga.y,prob.ridge.res.tcga)

#generate figures for the publication
ridge.tcga.roc.obj<-roc(tcga.y, as.vector(prob.ridge.res.tcga), plot=TRUE, 
                       legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)
ridge.tcga.roc.auc<-round(auc(tcga.y,as.numeric(prob.ridge.res.tcga)),3)


ridge.tcga.auroc<-ggroc(ridge.tcga.roc.obj, colour = 'steelblue', size = 2,legacy.axes = T) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', ridge.tcga.roc.auc, ')')) +
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+geom_abline(intercept =0 , slope = 1,linetype=2)+coord_equal()

ridge.tcga.auprc<-ggplot(data.frame(PR_tcga.ridge.0$curve),aes(x=X1,y=X2)) + 
  geom_path(colour = 'steelblue',size=2) +
  labs(x="Recall",y="Precision",
       title=format(PR_tcga.ridge.0$auc.integral,digits=3),
       colour="Threshold") + 
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+coord_equal()

setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/")

jpeg('TCGA_ridge_auroc.jpeg',width=4,height=4,units='in',res=300)
ridge.tcga.auroc
dev.off()
jpeg('TCGA_ridge_auprc.jpeg',width=4,height=4,units='in',res=300)
ridge.tcga.auprc
dev.off()


gc()

###################################################################
#skip KM screening for now.

#################################################################
#new GLM for TCGA Mut-like: old version 

setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP")
#save.image('GLM_CV_GEP005_Mut_stat.Rdata')#one that before using cv$fit.preval for train data score
load('GLM_CV_GEP005_Mut_stat.Rdata')
gc()


################################################
#get TCGA clinical data - not many data to look at
TCGA_clinical.final
LSC17_score.tcga$TCGA.Patient.ID.y<-sub('\\.03.*','',LSC17_score.tcga$LabId)
LSC17_score.comb.tcga$TCGA.Patient.ID.y<-sub('\\.03.*','',LSC17_score.comb.tcga$LabId)

TCGA_clinical.final.n<-right_join(TCGA_clinical.final,LSC17_score.tcga,by='TCGA.Patient.ID.y')
TCGA_clinical.final.n<-right_join(TCGA_clinical.final.n,LSC17_score.comb.tcga,by='TCGA.Patient.ID.y')

#fix Mut-like cases using KM results -lasso and elastic only:
#1.unknown-MUT-like+TP53 vs.unknown-WT-like 
KM04.01.df#lasso
KM04.02.df#ridge
# KM01.lasso.id<-KM01.lasso$LabId[1:18]
# KM01.elastic.id<-KM01.elastic$LabId[1:20]

#2.unknown-MUT-like vs. unknown-WT-like
KM05.01.df#lasso only 1
KM05.02.df#elastic 8 unknown(other-5 and 6)
KM05.lasso.id<-KM05.01.df$LabId[1]
KM05.elastic.id<-KM05.02.df$LabId[1:8]
#GLM_id<-c()
GLM_id.tcga<-KM05.elastic.id############
cls.rep.tcga<-(TCGA_clinical.final.n$TP53_mut_stat)
#View(cbind(prob.elastic.res.tcga,tcga.y))

###################################################################################
#have not updated to newer version yet
#2022.05.19
#New classification analysis: # this is hierarchcical model
#1.additional GLM with ridge identified TP53MUT+TP53MUT-like WT(based on the KM screening) vs. WT
#to identify TP53MUT-like WT in TCGA 
#2. additional GLM with ridge identified TP53MUT-like WT(based on the KM screening) vs. WT
#to identify TP53MUT-like WT in TCGA (Wt only)
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP")
#save.image('GLM_CV_GEP003.Rdata')
#load('GLM_CV_GEP003.Rdata')
load('GLM_CV_GEP004.Rdata')
gc()

source("C:/Users/yklee/Desktop/Sachs Lab/R_project/TP53/TP53GEP_functions.r")
packages<-c('edgeR','limma','dplyr','tidyr','sva','tidyverse',
            'ggplot2','glmnet','caret','pROC','survival','survminer')

for(p in packages){
  library(p,character.only = T)
}


TP53MT_WT_ix<-which(BEAT_GEP.Data.corrected$TP53_mut_stat != 'unknown')
TP53_MT_nMT<-rep(0,dim(BEAT_GEP.Data.corrected)[1])
TP53_MT_nMT[which(BEAT_GEP.Data.corrected$TP53_mut_stat=='TP53_MUT')]<-1
TP53_MT_WT<-TP53_MT_nMT[TP53MT_WT_ix]
TP53_WT_ix<-which(BEAT_GEP.Data.corrected$TP53_mut_stat == 'TP53_WT')
#MT vs. WT
KM01.ridge.42id<-KM01.ridge$LabId[1:42]

#2022.06.07
#saveRDS(KM01.ridge,'class ids/BEAT_WT_KM01_ridge.rds')
#saveRDS(KM02.ridge,'class ids/BEAT_WTUnk_KM02_ridge.rds')

##########################################################################################
#2. additional GLM with ridge identified TP53MUT-like WT(based on the KM screening) vs. WT
#to identify TP53MUT-like WT in TCGA (Wt only)
#use above MUT-like WT + MUT as MUT classes to retrain BEAT AML data
class_id<-KM01.ridge.42id


##1.4. run LASSO train/test  
BEAT_GEP.Data.corrected$TP53_mutlike_stat<-rep(0,dim(BEAT_GEP.Data.corrected)[1])
BEAT_GEP.Data.corrected$TP53_mutlike_stat[BEAT_GEP.Data.corrected$sampleID %in% class_id]<-1
TP53_group_sub.y<-as.factor(BEAT_GEP.Data.corrected$TP53_mutlike_stat[TP53_WT_ix])

#subset data for MT_WT and unknown and with all vs. DEGs
#1.all genes
#get only WT
#BEAT_GEP_TP53MT_WT.all<-data.matrix(BEAT_GEP.Data.corrected[TP53MT_WT_ix,-(21619:21622)])
BEAT_GEP_TP53WT.all<-data.matrix(BEAT_GEP.Data.corrected[TP53_WT_ix,-(21619:21622)])
BEAT_GEP_TP53nMT.all<-data.matrix(BEAT_GEP.Data.corrected[-TP53MT_WT_ix,-(21619:21622)])

setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP")



#set.seed(9011230)
#set.seed(90112312)
#set.seed(9011233)
set.seed(11233)
#set.seed(901123)
trn.idx<-TP53_group_sub.y %>% createDataPartition(p=0.6, list = FALSE)
trn.data<-data.matrix(BEAT_GEP_TP53WT.all[trn.idx,])
trn.y<-(TP53_group_sub.y[trn.idx])
test.data<-data.matrix(BEAT_GEP_TP53WT.all[-trn.idx,])
test.y<-(TP53_group_sub.y[-trn.idx])

#set 10-fold cv for the model comparison
set.seed(2843)
#set.seed(284)
#set.seed(283)
flds <- createFolds(trn.y, k = 10, list = TRUE, returnTrain = FALSE)
foldids = rep(1,length(trn.y))
foldids[flds$Fold02] = 2
foldids[flds$Fold03] = 3
foldids[flds$Fold04] = 4
foldids[flds$Fold05] = 5
foldids[flds$Fold06] = 6
foldids[flds$Fold07] = 7
foldids[flds$Fold08] = 8
foldids[flds$Fold09] = 9
foldids[flds$Fold10] = 10
#try ridge regression and elastic net 
cv.lasso<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=1,keep = TRUE,
                    type.measure = 'class',foldid=foldids)
cv.ridge<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=0,keep = TRUE,
                    type.measure = 'class',foldid=foldids)
cv.elastic<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=0.5,keep = TRUE,
                      type.measure = 'class',foldid=foldids)
# cv.lasso<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=1,keep = TRUE,
#                     type.measure = 'class')
# cv.ridge<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=0,keep = TRUE,
#                     type.measure = 'class')
# cv.elastic<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=0.5,keep = TRUE,
#                       type.measure = 'class')

# plot(cv.lasso)
# plot(cv.ridge)
# plot(cv.elastic)
model.lasso <- glmnet(trn.data, trn.y, alpha = 1, family = "binomial",
                      lambda = cv.lasso$lambda.1se)
model.ridge <- glmnet(trn.data, trn.y, alpha = 0, family = "binomial",
                      lambda = cv.ridge$lambda.1se)
model.elastic <- glmnet(trn.data, trn.y, alpha = 0.5, family = "binomial",
                        lambda = cv.elastic$lambda.1se)

# 
# model.lasso <- glmnet(trn.data, trn.y, alpha = 1, family = "binomial",
#                       lambda = cv.lasso$lambda.min)
# model.ridge <- glmnet(trn.data, trn.y, alpha = 0, family = "binomial",
#                       lambda = cv.ridge$lambda.min)
# model.elastic <- glmnet(trn.data, trn.y, alpha = 0.5, family = "binomial",
#                         lambda = cv.elastic$lambda.min)

prob.lasso.res <- model.lasso %>% predict(newx = trn.data,s=cv.lasso$lambda.1se,type='response')
prob.ridge.res <- model.ridge %>% predict(newx = trn.data,s=cv.ridge$lambda.1se,type='response')
prob.elastic.res <- model.elastic %>% predict(newx = trn.data,s=cv.elastic$lambda.1se,type='response')


prob.lasso.res.test <- model.lasso %>% predict(newx = test.data,type='response')
prob.ridge.res.test <- model.ridge %>% predict(newx = test.data,type='response')
prob.elastic.res.test <- model.elastic %>% predict(newx = test.data,type='response')


prob.lasso.res.tcga <- model.lasso %>% predict(newx = test.tcga,s=cv.lasso$lambda.1se,type='response')
prob.ridge.res.tcga <- model.ridge %>% predict(newx = test.tcga,s=cv.ridge$lambda.1se,type='response')
prob.elastic.res.tcga <- model.elastic %>% predict(newx = test.tcga,s=cv.elastic$lambda.1se,type='response')

out.performance<-get_AUC_summary(trn.y,test.y,tcga.y,prob.lasso.res,prob.ridge.res,prob.elastic.res,
                                 prob.lasso.res.test,prob.ridge.res.test,prob.elastic.res.test,
                                 prob.lasso.res.tcga,prob.ridge.res.tcga,prob.elastic.res.tcga)
#write.csv(out.performance,'out.performance.csv')


roc_trn.ridge<-pROC_summary(trn.y,prob.ridge.res)
PR_trn.ridge<-Prec_Recall_summary(trn.y,prob.ridge.res)

roc_test.ridge<-pROC_summary(test.y,prob.ridge.res.test)
PR_test.ridge<-Prec_Recall_summary(test.y,prob.ridge.res.test)

roc_test.ridge<-pROC_summary(tcga.y,prob.ridge.res.tcga)
PR_test.ridge<-Prec_Recall_summary(tcga.y,prob.ridge.res.tcga)


#generate perfomance plot

#get plot and manipulate it
roc_trn.ridge.0<-pROC_summary.0(trn.y,prob.ridge.res)
PR_trn.ridge.0<-Prec_Recall_summary.0(trn.y,prob.ridge.res)

roc_test.ridge.0<-pROC_summary.0(test.y,prob.ridge.res.test)
PR_test.ridge.0<-Prec_Recall_summary.0(test.y,prob.ridge.res.test)

roc_tcga.ridge.0<-pROC_summary.0(tcga.y,prob.ridge.res.tcga)
PR_tcga.ridge.0<-Prec_Recall_summary.0(tcga.y,prob.ridge.res.tcga)



#2022.12.12
#generate figures for the publication
ridge.trn.roc.obj<-roc(trn.y, as.vector(prob.ridge.res), plot=TRUE, 
                       legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)
ridge.trn.roc.auc<-round(auc(trn.y,as.numeric(prob.ridge.res)),3)

ridge.test.roc.obj<-roc(test.y, as.vector(prob.ridge.res.test), plot=TRUE, 
                        legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)
ridge.test.roc.auc<-round(auc(test.y,as.numeric(prob.ridge.res.test)),3)

ridge.tcga.roc.obj<-roc(tcga.y, as.vector(prob.ridge.res.tcga), plot=TRUE, 
                        legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)
ridge.tcga.roc.auc<-round(auc(tcga.y,as.numeric(prob.ridge.res.tcga)),3)


ridge.trn.auroc<-ggroc(ridge.trn.roc.obj, colour = 'steelblue', size = 2,legacy.axes = T) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', ridge.trn.roc.auc, ')')) +
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+geom_abline(intercept =0 , slope = 1,linetype=2)+coord_equal()

ridge.trn.auprc<-ggplot(data.frame(PR_trn.ridge.0$curve),aes(x=X1,y=X2)) + 
  geom_path(colour = 'steelblue',size=2) +
  labs(x="Recall",y="Precision",
       title=format(PR_trn.ridge.0$auc.integral,digits=3),
       colour="Threshold") + 
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+coord_equal()

ridge.test.auroc<-ggroc(ridge.test.roc.obj, colour = 'steelblue', size = 2,legacy.axes = T) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', ridge.test.roc.auc, ')')) +
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+geom_abline(intercept =0 , slope = 1,linetype=2)+coord_equal()

ridge.test.auprc<-ggplot(data.frame(PR_test.ridge.0$curve),aes(x=X1,y=X2)) + 
  geom_path(colour = 'steelblue',size=2) +
  labs(x="Recall",y="Precision",
       title=format(PR_test.ridge.0$auc.integral,digits=3),
       colour="Threshold") + 
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+coord_equal()

#tcga
ridge.tcga.auroc<-ggroc(ridge.tcga.roc.obj, colour = 'steelblue', size = 2,legacy.axes = T) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', ridge.tcga.roc.auc, ')')) +
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+geom_abline(intercept =0 , slope = 1,linetype=2)+coord_equal()

ridge.tcga.auprc<-ggplot(data.frame(PR_tcga.ridge.0$curve),aes(x=X1,y=X2)) + 
  geom_path(colour = 'steelblue',size=2) +
  labs(x="Recall",y="Precision",
       title=format(PR_tcga.ridge.0$auc.integral,digits=3),
       colour="Threshold") + 
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+coord_equal()

setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/")

jpeg('newGLM_BEAT_ridge_trn_auroc_mtl.jpeg',width=4,height=4,units='in',res=300)
ridge.trn.auroc
dev.off()
jpeg('newGLM_BEAT_ridge_trn_auprc_mtl.jpeg',width=4,height=4,units='in',res=300)
ridge.trn.auprc
dev.off()
jpeg('newGLM_BEAT_ridge_test_auroc_mtl.jpeg',width=4,height=4,units='in',res=300)
ridge.test.auroc
dev.off()
jpeg('newGLM_BEAT_ridge_test_auprc_mtl.jpeg',width=4,height=4,units='in',res=300)
ridge.test.auprc
dev.off()
jpeg('newGLM_TCGA_ridge_auroc_mut_cls_instead.jpeg',width=4,height=4,units='in',res=300)
ridge.tcga.auroc
dev.off()
jpeg('newGLM_TCGA_ridge_auprc_mut_cls_instead.jpeg',width=4,height=4,units='in',res=300)
ridge.tcga.auprc
dev.off()





















#############################################################################################
#check for the Overall survival as a TP53Mut-like cutoff in the unknown settings 
#BEAT AML
# BEAT_ct_clinc.RNA_avail.all.01<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_AML/BEAT_AML_preprocessed/BEAT_ct_clinc_RNA_avail_all_01.rds')
BEAT_ct_clinc.RNA_avail.all.01.WES<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_AML/BEAT_AML_preprocessed/BEAT_ct_RNA_01.WES.allmut.IP.rds')
Mut_list<-c('DNMT3A','TET2','IDH1','IDH2','WT1','EZH2','ASXL1','SRSF2','STAG2','RAD21',
            'FLT3','KIT','KRAS','NRAS','PTPN11','NF1','JAK2','BCOR','GATA2','ETV6',
            'CEBPA','RUNX1','CDKN2A','NPM1','MYC','U2AF1','PHF6','SMC3','SMC1A','CSF3R',
            'PDS5B','SF3B1','BRINP3','CROCC','PLCE1','SPEN')

Mut_list.n<-str_c('WES',Mut_list,sep='_')
#use clinical summary data for NPM1 and FLT3-ITD - these two are complete
Mut_list.n[Mut_list.n =='WES_NPM1']<-'NPM1'
Mut_list.n<-c(Mut_list.n,'FLT3-ITD')

mut_ix<-which(colnames(BEAT_ct_clinc.RNA_avail.all.01.WES)%in% Mut_list.n)
#clinical data that we are interested in 
clinical_summary_var<-BEAT_ct_clinc.RNA_avail.all.01.WES[ ,c(1,2,8,9,31,47,48,49,54,55,62,63,67,74,84,mut_ix)]

library(readxl)
ZS_anno_brief<-read_excel('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/BEAT AML dataset Karyotype annotation_brief.xlsx')
#ZS_anno<-read_excel('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/BEAT AML dataset Karyotype annotation.xlsx')
rna_ix<-which(ZS_anno_brief$LabId%in% BEAT_ct_clinc.RNA_avail.all.01.WES$LabId)
ZS_anno_brief_rna_count<-ZS_anno_brief[rna_ix,]
#ZS_anno_cytogenetics

ZS_anno_group<-ZS_anno_brief_rna_count$`ZS Cytogenetics Annotation`
ZS_anno_group[which(ZS_anno_group=='not avial')]<-'not avail'
ZS_anno_group[which(ZS_anno_group=='noncomplex-nonmono-nonnormal')]<-'non-comp/mono/norm'
unique(ZS_anno_group)
ZS_anno_brief_rna_count$ZS_anno_group<-ZS_anno_group
ZS_anno_brief_rna_count.use<-ZS_anno_brief_rna_count[,c(1,5)]

clinical_summary_var<-right_join(clinical_summary_var,ZS_anno_brief_rna_count.use,by='LabId')

#apply model to all the BEAT AML data - regardless of it is training data or not
#2022.05.31:however, use scores from the cv$predmat each fold's prediction score
#generated from other 9 fold trained from model
#a. compare with the previous version(apply model to all the data)
# trn.data<-data.matrix(BEAT_GEP_TP53MT_WT.all[trn.idx,])
# test.data<-data.matrix(BEAT_GEP_TP53MT_WT.all[-trn.idx,])

BEAT.all<-data.matrix(BEAT_GEP.Data.corrected[,-(21619:21621)])
prob.ridge.res.all.old <- model.ridge %>% predict(newx = BEAT.all,s=cv.ridge$lambda.min,type='response')

#merge using fit.preval in cv + model:test +model:unknown
ref_id<-data.frame("ID"=BEAT_GEP.Data.corrected$sampleID)
library(boot)
prob.lasso.res.all.trn<-rownames_to_column(data.frame('s0'=inv.logit(cv.lasso$fit.preval[,14])),var='ID')#min.lambda index 14
prob.ridge.res.all.trn<-rownames_to_column(data.frame('s0'=inv.logit(cv.ridge$fit.preval[,81])),var='ID')#min.lambda index 81
prob.elastic.res.all.trn<-rownames_to_column(data.frame('s0'=inv.logit(cv.elastic$fit.preval[,44])),var='ID')#min.lambda index 44

prob.lasso.res.all.test<-rownames_to_column(as.data.frame(predict(model.lasso,newx=test.data,type='response')),var='ID')
prob.ridge.res.all.test<-rownames_to_column(as.data.frame(predict(model.ridge,newx=test.data,type='response')),var='ID')
prob.elastic.res.all.test<-rownames_to_column(as.data.frame(predict(model.elastic,newx=test.data,type='response')),var='ID')

prob.lasso.res.all.unk<-rownames_to_column(as.data.frame(predict(model.lasso,newx=New.data,type='response')),var='ID')
prob.ridge.res.all.unk<-rownames_to_column(as.data.frame(predict(model.ridge,newx=New.data,type='response')),var='ID')
prob.elastic.res.all.unk<-rownames_to_column(as.data.frame(predict(model.elastic,newx=New.data,type='response')),var='ID')

BEAT.lasso.rbind<-rbind(prob.lasso.res.all.trn,prob.lasso.res.all.test,prob.lasso.res.all.unk)
BEAT.ridge.rbind<-rbind(prob.ridge.res.all.trn,prob.ridge.res.all.test,prob.ridge.res.all.unk)
BEAT.elastic.rbind<-rbind(prob.elastic.res.all.trn,prob.elastic.res.all.test,prob.elastic.res.all.unk)

prob.lasso.res.all<-right_join(ref_id,BEAT.lasso.rbind,by='ID')%>% column_to_rownames(var='ID')
prob.ridge.res.all<-right_join(ref_id,BEAT.ridge.rbind,by='ID')%>% column_to_rownames(var='ID')
prob.elastic.res.all<-right_join(ref_id,BEAT.elastic.rbind,by='ID')%>% column_to_rownames(var='ID')

plot(prob.ridge.res.all$s0,prob.ridge.res.all.old)
cor(prob.ridge.res.all$s0,prob.ridge.res.all.old,method='pearson')#0.6311
cor(prob.ridge.res.all$s0,prob.ridge.res.all.old,method='spearman')#rank 0.9438


#combine lasso results (MT_WT_unknown(lasso identified class))+ clinical data
#TP53nMT_pred_cyto.df #predicted class in unknown + group+sample id
#BEAT_GEP.Data.corrected
#TP53_MT_class_cyto<-BEAT_GEP.Data.corrected[TP53MT_WT_ix,21619:21621]
# TP53_MT_class_cyto$lasso.s1<-ifelse(TP53_MT_class_cyto$TP53_mut_stat=='TP53_MUT',1,0)
# TP53_MT_class_cyto$ridge.s1<-ifelse(TP53_MT_class_cyto$TP53_mut_stat=='TP53_MUT',1,0)
# TP53_MT_class_cyto$elastic.s1<-ifelse(TP53_MT_class_cyto$TP53_mut_stat=='TP53_MUT',1,0)
#TP53_MT_class_cyto<-TP53_MT_class_cyto[,c(4,5,6,1,2,3)]
#should I merge?
#BEAT_LASSO_final_res<-bind_rows(TP53_MT_class_cyto,TP53nMT_pred_cyto.df)

TP53_MT_class_cyto<-BEAT_GEP.Data.corrected[,21619:21621]
TP53_MT_class_cyto$lasso.s1<-as.numeric(prob.lasso.res.all$s0)
TP53_MT_class_cyto$ridge.s1<-as.numeric(prob.ridge.res.all$s0)
TP53_MT_class_cyto$elastic.s1<-as.numeric(prob.elastic.res.all$s0)
TP53_MT_class_cyto<-TP53_MT_class_cyto[,c(4,5,6,1,2,3)]
BEAT_LASSO_final_res<-TP53_MT_class_cyto


#merge with clinical data based on the sample_id
colnames(BEAT_LASSO_final_res)[4]<-'LabId'
BEAT_LASSO_final_res.clinical<-right_join(BEAT_LASSO_final_res,
                                          clinical_summary_var,by = 'LabId')
View(BEAT_LASSO_final_res.clinical)

dat.full<-BEAT_LASSO_final_res.clinical
dat.full$OS_status<-ifelse(dat.full$vitalStatus=='Dead',1,0)
dat.full$TP53.z<-ifelse(dat.full$TP53_mut_stat=='TP53_MUT',1,0)
#TP53MTvsWT
dat.sub01.mtwt<-dat.full[which(dat.full$TP53_mut_stat=='TP53_MUT' | dat.full$TP53_mut_stat=='TP53_WT'),]
dat.sub01.wt<-dat.full[which(dat.full$TP53_mut_stat=='TP53_WT'),]

###below will not be used for KM analysis
dat.sub01.unk<-dat.full[-which(dat.full$TP53_mut_stat=='TP53_MUT' | dat.full$TP53_mut_stat=='TP53_WT'),]
dat.sub01.wtunk<-dat.full[-which(dat.full$TP53_mut_stat=='TP53_MUT'),]
###

#TP53mut vs. all(WT+unknown)
dat.full.obj<-Surv(time = dat.full$overallSurvival,
                   event = dat.full$OS_status)
#TP53mut vs. all(WT)
dat.sub01.mtwt.obj<-Surv(time = dat.sub01.mtwt$overallSurvival,
                         event = dat.sub01.mtwt$OS_status)

fit.all <- survfit(dat.full.obj ~ TP53.z, data = dat.full)
fit.mtwt <- survfit(dat.sub01.mtwt.obj ~ TP53.z, data = dat.sub01.mtwt)

g.all<-ggsurvplot(fit.all, data = dat.full, pval = TRUE)
g.mtwt<-ggsurvplot(fit.mtwt, data = dat.sub01.mtwt, pval = TRUE)
pdf('BEAT_AML_TP53mut_wtunk.pdf')
print(g.all)
dev.off()
pdf('BEAT_AML_TP53mutwt.pdf')
print(g.mtwt)
dev.off()

# surv_pvalue(fit.mtwt)[2]
# surv_pvalue(fit.mtnmt)[2]
fit.coxph.all <- coxph(dat.full.obj ~ TP53.z,data=dat.full)
fit.coxph.mtwt <- coxph(dat.sub01.mtwt.obj ~ TP53.z,data=dat.sub01.mtwt)
#pvalue in cox
#summary(fit.coxph.mtwt)$coefficients[5]
#summary(fit.coxph.mtnmt)$coefficients[5]

#find unknown cutoffs 2 scenarios
#1)MUT_like vs. other unknown and 2) MUT_like vs. WT
#data without TP53MUT and WT
dat.sub01.unk
#data without TP53MUT
dat.sub01.wtunk
#data with only MT and WT
dat.sub01.mtwt

dat.unk.obj<-Surv(time = dat.sub01.unk$overallSurvival,
                  event = dat.sub01.unk$OS_status)
dat.wtunk.obj<-Surv(time = dat.sub01.wtunk$overallSurvival,
                    event = dat.sub01.wtunk$OS_status)
dat.wt.obj<-Surv(time = dat.sub01.wt$overallSurvival,
                 event = dat.sub01.wt$OS_status)
# dat.sub01.mtwt<-dat.full[which(dat.full$TP53_mut_stat=='TP53_MUT' | dat.full$TP53_mut_stat=='TP53_WT'),]
dat.mtwt.obj<-Surv(time = dat.sub01.mtwt$overallSurvival,
                   event = dat.sub01.mtwt$OS_status)


ix01.lasso<-order(dat.sub01.wt$lasso.s1,decreasing = T)
ix01.ridge<-order(dat.sub01.wt$ridge.s1,decreasing = T)
ix01.elastic<-order(dat.sub01.wt$elastic.s1,decreasing = T)

KM01.lasso<-data.frame('prob_cutoff'=dat.sub01.wt[ix01.lasso,]$lasso.s1,
                       'KM_pval'=rep(0,dim(dat.sub01.wt)[1]),
                       'LabId'=dat.sub01.wt[ix01.lasso,]$LabId)

KM01.ridge<-data.frame('prob_cutoff'=dat.sub01.wt[ix01.ridge,]$ridge.s1,
                       'KM_pval'=rep(0,dim(dat.sub01.wt)[1]),
                       'LabId'=dat.sub01.wt[ix01.ridge,]$LabId)

KM01.elastic<-data.frame('prob_cutoff'=dat.sub01.wt[ix01.elastic,]$elastic.s1,
                         'KM_pval'=rep(0,dim(dat.sub01.wt)[1]),
                         'LabId'=dat.sub01.wt[ix01.elastic,]$LabId)

for(i in 1:length(ix01.lasso)){
  #for(i in 1:5){
  
  ixx.lasso<-ix01.lasso[1:i]
  ixx.ridge<-ix01.ridge[1:i]
  ixx.elastic<-ix01.elastic[1:i]
  
  # print(ixx.lasso)
  # print(ixx.ridge)
  # print(ixx.elastic)
  
  mut_like_status<-rep(0,dim(dat.sub01.wt)[1])
  mut_like_status.lasso<-mut_like_status
  mut_like_status.ridge<-mut_like_status
  mut_like_status.elastic<-mut_like_status
  
  mut_like_status.lasso[ixx.lasso]<-1
  mut_like_status.ridge[ixx.ridge]<-1
  mut_like_status.elastic[ixx.elastic]<-1
  dat.sub01.wt$TP53.lasso<-mut_like_status.lasso
  dat.sub01.wt$TP53.ridge<-mut_like_status.ridge
  dat.sub01.wt$TP53.elastic<-mut_like_status.elastic
  
  fit.km.lasso<-survfit(dat.wt.obj ~TP53.lasso, data=dat.sub01.wt)
  #fit.cph.lasso<-coxph(dat.sub02.obj ~TP53.lasso, data=dat.sub02)
  fit.km.ridge<-survfit(dat.wt.obj ~TP53.ridge, data=dat.sub01.wt)
  #fit.cph.ridge<-coxph(dat.sub02.obj ~TP53.ridge, data=dat.sub02)
  fit.km.elastic<-survfit(dat.wt.obj ~TP53.elastic, data=dat.sub01.wt)
  #fit.cph.elastic<-coxph(dat.sub02.obj ~TP53.elastic, data=dat.sub02)
  
  KM01.lasso[i,2]<-surv_pvalue(fit.km.lasso)[2]
  #CPH01.lasso[i,2]<-summary(fit.cph.lasso)$coefficients[5]
  KM01.ridge[i,2]<-surv_pvalue(fit.km.ridge)[2]
  #CPH01.ridge[i,2]<-summary(fit.cph.ridge)$coefficients[5]
  KM01.elastic[i,2]<-surv_pvalue(fit.km.elastic)[2]
  #CPH01.elastic[i,2]<-summary(fit.cph.elastic)$coefficients[5]
  
}
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/Survival analysis/")
pdf('BEAT_regressions_cutoff_unknown_KM.pdf')
plot(1:377,KM01.lasso$KM_pval)
abline(h=0.05,col='red')
axis(1,1:377,round(KM01.lasso$prob_cutoff,digits=2),line=2)
plot(1:377,KM01.ridge$KM_pval)
abline(h=0.05,col='red')
axis(1,1:377,round(KM01.ridge$prob_cutoff,digits=2),line=2)
plot(1:377,KM01.elastic$KM_pval)
abline(h=0.05,col='red')
axis(1,1:377,round(KM01.elastic$prob_cutoff,digits=2),line=2)
dev.off()


###
ix02.lasso<-order(dat.sub01.wtunk$lasso.s1,decreasing = T)
ix02.ridge<-order(dat.sub01.wtunk$ridge.s1,decreasing = T)
ix02.elastic<-order(dat.sub01.wtunk$elastic.s1,decreasing = T)


KM02.lasso<-data.frame('prob_cutoff'=dat.sub01.wtunk[ix02.lasso,]$lasso.s1,
                       'KM_pval'=rep(0,dim(dat.sub01.wtunk)[1]),
                       'LabId'=dat.sub01.wtunk[ix02.lasso,]$LabId)

KM02.ridge<-data.frame('prob_cutoff'=dat.sub01.wtunk[ix02.ridge,]$ridge.s1,
                       'KM_pval'=rep(0,dim(dat.sub01.wtunk)[1]),
                       'LabId'=dat.sub01.wtunk[ix02.ridge,]$LabId)
KM02.elastic<-data.frame('prob_cutoff'=dat.sub01.wtunk[ix02.elastic,]$elastic.s1,
                         'KM_pval'=rep(0,dim(dat.sub01.wtunk)[1]),
                         'LabId'=dat.sub01.wtunk[ix02.elastic,]$LabId)


for(i in 1:length(ix02.lasso)){
  #for(i in 1:5){
  
  ixx.lasso<-ix02.lasso[1:i]
  ixx.ridge<-ix02.ridge[1:i]
  ixx.elastic<-ix02.elastic[1:i]
  
  # print(ixx.lasso)
  # print(ixx.ridge)
  # print(ixx.elastic)
  
  mut_like_status<-rep(0,dim(dat.sub01.wtunk)[1])
  mut_like_status.lasso<-mut_like_status
  mut_like_status.ridge<-mut_like_status
  mut_like_status.elastic<-mut_like_status
  
  mut_like_status.lasso[ixx.lasso]<-1
  mut_like_status.ridge[ixx.ridge]<-1
  mut_like_status.elastic[ixx.elastic]<-1
  dat.sub01.wtunk$TP53.lasso<-mut_like_status.lasso
  dat.sub01.wtunk$TP53.ridge<-mut_like_status.ridge
  dat.sub01.wtunk$TP53.elastic<-mut_like_status.elastic
  
  fit.km.lasso<-survfit(dat.wtunk.obj ~TP53.lasso, data=dat.sub01.wtunk)
  #fit.cph.lasso<-coxph(dat.sub03.obj ~TP53.lasso, data=dat.sub03)
  fit.km.ridge<-survfit(dat.wtunk.obj ~TP53.ridge, data=dat.sub01.wtunk)
  #fit.cph.ridge<-coxph(dat.sub03.obj ~TP53.ridge, data=dat.sub03)
  fit.km.elastic<-survfit(dat.wtunk.obj ~TP53.elastic, data=dat.sub01.wtunk)
  #fit.cph.elastic<-coxph(dat.sub03.obj ~TP53.elastic, data=dat.sub03)
  
  KM02.lasso[i,2]<-surv_pvalue(fit.km.lasso)[2]
  #CPH02.lasso[i,2]<-summary(fit.cph.lasso)$coefficients[5]
  KM02.ridge[i,2]<-surv_pvalue(fit.km.ridge)[2]
  #CPH02.ridge[i,2]<-summary(fit.cph.ridge)$coefficients[5]
  KM02.elastic[i,2]<-surv_pvalue(fit.km.elastic)[2]
  #CPH02.elastic[i,2]<-summary(fit.cph.elastic)$coefficients[5]
  
}
#pdf('BEAT_regressions_cutoff_WTunknown_KM.pdf')
plot(1:425,KM02.lasso$KM_pval)
abline(h=0.05,col='red')
axis(1,1:425,round(KM02.lasso$prob_cutoff,digits=2),line=2)
plot(1:425,KM02.ridge$KM_pval)
abline(h=0.05,col='red')
axis(1,1:425,round(KM02.ridge$prob_cutoff,digits=2),line=2)
plot(1:425,KM02.elastic$KM_pval)
abline(h=0.05,col='red')
axis(1,1:425,round(KM02.elastic$prob_cutoff,digits=2),line=2)
#dev.off()

ix03.lasso<-order(dat.full$lasso.s1,decreasing = T)
ix03.ridge<-order(dat.full$ridge.s1,decreasing = T)
ix03.elastic<-order(dat.full$elastic.s1,decreasing = T)
#add random index to test the inclusive criteria (leave 27 MUT status)
mut_ixx<-which(dat.full$TP53_mut_stat=='TP53_MUT')
#set.seed(1234)
set.seed(12345)
random_ix01<-sample(dim(dat.full)[1],dim(dat.full)[1])
random_ix02<-sample(dim(dat.full)[1],dim(dat.full)[1])
random_ix03<-sample(dim(dat.full)[1],dim(dat.full)[1])

set.seed(123)
ix03.random01<-c(sample(mut_ixx),random_ix01[!random_ix01 %in% mut_ixx])
set.seed(12345)
ix03.random02<-c(sample(mut_ixx),random_ix02[!random_ix02 %in% mut_ixx])
set.seed(1423)
ix03.random03<-c(sample(mut_ixx),random_ix03[!random_ix03 %in% mut_ixx])


KM03.lasso<-data.frame('prob_cutoff'=dat.full[ix03.lasso,]$lasso.s1,
                       'KM_pval'=rep(0,dim(dat.full)[1]),
                       'LabId'=dat.full[ix03.lasso,]$LabId)
KM03.ridge<-data.frame('prob_cutoff'=dat.full[ix03.ridge,]$ridge.s1,
                       'KM_pval'=rep(0,dim(dat.full)[1]),
                       'LabId'=dat.full[ix03.ridge,]$LabId)
KM03.elastic<-data.frame('prob_cutoff'=dat.full[ix03.elastic,]$elastic.s1,
                         'KM_pval'=rep(0,dim(dat.full)[1]),
                         'LabId'=dat.full[ix03.elastic,]$LabId)
#create random 
KM03.random01<-data.frame('prob_cutoff'=1:461,
                          'KM_pval'=rep(0,dim(dat.full)[1]))
KM03.random02<-data.frame('prob_cutoff'=1:461,
                          'KM_pval'=rep(0,dim(dat.full)[1]))
KM03.random03<-data.frame('prob_cutoff'=1:461,
                          'KM_pval'=rep(0,dim(dat.full)[1]))


for(i in 1:length(ix03.lasso)){
  #for(i in 1:5){
  
  ixx.lasso<-ix03.lasso[1:i]
  ixx.ridge<-ix03.ridge[1:i]
  ixx.elastic<-ix03.elastic[1:i]
  ixx.random01<-ix03.random01[1:i]
  ixx.random02<-ix03.random02[1:i]
  ixx.random03<-ix03.random03[1:i]
  
  # print(ixx.lasso)
  # print(ixx.ridge)
  # print(ixx.elastic)
  
  mut_like_status<-rep(0,dim(dat.full)[1])
  mut_like_status.lasso<-mut_like_status
  mut_like_status.ridge<-mut_like_status
  mut_like_status.elastic<-mut_like_status
  mut_like_status.random01<-mut_like_status
  mut_like_status.random02<-mut_like_status
  mut_like_status.random03<-mut_like_status
  
  mut_like_status.lasso[ixx.lasso]<-1
  mut_like_status.ridge[ixx.ridge]<-1
  mut_like_status.elastic[ixx.elastic]<-1
  mut_like_status.random01[ixx.random01]<-1
  mut_like_status.random02[ixx.random02]<-1
  mut_like_status.random03[ixx.random03]<-1
  
  dat.full$TP53.lasso<-mut_like_status.lasso
  dat.full$TP53.ridge<-mut_like_status.ridge
  dat.full$TP53.elastic<-mut_like_status.elastic
  dat.full$TP53.random01<-mut_like_status.random01
  dat.full$TP53.random02<-mut_like_status.random02
  dat.full$TP53.random03<-mut_like_status.random03
  
  
  fit.km.lasso<-survfit(dat.full.obj ~TP53.lasso, data=dat.full)
  
  fit.km.ridge<-survfit(dat.full.obj ~TP53.ridge, data=dat.full)
  
  fit.km.elastic<-survfit(dat.full.obj ~TP53.elastic, data=dat.full)
  
  fit.km.random01<-survfit(dat.full.obj ~TP53.random01, data=dat.full)
  
  fit.km.random02<-survfit(dat.full.obj ~TP53.random02, data=dat.full)
  
  fit.km.random03<-survfit(dat.full.obj ~TP53.random03, data=dat.full)
  
  
  KM03.lasso[i,2]<-surv_pvalue(fit.km.lasso)[2]
  KM03.ridge[i,2]<-surv_pvalue(fit.km.ridge)[2]
  KM03.elastic[i,2]<-surv_pvalue(fit.km.elastic)[2]
  #
  KM03.random01[i,2]<-surv_pvalue(fit.km.random01)[2]
  KM03.random02[i,2]<-surv_pvalue(fit.km.random02)[2]
  KM03.random03[i,2]<-surv_pvalue(fit.km.random03)[2]
  
}

pdf('BEAT_regressions_cutoff_MTWTunknown_KM.pdf')
plot(1:461,KM03.lasso$KM_pval)
abline(h=0.05,col='red')
axis(1,1:461,round(KM03.lasso$prob_cutoff,digits=2),line=2)
plot(1:461,KM03.ridge$KM_pval)
abline(h=0.05,col='red')
axis(1,1:461,round(KM03.ridge$prob_cutoff,digits=2),line=2)
plot(1:461,KM03.elastic$KM_pval)
abline(h=0.05,col='red')
axis(1,1:461,round(KM03.elastic$prob_cutoff,digits=2),line=2)

plot(KM03.random01$prob_cutoff,KM03.random01$KM_pval)
abline(h=0.05,col='red')
plot(KM03.random02$prob_cutoff,KM03.random02$KM_pval)
abline(h=0.05,col='red')
plot(KM03.random03$prob_cutoff,KM03.random03$KM_pval)
abline(h=0.05,col='red')
dev.off()



#save.image('TP53GEP_other_regression_models_APR08.RData')

######################################################################################
######################################################################################
#TCGA_clinical.rna.prep.178<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/TCGA/TCGA_preprocessed/clinical_rna_prep_178.rds')
#use new version of clinical data from GDC
#TCGA_clinical.rna.GDC.178<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/TCGA/TCGA_preprocessed/TCGA_clinical_GDC_rna_178_final.rds')
TCGA_clinical.rna.GDC.178
#TCGA_GEP.Data.corrected
tcga.y<-ifelse(TCGA_GEP.Data.corrected$TP53_mut_stat=='TP53_MUT',1,0)
prob.lasso.res.tcga #lasso identified prob
prob.ridge.res.tcga #ridge -the best performance 
#prob.elastic.res.tcga#elastic identified prob

TCGA_group<-cbind(TCGA_GEP.Data.corrected[,21619:21621],prob.lasso.res.tcga,prob.ridge.res.tcga)
TCGA_group$sampleID<-gsub('.03A.*|.03B.*','',TCGA_group$sampleID)
colnames(TCGA_group)<-c("patient_RNA_id","group","TP53_mut_stat","lasso.s1","ridge.s1")
TCGA_group$TP53_stat<-ifelse(TCGA_group$TP53_mut_stat=='TP53_MUT',1,0)
# TCGA_clinical.final<-right_join(TCGA_group,
#                                 TCGA_clinical.rna.prep.178,by = 'patient_RNA_id')
TCGA_clinical.final<-right_join(TCGA_group,
                                TCGA_clinical.rna.GDC.178,by = 'patient_RNA_id')
# 
# TCGA_clinical.final$OS_m<-ifelse(TCGA_clinical.final$`Overall Survival (Months)`=='[Not Available]',NA,
#                                  as.numeric(TCGA_clinical.final$`Overall Survival (Months)`))
# TCGA_clinical.final$OS_d<-as.numeric(TCGA_clinical.final$OS_m)*(365.24/12)
# TCGA_clinical.final$OS_status<-ifelse(TCGA_clinical.final$`Overall Survival Status`=='1:DECEASED',1,0)

# OS_RNA_df02<-data.frame('patient'=gsub('.*_','',BEAT_ct_clinc.RNA_avail.sample_patient.comb02$`#Patient Identifier`),
#                         'OS_m'=as.numeric(BEAT_ct_clinc.RNA_avail.sample_patient.comb02$`Overall Survival (Months)`),
#                         'OS_d'=floor(as.numeric(BEAT_ct_clinc.RNA_avail.sample_patient.comb02$`Overall Survival (Months)`)*365.24/12),
#                         'OS_status'=BEAT_ct_clinc.RNA_avail.sample_patient.comb02$`Overall Survival Status`,
#                         'TP53'=BEAT_ct_clinc.RNA_avail.sample_patient.comb02$`TP53 Mutation`)

dat.tcga<-TCGA_clinical.final
# dat.full$OS_status<-ifelse(dat.full$vitalStatus=='Dead',1,0)
# dat.full$TP53.z<-ifelse(dat.full$TP53.x=='TP53_MUT',1,0)
dat.tcga.obj<-Surv(time = dat.tcga$OS_d,
                   event = dat.tcga$OS_status)
mut_ixx.tcga<-which(dat.tcga$TP53_mut_stat=='TP53_MUT')

#check for the original prediction and prediction + TP53mut status
ix04.01<-order(dat.tcga$lasso.s1,decreasing = T)
ix04.02<-order(dat.tcga$ridge.s1,decreasing = T)

mut_ixx.tcga<-which(dat.tcga$TP53_mut_stat=='TP53_MUT')
set.seed(12345)
random_ix04.01<-sample(dim(dat.tcga)[1],dim(dat.tcga)[1])
random_ix04.02<-sample(dim(dat.tcga)[1],dim(dat.tcga)[1])
random_ix04.03<-sample(dim(dat.tcga)[1],dim(dat.tcga)[1])
set.seed(123)
ix04.random01<-c(sample(mut_ixx.tcga),random_ix04.01[!random_ix04.01 %in% mut_ixx.tcga])
set.seed(1235)
ix04.random02<-c(sample(mut_ixx.tcga),random_ix04.02[!random_ix04.02 %in% mut_ixx.tcga])
set.seed(1237)
ix04.random03<-c(sample(mut_ixx.tcga),random_ix04.03[!random_ix04.03 %in% mut_ixx.tcga])

KM04.01.df<-data.frame('prob_cutoff'=dat.tcga[ix04.01,]$lasso.s1,
                       'KM_pval'=rep(0,dim(dat.tcga)[1]),
                       'LabId'=dat.tcga[ix04.01,]$patient_RNA_id)

KM04.02.df<-data.frame('prob_cutoff'=dat.tcga[ix04.02,]$ridge.s1,
                       'KM_pval'=rep(0,dim(dat.tcga)[1]),
                       'LabId'=dat.tcga[ix04.02,]$patient_RNA_id)
##create random 
KM04.random01<-data.frame('prob_cutoff'=1:178,
                          'KM_pval'=rep(0,dim(dat.tcga)[1]))
KM04.random02<-data.frame('prob_cutoff'=1:178,
                          'KM_pval'=rep(0,dim(dat.tcga)[1]))
KM04.random03<-data.frame('prob_cutoff'=1:178,
                          'KM_pval'=rep(0,dim(dat.tcga)[1]))
for(i in 1:length(ix04.01)){
  ixx01<-ix04.01[1:i]
  ixx02<-ix04.02[1:i]
  
  ixx.random01<-ix04.random01[1:i]
  ixx.random02<-ix04.random02[1:i]
  ixx.random03<-ix04.random03[1:i]
  #print(ixx)
  
  mut_like_status<-rep(0,dim(dat.tcga)[1])
  mut_like_status.01<-mut_like_status
  mut_like_status.02<-mut_like_status
  mut_like_status.random01<-mut_like_status
  mut_like_status.random02<-mut_like_status
  mut_like_status.random03<-mut_like_status
  
  mut_like_status.01[ixx01]<-1
  mut_like_status.02[ixx02]<-1
  mut_like_status.random01[ixx.random01]<-1
  mut_like_status.random02[ixx.random02]<-1
  mut_like_status.random03[ixx.random03]<-1
  
  dat.tcga$TP53.cutoff.01<-mut_like_status.01
  dat.tcga$TP53.cutoff.02<-mut_like_status.02
  dat.tcga$TP53.random01<-mut_like_status.random01
  dat.tcga$TP53.random02<-mut_like_status.random02
  dat.tcga$TP53.random03<-mut_like_status.random03
  
  fit.km.01<-survfit(dat.tcga.obj ~TP53.cutoff.01, data=dat.tcga)
  fit.km.02<-survfit(dat.tcga.obj ~TP53.cutoff.02, data=dat.tcga)
  
  fit.km.random01<-survfit(dat.tcga.obj ~TP53.random01, data=dat.tcga)
  fit.km.random02<-survfit(dat.tcga.obj ~TP53.random02, data=dat.tcga)
  fit.km.random03<-survfit(dat.tcga.obj ~TP53.random03, data=dat.tcga)
  
  KM04.01.df[i,2]<-surv_pvalue(fit.km.01)[2]
  KM04.02.df[i,2]<-surv_pvalue(fit.km.02)[2]
  
  KM04.random01[i,2]<-surv_pvalue(fit.km.random01)[2]
  KM04.random02[i,2]<-surv_pvalue(fit.km.random02)[2]
  KM04.random03[i,2]<-surv_pvalue(fit.km.random03)[2]
}

pdf('TCGA_regressions_cutoff_MTWT_KM_GDCclinic.pdf')
plot(1:178,KM04.01.df$KM_pval)
abline(h=0.05,col='red')
axis(1,1:178,round(KM04.01.df$prob_cutoff,digits=2),line=2)
plot(1:178,KM04.02.df$KM_pval)
abline(h=0.05,col='red')
axis(1,1:178,round(KM04.02.df$prob_cutoff,digits=2),line=2)

plot(1:178,KM04.random01$KM_pval)
abline(h=0.05,col='red')
axis(1,1:178,round(KM04.random01$prob_cutoff,digits=2),line=2)
plot(1:178,KM04.random02$KM_pval)
abline(h=0.05,col='red')
axis(1,1:178,round(KM04.random02$prob_cutoff,digits=2),line=2)
plot(1:178,KM04.random03$KM_pval)
abline(h=0.05,col='red')
axis(1,1:178,round(KM04.random03$prob_cutoff,digits=2),line=2)
dev.off()

#Survival analysis without TP53MUT cases in TCGA
#mut_ixx.tcga<-which(dat.tcga$TP53_mut_stat=='TP53_MUT')
dat.tcga.noTP53mut<-dat.tcga[-mut_ixx.tcga,]
dat.tcga.noTP53mut.obj<-Surv(time = dat.tcga.noTP53mut$OS_d,
                             event = dat.tcga.noTP53mut$OS_status)

#check for the original prediction and prediction + TP53mut status
ix05.01<-order(dat.tcga.noTP53mut$lasso.s1,decreasing = T)
ix05.02<-order(dat.tcga.noTP53mut$ridge.s1,decreasing = T)

KM05.01.df<-data.frame('prob_cutoff'=dat.tcga.noTP53mut[ix05.01,]$lasso.s1,
                       'KM_pval'=rep(0,dim(dat.tcga.noTP53mut)[1]),
                       'LabId'=dat.tcga.noTP53mut[ix05.01,]$patient_RNA_id)

KM05.02.df<-data.frame('prob_cutoff'=dat.tcga.noTP53mut[ix05.02,]$ridge.s1,
                       'KM_pval'=rep(0,dim(dat.tcga.noTP53mut)[1]),
                       'LabId'=dat.tcga.noTP53mut[ix05.02,]$patient_RNA_id)

for(i in 1:length(ix05.01)){
  ixx01<-ix05.01[1:i]
  ixx02<-ix05.02[1:i]
  
  mut_like_status<-rep(0,dim(dat.tcga.noTP53mut)[1])
  mut_like_status.01<-mut_like_status
  mut_like_status.02<-mut_like_status
  
  mut_like_status.01[ixx01]<-1
  mut_like_status.02[ixx02]<-1
  
  dat.tcga.noTP53mut$TP53.cutoff.01<-mut_like_status.01
  dat.tcga.noTP53mut$TP53.cutoff.02<-mut_like_status.02
  
  fit.km.01<-survfit(dat.tcga.noTP53mut.obj ~TP53.cutoff.01, data=dat.tcga.noTP53mut)
  fit.km.02<-survfit(dat.tcga.noTP53mut.obj ~TP53.cutoff.02, data=dat.tcga.noTP53mut)
  
  KM05.01.df[i,2]<-surv_pvalue(fit.km.01)[2]
  KM05.02.df[i,2]<-surv_pvalue(fit.km.02)[2]
}
#pdf('TCGA_regressions_cutoff_WT_KM_GDCclinc.pdf')
plot(1:163,KM05.01.df$KM_pval)
abline(h=0.05,col='red')
axis(1,1:163,round(KM05.01.df$prob_cutoff,digits=2),line=2)
plot(1:163,KM05.02.df$KM_pval)
abline(h=0.05,col='red')
axis(1,1:163,round(KM05.02.df$prob_cutoff,digits=2),line=2)
#dev.off()

#recheck TCGA cases

#setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis/LASSO/GEP")


#2022.05.19
#test for the new KM screening: set comparison group A from top 0% to 50% of the data while
#holding the bottom 50% of the score group as the comparison group B (every screening case)

#1.test TCGA lasso and ridge first 
#2. BEAT AML test for LASSO and Ridge

#1.test TCGA lasso and ridge first 

dat.tcga.noTP53mut<-dat.tcga[-mut_ixx.tcga,]
dat.tcga.noTP53mut.obj<-Surv(time = dat.tcga.noTP53mut$OS_d,
                             event = dat.tcga.noTP53mut$OS_status)

#check for the original prediction and prediction + TP53mut status
ix05.01<-order(dat.tcga.noTP53mut$lasso.s1,decreasing = T)
ix05.02<-order(dat.tcga.noTP53mut$ridge.s1,decreasing = T)

KM06.01.df<-data.frame('prob_cutoff'=dat.tcga.noTP53mut[ix05.01,]$lasso.s1,
                       'KM_pval'=rep(NA,dim(dat.tcga.noTP53mut)[1]),
                       'LabId'=dat.tcga.noTP53mut[ix05.01,]$patient_RNA_id)

KM06.02.df<-data.frame('prob_cutoff'=dat.tcga.noTP53mut[ix05.02,]$ridge.s1,
                       'KM_pval'=rep(NA,dim(dat.tcga.noTP53mut)[1]),
                       'LabId'=dat.tcga.noTP53mut[ix05.02,]$patient_RNA_id)

for(i in 1:floor(length(ix05.01)*0.5)){# walk only until the top 50% of the data: 1 to 81
  ixx01<-ix05.01[1:i]
  ixx02<-ix05.02[1:i]
  ixx01.base<-ix05.01[(floor(length(ix05.01)*0.65)+1):length(ix05.01)]
  ixx02.base<-ix05.02[(floor(length(ix05.01)*0.65)+1):length(ix05.01)]
  
  mut_like_status<-rep(0,dim(dat.tcga.noTP53mut)[1])
  mut_like_status.01<-mut_like_status
  mut_like_status.02<-mut_like_status
  
  mut_like_status.01[ixx01]<-1
  mut_like_status.02[ixx02]<-1
  
  dat.tcga.noTP53mut$TP53.cutoff.01<-mut_like_status.01
  dat.tcga.noTP53mut$TP53.cutoff.02<-mut_like_status.02
  
  dat.tcga.noTP53mut.obj.01<-Surv(time = dat.tcga.noTP53mut$OS_d[c(ixx01,ixx01.base)],
                                  event = dat.tcga.noTP53mut$OS_status[c(ixx01,ixx01.base)])
  dat.tcga.noTP53mut.obj.02<-Surv(time = dat.tcga.noTP53mut$OS_d[c(ixx02,ixx02.base)],
                                  event = dat.tcga.noTP53mut$OS_status[c(ixx02,ixx02.base)])
  fit.km.01<-survfit(dat.tcga.noTP53mut.obj.01 ~TP53.cutoff.01, data=dat.tcga.noTP53mut[c(ixx01,ixx01.base),])
  fit.km.02<-survfit(dat.tcga.noTP53mut.obj.02 ~TP53.cutoff.02, data=dat.tcga.noTP53mut[c(ixx02,ixx02.base),])
  
  KM06.01.df[i,2]<-surv_pvalue(fit.km.01)[2]
  KM06.02.df[i,2]<-surv_pvalue(fit.km.02)[2]
}
pdf('TCGA_regressions_cutoff_WT_KM_topbottom.pdf')
plot(1:163,KM06.01.df$KM_pval)
abline(h=0.05,col='red')
axis(1,1:163,round(KM05.01.df$prob_cutoff,digits=2),line=2)
plot(1:163,KM06.02.df$KM_pval)
abline(h=0.05,col='red')
axis(1,1:163,round(KM06.02.df$prob_cutoff,digits=2),line=2)
dev.off()

#2. test for BEAT AML case with fixed bottom 50% (onl WT and MT)
ix01.lasso<-order(dat.sub01.wt$lasso.s1,decreasing = T)
ix01.ridge<-order(dat.sub01.wt$ridge.s1,decreasing = T)
ix01.elastic<-order(dat.sub01.wt$elastic.s1,decreasing = T)


KM07.lasso<-data.frame('prob_cutoff'=dat.sub01.wt[ix01.lasso,]$lasso.s1,
                       'KM_pval'=rep(NA,dim(dat.sub01.wt)[1]),
                       'LabId'=dat.sub01.wt[ix01.lasso,]$LabId)
KM07.ridge<-data.frame('prob_cutoff'=dat.sub01.wt[ix01.ridge,]$ridge.s1,
                       'KM_pval'=rep(NA,dim(dat.sub01.wt)[1]),
                       'LabId'=dat.sub01.wt[ix01.ridge,]$LabId)
KM07.elastic<-data.frame('prob_cutoff'=dat.sub01.wt[ix01.elastic,]$elastic.s1,
                         'KM_pval'=rep(NA,dim(dat.sub01.wt)[1]),
                         'LabId'=dat.sub01.wt[ix01.elastic,]$LabId)

top_cutoff<-0.4
bottom_cutoff<-0.5 #bottom 
for(i in 1:floor(length(ix01.elastic)*top_cutoff)){
  #for(i in 1:5){
  
  ixx.lasso<-ix01.lasso[1:i]
  ixx.ridge<-ix01.ridge[1:i]
  ixx.elastic<-ix01.elastic[1:i]
  ixx.lasso.base<-ix01.lasso[(floor(length(ix01.lasso)*bottom_cutoff)+1):length(ix01.lasso)]
  ixx.ridge.base<-ix01.ridge[(floor(length(ix01.ridge)*bottom_cutoff)+1):length(ix01.ridge)]
  ixx.elastic.base<-ix01.elastic[(floor(length(ix01.elastic)*bottom_cutoff)+1):length(ix01.elastic)]
  
  mut_like_status<-rep(0,dim(dat.sub01.wt)[1])
  mut_like_status.lasso<-mut_like_status
  mut_like_status.ridge<-mut_like_status
  mut_like_status.elastic<-mut_like_status
  
  mut_like_status.lasso[ixx.lasso]<-1
  mut_like_status.ridge[ixx.ridge]<-1
  mut_like_status.elastic[ixx.elastic]<-1
  dat.sub01.wt$TP53.lasso<-mut_like_status.lasso
  dat.sub01.wt$TP53.ridge<-mut_like_status.ridge
  dat.sub01.wt$TP53.elastic<-mut_like_status.elastic
  
  dat.wt.obj.lasso<-Surv(time = dat.sub01.wt$overallSurvival[c(ixx.lasso,ixx.lasso.base)],
                         event = dat.sub01.wt$OS_status[c(ixx.lasso,ixx.lasso.base)])
  dat.wt.obj.ridge<-Surv(time = dat.sub01.wt$overallSurvival[c(ixx.ridge,ixx.ridge.base)],
                         event = dat.sub01.wt$OS_status[c(ixx.ridge,ixx.ridge.base)])
  dat.wt.obj.elastic<-Surv(time = dat.sub01.wt$overallSurvival[c(ixx.elastic,ixx.elastic.base)],
                           event = dat.sub01.wt$OS_status[c(ixx.elastic,ixx.elastic.base)])
  
  fit.km.lasso<-survfit(dat.wt.obj.lasso ~TP53.lasso, data=dat.sub01.wt[c(ixx.lasso,ixx.lasso.base),])
  fit.km.ridge<-survfit(dat.wt.obj.ridge ~TP53.ridge, data=dat.sub01.wt[c(ixx.ridge,ixx.ridge.base),])
  fit.km.elastic<-survfit(dat.wt.obj.elastic ~TP53.elastic, data=dat.sub01.wt[c(ixx.elastic,ixx.elastic.base),])
  
  KM07.lasso[i,2]<-surv_pvalue(fit.km.lasso)[2]
  KM07.ridge[i,2]<-surv_pvalue(fit.km.ridge)[2]
  KM07.elastic[i,2]<-surv_pvalue(fit.km.elastic)[2]
  
}
#pdf('BEAT_regressions_cutoff_WTunknown_KM_topbottom.pdf')
plot(1:377,KM07.lasso$KM_pval)
abline(h=0.05,col='red')
axis(1,1:377,round(KM07.lasso$prob_cutoff,digits=2),line=2)
plot(1:377,KM07.ridge$KM_pval)
abline(h=0.05,col='red')
axis(1,1:377,round(KM07.ridge$prob_cutoff,digits=2),line=2)
plot(1:377,KM07.elastic$KM_pval)
abline(h=0.05,col='red')
axis(1,1:377,round(KM07.elastic$prob_cutoff,digits=2),line=2)
#dev.off()

#check KM02.ridge.46id to the KM01.ridge cutoff
KM02.ridge.46id<-KM02.ridge$LabId[1:46]
which(KM01.ridge$LabId %in% KM02.ridge.46id)

##############################################################

#fix Mut-like cases using KM results -lasso and elastic only:
#1.WT-MUT-like vs.WT-like 
KM01.lasso#46 23(0.0647)- skip
KM01.ridge#3,4,5,6,7
KM01.ridge.42id<-KM01.ridge$LabId[1:42]

#2.unknown-MUT-like vs. unknown-WT-like +WT
KM02.lasso#1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19,20,21,23,24,25,26,27,28,30,31,32,33,34,35,38,420,421,422,423
KM02.ridge#1,2,3,19,20,21,27,28,29,30,31,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61
KM02.lasso.id<-KM02.lasso$LabId[1:46]

KM07.lasso#1,2,3,4,5,6,7,8,9,10,11,12,13,19,20,23,24,25,26,27,28,29,30
KM07.ridge#1,2,3,4,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61
#KM07.elastic#1,2,3,4,20

setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP")
#save.image('GLM_CV_GEP003.Rdata')#one that before using cv$fit.preval for train data score
#save.image('GLM_CV_GEP004.Rdata')#one that after using cv$fit.preval for train data score(2022.05.31)
#load('GLM_CV_GEP003.Rdata')
load('GLM_CV_GEP004.Rdata')
gc()
#3.TP53MUT+unknown-MUT-like vs. unknown-WT-like +WT - skip this for now
###########################
###########################
BEAT_LASSO_final_res.clinical.n
cls.rep<-(BEAT_LASSO_final_res.clinical.n$TP53_mut_stat)
#make class for the screening
#MT vs. WT
KM01.ridge.42id<-KM01.ridge$LabId[1:42]

#MT vs. WT+unknown
KM02.lasso.id
KM02.ridge.id
KM02.ridge.46.id

KM07.lasso.id
KM07.ridge.id.30
KM07.ridge.id.61
KM07.elastic.id

#GLM_id<-c()
GLM_id<-KM01.ridge.42id############

cls.rep<-(BEAT_LASSO_final_res.clinical.n$TP53_mut_stat)
cls.GLM<-cls.rep;cls.GLM[which(BEAT_LASSO_final_res.clinical.n$LabId %in% GLM_id)]<-'TP53_WTunk_MT_like';
cls.GLM<-ifelse(cls.GLM=='unknown','TP53_unk_WT_like',cls.GLM)
BEAT_LASSO_final_res.clinical.n$cls.GLM<-cls.GLM

#quickly check class and OS in the MUT-like WT
#MTWTunk 
ss<-right_join(KM02.ridge,BEAT_LASSO_final_res.clinical.n,by='LabId')
View(ss[,c(1,2,63,8,17)])
#only MT and WT
sss<-ss[which(ss$TP53_mut_stat=='TP53_MUT' | ss$TP53_mut_stat=='TP53_WT'),c(1,2,63,8,17)]
View(sss)
####################################################################################
#survival plot
#save.image('GLM_CV_GEP005.Rdata')#one that before using cv$fit.preval for train data score
load('GLM_CV_GEP005.Rdata')
#do chi-square/pairwise-fisher exact testing for one vairable for now
#how to extract count data?make a separate data to do counting data

setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP")
#save.image('GLM_CV_GEP005_Mut_stat.Rdata')#one that before using cv$fit.preval for train data score
load('GLM_CV_GEP005_Mut_stat.Rdata')
gc()


################################################
#get TCGA clinical data - not many data to look at
TCGA_clinical.final
LSC17_score.tcga$TCGA.Patient.ID.y<-sub('\\.03.*','',LSC17_score.tcga$LabId)
LSC17_score.comb.tcga$TCGA.Patient.ID.y<-sub('\\.03.*','',LSC17_score.comb.tcga$LabId)

TCGA_clinical.final.n<-right_join(TCGA_clinical.final,LSC17_score.tcga,by='TCGA.Patient.ID.y')
TCGA_clinical.final.n<-right_join(TCGA_clinical.final.n,LSC17_score.comb.tcga,by='TCGA.Patient.ID.y')

#fix Mut-like cases using KM results -lasso and elastic only:
#1.unknown-MUT-like+TP53 vs.unknown-WT-like 
KM04.01.df#lasso
KM04.02.df#ridge
# KM01.lasso.id<-KM01.lasso$LabId[1:18]
# KM01.elastic.id<-KM01.elastic$LabId[1:20]

#2.unknown-MUT-like vs. unknown-WT-like
KM05.01.df#lasso only 1
KM05.02.df#elastic 8 unknown(other-5 and 6)
KM05.lasso.id<-KM05.01.df$LabId[1]
KM05.elastic.id<-KM05.02.df$LabId[1:8]
#GLM_id<-c()
GLM_id.tcga<-KM05.elastic.id############
cls.rep.tcga<-(TCGA_clinical.final.n$TP53_mut_stat)
#View(cbind(prob.elastic.res.tcga,tcga.y))

###################################################################################
#have not updated to newer version yet
#2022.05.19
#New classification analysis: # this is hierarchcical model
#1.additional GLM with ridge identified TP53MUT+TP53MUT-like WT(based on the KM screening) vs. WT
#to identify TP53MUT-like WT in TCGA 
#2. additional GLM with ridge identified TP53MUT-like WT(based on the KM screening) vs. WT
#to identify TP53MUT-like WT in TCGA (Wt only)
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP")
#save.image('GLM_CV_GEP003.Rdata')
#load('GLM_CV_GEP003.Rdata')
load('GLM_CV_GEP004.Rdata')
gc()

source("C:/Users/yklee/Desktop/Sachs Lab/R_project/TP53/TP53GEP_functions.r")
packages<-c('edgeR','limma','dplyr','tidyr','sva','tidyverse',
            'ggplot2','glmnet','caret','pROC','survival','survminer')

for(p in packages){
  library(p,character.only = T)
}


TP53MT_WT_ix<-which(BEAT_GEP.Data.corrected$TP53_mut_stat != 'unknown')
TP53_MT_nMT<-rep(0,dim(BEAT_GEP.Data.corrected)[1])
TP53_MT_nMT[which(BEAT_GEP.Data.corrected$TP53_mut_stat=='TP53_MUT')]<-1
TP53_MT_WT<-TP53_MT_nMT[TP53MT_WT_ix]
TP53_WT_ix<-which(BEAT_GEP.Data.corrected$TP53_mut_stat == 'TP53_WT')
#MT vs. WT
KM01.ridge.42id<-KM01.ridge$LabId[1:42]

#2022.06.07
#saveRDS(KM01.ridge,'class ids/BEAT_WT_KM01_ridge.rds')
#saveRDS(KM02.ridge,'class ids/BEAT_WTUnk_KM02_ridge.rds')

##########################################################################################
#2. additional GLM with ridge identified TP53MUT-like WT(based on the KM screening) vs. WT
#to identify TP53MUT-like WT in TCGA (Wt only)
#use above MUT-like WT + MUT as MUT classes to retrain BEAT AML data
class_id<-KM01.ridge.42id


##1.4. run LASSO train/test  
BEAT_GEP.Data.corrected$TP53_mutlike_stat<-rep(0,dim(BEAT_GEP.Data.corrected)[1])
BEAT_GEP.Data.corrected$TP53_mutlike_stat[BEAT_GEP.Data.corrected$sampleID %in% class_id]<-1
TP53_group_sub.y<-as.factor(BEAT_GEP.Data.corrected$TP53_mutlike_stat[TP53_WT_ix])

#subset data for MT_WT and unknown and with all vs. DEGs
#1.all genes
#get only WT
#BEAT_GEP_TP53MT_WT.all<-data.matrix(BEAT_GEP.Data.corrected[TP53MT_WT_ix,-(21619:21622)])
BEAT_GEP_TP53WT.all<-data.matrix(BEAT_GEP.Data.corrected[TP53_WT_ix,-(21619:21622)])
BEAT_GEP_TP53nMT.all<-data.matrix(BEAT_GEP.Data.corrected[-TP53MT_WT_ix,-(21619:21622)])

setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP")



#set.seed(9011230)
#set.seed(90112312)
#set.seed(9011233)
set.seed(11233)
#set.seed(901123)
trn.idx<-TP53_group_sub.y %>% createDataPartition(p=0.6, list = FALSE)
trn.data<-data.matrix(BEAT_GEP_TP53WT.all[trn.idx,])
trn.y<-(TP53_group_sub.y[trn.idx])
test.data<-data.matrix(BEAT_GEP_TP53WT.all[-trn.idx,])
test.y<-(TP53_group_sub.y[-trn.idx])

#set 10-fold cv for the model comparison
set.seed(2843)
#set.seed(284)
#set.seed(283)
flds <- createFolds(trn.y, k = 10, list = TRUE, returnTrain = FALSE)
foldids = rep(1,length(trn.y))
foldids[flds$Fold02] = 2
foldids[flds$Fold03] = 3
foldids[flds$Fold04] = 4
foldids[flds$Fold05] = 5
foldids[flds$Fold06] = 6
foldids[flds$Fold07] = 7
foldids[flds$Fold08] = 8
foldids[flds$Fold09] = 9
foldids[flds$Fold10] = 10
#try ridge regression and elastic net 
cv.lasso<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=1,keep = TRUE,
                    type.measure = 'class',foldid=foldids)
cv.ridge<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=0,keep = TRUE,
                    type.measure = 'class',foldid=foldids)
cv.elastic<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=0.5,keep = TRUE,
                      type.measure = 'class',foldid=foldids)
# cv.lasso<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=1,keep = TRUE,
#                     type.measure = 'class')
# cv.ridge<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=0,keep = TRUE,
#                     type.measure = 'class')
# cv.elastic<-cv.glmnet(trn.data,trn.y,family='binomial',alpha=0.5,keep = TRUE,
#                       type.measure = 'class')

# plot(cv.lasso)
# plot(cv.ridge)
# plot(cv.elastic)
model.lasso <- glmnet(trn.data, trn.y, alpha = 1, family = "binomial",
                      lambda = cv.lasso$lambda.1se)
model.ridge <- glmnet(trn.data, trn.y, alpha = 0, family = "binomial",
                      lambda = cv.ridge$lambda.1se)
model.elastic <- glmnet(trn.data, trn.y, alpha = 0.5, family = "binomial",
                        lambda = cv.elastic$lambda.1se)

# 
# model.lasso <- glmnet(trn.data, trn.y, alpha = 1, family = "binomial",
#                       lambda = cv.lasso$lambda.min)
# model.ridge <- glmnet(trn.data, trn.y, alpha = 0, family = "binomial",
#                       lambda = cv.ridge$lambda.min)
# model.elastic <- glmnet(trn.data, trn.y, alpha = 0.5, family = "binomial",
#                         lambda = cv.elastic$lambda.min)

prob.lasso.res <- model.lasso %>% predict(newx = trn.data,s=cv.lasso$lambda.1se,type='response')
prob.ridge.res <- model.ridge %>% predict(newx = trn.data,s=cv.ridge$lambda.1se,type='response')
prob.elastic.res <- model.elastic %>% predict(newx = trn.data,s=cv.elastic$lambda.1se,type='response')


prob.lasso.res.test <- model.lasso %>% predict(newx = test.data,type='response')
prob.ridge.res.test <- model.ridge %>% predict(newx = test.data,type='response')
prob.elastic.res.test <- model.elastic %>% predict(newx = test.data,type='response')


prob.lasso.res.tcga <- model.lasso %>% predict(newx = test.tcga,s=cv.lasso$lambda.1se,type='response')
prob.ridge.res.tcga <- model.ridge %>% predict(newx = test.tcga,s=cv.ridge$lambda.1se,type='response')
prob.elastic.res.tcga <- model.elastic %>% predict(newx = test.tcga,s=cv.elastic$lambda.1se,type='response')

out.performance<-get_AUC_summary(trn.y,test.y,tcga.y,prob.lasso.res,prob.ridge.res,prob.elastic.res,
                                 prob.lasso.res.test,prob.ridge.res.test,prob.elastic.res.test,
                                 prob.lasso.res.tcga,prob.ridge.res.tcga,prob.elastic.res.tcga)
#write.csv(out.performance,'out.performance.csv')


roc_trn.ridge<-pROC_summary(trn.y,prob.ridge.res)
PR_trn.ridge<-Prec_Recall_summary(trn.y,prob.ridge.res)

roc_test.ridge<-pROC_summary(test.y,prob.ridge.res.test)
PR_test.ridge<-Prec_Recall_summary(test.y,prob.ridge.res.test)

roc_test.ridge<-pROC_summary(tcga.y,prob.ridge.res.tcga)
PR_test.ridge<-Prec_Recall_summary(tcga.y,prob.ridge.res.tcga)

# tmp_coeffs<-coef(model.ridge)
# length(coef(model.ridge))
# predict(model.ridge,type="coef")
# coef.df<-data.frame(name = tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
# View(varImp(model.ridge,lambda = model.ridge$lambda))
# 
# plot(model.ridge,xvar='lambda')
# plot(cv.ridge)


#########################################################################
#see if the models work in TCGA data first
tcga.y<-ifelse(TCGA_GEP.Data.corrected$TP53_mut_stat=='TP53_MUT',1,0)
test.tcga<-data.matrix(TCGA_GEP.Data.corrected[,-(21619:21621)])

prob.lasso.res.tcga <- model.lasso %>% predict(newx = test.tcga,s=cv.lasso$lambda.1se,type='response')
prob.ridge.res.tcga <- model.ridge %>% predict(newx = test.tcga,s=cv.ridge$lambda.1se,type='response')
prob.elastic.res.tcga <- model.elastic %>% predict(newx = test.tcga,s=cv.elastic$lambda.1se,type='response')

# hist(prob.lasso.res.tcga)
# hist(prob.ridge.res.tcga)
# hist(prob.elastic.res.tcga)

roc_tcga.lasso<-pROC_summary(tcga.y,prob.lasso.res.tcga)
PR_tcga.lasso<-Prec_Recall_summary(tcga.y,prob.lasso.res.tcga)

roc_tcga.ridge<-pROC_summary(tcga.y,prob.ridge.res.tcga)
PR_tcga.ridge<-Prec_Recall_summary(tcga.y,prob.ridge.res.tcga)

roc_tcga.elastic<-pROC_summary(tcga.y,prob.elastic.res.tcga)
PR_tcga.elastic<-Prec_Recall_summary(tcga.y,prob.elastic.res.tcga)

###################################################################
#check on TCGA KM screening


TCGA_group<-cbind(TCGA_GEP.Data.corrected[,21619:21621],prob.lasso.res.tcga,prob.ridge.res.tcga,prob.elastic.res.tcga)
TCGA_group$sampleID<-gsub('.03A.*|.03B.*','',TCGA_group$sampleID)
colnames(TCGA_group)<-c("patient_RNA_id","group","TP53_mut_stat","lasso.s1","ridge.s1","elastic.s1")
TCGA_group$TP53_stat<-ifelse(TCGA_group$TP53_mut_stat=='TP53_MUT',1,0)
# TCGA_clinical.final<-right_join(TCGA_group,
#                                 TCGA_clinical.rna.prep.178,by = 'patient_RNA_id')
TCGA_clinical.final<-right_join(TCGA_group,
                                TCGA_clinical.rna.GDC.178,by = 'patient_RNA_id')

dat.tcga<-TCGA_clinical.final
# dat.full$OS_status<-ifelse(dat.full$vitalStatus=='Dead',1,0)
# dat.full$TP53.z<-ifelse(dat.full$TP53.x=='TP53_MUT',1,0)
dat.tcga.obj<-Surv(time = dat.tcga$OS_d,
                   event = dat.tcga$OS_status)

#check for the original prediction and prediction + TP53mut status
mut_ixx.tcga<-which(dat.tcga$TP53_mut_stat=='TP53_MUT')

dat.tcga.noTP53mut<-dat.tcga[-mut_ixx.tcga,]
dat.tcga.noTP53mut.obj<-Surv(time = dat.tcga.noTP53mut$OS_d,
                             event = dat.tcga.noTP53mut$OS_status)

#check for the original prediction and prediction + TP53mut status
ix05.01<-order(dat.tcga.noTP53mut$lasso.s1,decreasing = T)
ix05.02<-order(dat.tcga.noTP53mut$ridge.s1,decreasing = T)
ix05.03<-order(dat.tcga.noTP53mut$elastic.s1,decreasing = T)

KM08.01.df<-data.frame('prob_cutoff'=dat.tcga.noTP53mut[ix05.01,]$lasso.s1,
                       'KM_pval'=rep(0,dim(dat.tcga.noTP53mut)[1]),
                       'LabId'=dat.tcga.noTP53mut[ix05.01,]$patient_RNA_id)

KM08.02.df<-data.frame('prob_cutoff'=dat.tcga.noTP53mut[ix05.02,]$ridge.s1,
                       'KM_pval'=rep(0,dim(dat.tcga.noTP53mut)[1]),
                       'LabId'=dat.tcga.noTP53mut[ix05.02,]$patient_RNA_id)
KM08.03.df<-data.frame('prob_cutoff'=dat.tcga.noTP53mut[ix05.03,]$elastic.s1,
                       'KM_pval'=rep(0,dim(dat.tcga.noTP53mut)[1]),
                       'LabId'=dat.tcga.noTP53mut[ix05.03,]$patient_RNA_id)

for(i in 1:length(ix05.01)){
  ixx01<-ix05.01[1:i]
  ixx02<-ix05.02[1:i]
  ixx03<-ix05.03[1:i]
  
  mut_like_status<-rep(0,dim(dat.tcga.noTP53mut)[1])
  mut_like_status.01<-mut_like_status
  mut_like_status.02<-mut_like_status
  mut_like_status.03<-mut_like_status
  
  
  mut_like_status.01[ixx01]<-1
  mut_like_status.02[ixx02]<-1
  mut_like_status.03[ixx03]<-1
  
  
  dat.tcga.noTP53mut$TP53.cutoff.01<-mut_like_status.01
  dat.tcga.noTP53mut$TP53.cutoff.02<-mut_like_status.02
  dat.tcga.noTP53mut$TP53.cutoff.03<-mut_like_status.03
  
  fit.km.01<-survfit(dat.tcga.noTP53mut.obj ~TP53.cutoff.01, data=dat.tcga.noTP53mut)
  fit.km.02<-survfit(dat.tcga.noTP53mut.obj ~TP53.cutoff.02, data=dat.tcga.noTP53mut)
  fit.km.03<-survfit(dat.tcga.noTP53mut.obj ~TP53.cutoff.03, data=dat.tcga.noTP53mut)
  
  
  KM08.01.df[i,2]<-surv_pvalue(fit.km.01)[2]
  KM08.02.df[i,2]<-surv_pvalue(fit.km.02)[2]
  KM08.03.df[i,2]<-surv_pvalue(fit.km.03)[2]
  
}
#pdf('TCGA_regressions_cutoff_WT_KM_GDCclinc.pdf')
plot(1:163,KM08.01.df$KM_pval)
abline(h=0.05,col='red')
axis(1,1:163,round(KM08.01.df$prob_cutoff,digits=2),line=2)
plot(1:163,KM08.02.df$KM_pval)
abline(h=0.05,col='red')
axis(1,1:163,round(KM08.02.df$prob_cutoff,digits=2),line=2)
plot(1:163,KM08.03.df$KM_pval)
abline(h=0.05,col='red')
axis(1,1:163,round(KM08.03.df$prob_cutoff,digits=2),line=2)
#dev.off()
#saveRDS(KM08.02.df,'class ids/TCGA_WTUnk_KM08.ridge.rds')


dat.tcga.noTP53mut<-dat.tcga[-mut_ixx.tcga,]
dat.tcga.noTP53mut.obj<-Surv(time = dat.tcga.noTP53mut$OS_d,
                             event = dat.tcga.noTP53mut$OS_status)

#check for the original prediction and prediction + TP53mut status
ix05.01<-order(dat.tcga.noTP53mut$lasso.s1,decreasing = T)
ix05.02<-order(dat.tcga.noTP53mut$ridge.s1,decreasing = T)
ix05.03<-order(dat.tcga.noTP53mut$elastic.s1,decreasing = T)

KM09.01.df<-data.frame('prob_cutoff'=dat.tcga.noTP53mut[ix05.01,]$lasso.s1,
                       'KM_pval'=rep(NA,dim(dat.tcga.noTP53mut)[1]),
                       'LabId'=dat.tcga.noTP53mut[ix05.01,]$patient_RNA_id)
KM09.02.df<-data.frame('prob_cutoff'=dat.tcga.noTP53mut[ix05.02,]$ridge.s1,
                       'KM_pval'=rep(NA,dim(dat.tcga.noTP53mut)[1]),
                       'LabId'=dat.tcga.noTP53mut[ix05.02,]$patient_RNA_id)
KM09.03.df<-data.frame('prob_cutoff'=dat.tcga.noTP53mut[ix05.03,]$elastic.s1,
                       'KM_pval'=rep(NA,dim(dat.tcga.noTP53mut)[1]),
                       'LabId'=dat.tcga.noTP53mut[ix05.03,]$patient_RNA_id)

top_cutoff<-0.4
bottom_cutoff<-0.75 #bottom 
for(i in 1:floor(length(ix05.01)*top_cutoff)){# walk only until the top 50% of the data: 1 to 81
  ixx01<-ix05.01[1:i]
  ixx02<-ix05.02[1:i]
  ixx03<-ix05.03[1:i]
  
  ixx01.base<-ix05.01[(floor(length(ix05.01)*bottom_cutoff)+1):length(ix05.01)]
  ixx02.base<-ix05.02[(floor(length(ix05.01)*bottom_cutoff)+1):length(ix05.01)]
  ixx03.base<-ix05.03[(floor(length(ix05.01)*bottom_cutoff)+1):length(ix05.01)]
  
  
  mut_like_status<-rep(0,dim(dat.tcga.noTP53mut)[1])
  mut_like_status.01<-mut_like_status
  mut_like_status.02<-mut_like_status
  mut_like_status.03<-mut_like_status
  
  
  mut_like_status.01[ixx01]<-1
  mut_like_status.02[ixx02]<-1
  mut_like_status.03[ixx02]<-1
  
  
  dat.tcga.noTP53mut$TP53.cutoff.01<-mut_like_status.01
  dat.tcga.noTP53mut$TP53.cutoff.02<-mut_like_status.02
  dat.tcga.noTP53mut$TP53.cutoff.03<-mut_like_status.03
  
  
  dat.tcga.noTP53mut.obj.01<-Surv(time = dat.tcga.noTP53mut$OS_d[c(ixx01,ixx01.base)],
                                  event = dat.tcga.noTP53mut$OS_status[c(ixx01,ixx01.base)])
  dat.tcga.noTP53mut.obj.02<-Surv(time = dat.tcga.noTP53mut$OS_d[c(ixx02,ixx02.base)],
                                  event = dat.tcga.noTP53mut$OS_status[c(ixx02,ixx02.base)])
  dat.tcga.noTP53mut.obj.03<-Surv(time = dat.tcga.noTP53mut$OS_d[c(ixx03,ixx03.base)],
                                  event = dat.tcga.noTP53mut$OS_status[c(ixx03,ixx03.base)])
  fit.km.01<-survfit(dat.tcga.noTP53mut.obj.01 ~TP53.cutoff.01, data=dat.tcga.noTP53mut[c(ixx01,ixx01.base),])
  fit.km.02<-survfit(dat.tcga.noTP53mut.obj.02 ~TP53.cutoff.02, data=dat.tcga.noTP53mut[c(ixx02,ixx02.base),])
  fit.km.03<-survfit(dat.tcga.noTP53mut.obj.03 ~TP53.cutoff.03, data=dat.tcga.noTP53mut[c(ixx03,ixx03.base),])
  
  KM09.01.df[i,2]<-surv_pvalue(fit.km.01)[2]
  KM09.02.df[i,2]<-surv_pvalue(fit.km.02)[2]
  KM09.03.df[i,2]<-surv_pvalue(fit.km.03)[2]
}
plot(1:163,KM09.01.df$KM_pval)
abline(h=0.05,col='red')
axis(1,1:163,round(KM09.01.df$prob_cutoff,digits=2),line=2)
plot(1:163,KM09.02.df$KM_pval)
abline(h=0.05,col='red')
axis(1,1:163,round(KM09.02.df$prob_cutoff,digits=2),line=2)
plot(1:163,KM09.03.df$KM_pval)
abline(h=0.05,col='red')
axis(1,1:163,round(KM09.03.df$prob_cutoff,digits=2),line=2)


setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP")

ixx<-ix05.02[1:23]
mut_like_status<-rep(0,dim(dat.tcga.noTP53mut)[1])
mut_like_status.02<-mut_like_status
mut_like_status.02[ixx]<-1
dat.tcga.noTP53mut$TP53.cutoff.02<-mut_like_status.02  
fit.km.02<-survfit(dat.tcga.noTP53mut.obj ~TP53.cutoff.02, data=dat.tcga.noTP53mut)
ggsurvplot(fit.km.02, data = dat.tcga.noTP53mut, pval = TRUE)
#sort(dat.tcga.noTP53mut$OS_d[dat.tcga.noTP53mut$TP53.cutoff.02==1])
#sort(dat.tcga.noTP53mut$OS_d[dat.tcga.noTP53mut$TP53.cutoff.02==0])


dat.tcga<-TCGA_clinical.final
dat.tcga.obj<-Surv(time = dat.tcga$OS_d,
                   event = dat.tcga$OS_status)
mut_ixx.tcga<-which(dat.tcga$TP53_mut_stat=='TP53_MUT')
wt_ixx.tcga<-which(dat.tcga$TP53_mut_stat=='TP53_WT')

ixx.id<-KM08.02.df$LabId[1:23]
mut_like_ixx.tcga<-which(dat.tcga$patient_RNA_id %in% ixx.id)

mut_like_status.03<-rep(0,dim(dat.tcga)[1])
mut_like_status.03[mut_ixx.tcga]<-'TP53_MUT';mut_like_status.03[wt_ixx.tcga]<-'TP53_WT';mut_like_status.03[mut_like_ixx.tcga]<-'TP53_MUTlike'
dat.tcga$TP53.final<-mut_like_status.03  
fit.km.03<-survfit(dat.tcga.obj ~TP53.final, data=dat.tcga)
ggsurvplot(fit.km.03, data = dat.tcga, pval = TRUE)

##############################################################################
##############################################################################




################################################################################
################################################################################
#percentages of each antigens in TP53MUT and WT
BEAT_ct_RNA_01.WES.allmut.IP<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_AML/BEAT_AML_preprocessed/BEAT_ct_RNA_01.WES.allmut.IP.rds')

TP53MUTWT_TP<-BEAT_ct_RNA_01.WES.allmut.IP[BEAT_ct_RNA_01.WES.allmut.IP$TP53_WES_class!='unknown',c(1,161,199:274)]
# TP53_IP<-BEAT_ct_RNA_01.WES.allmut.IP[BEAT_ct_RNA_01.WES.allmut.IP$TP53_WES_class=='TP53_MUT',c(1,161,199:274)]
# TP53WT_IP<-BEAT_ct_RNA_01.WES.allmut.IP[BEAT_ct_RNA_01.WES.allmut.IP$TP53_WES_class=='WT',c(1,161,199:274)]

library(reshape2)
shp01<-melt(TP53MUTWT_TP,measure.vars = colnames(IP.df))
shp01.sum<-shp01 %>% group_by(TP53_WES_class,variable,value) %>% summarize(n=n(),na.rm=TRUE) %>% mutate(freq=n/sum(n,na.rm=T))

shp01.sum$export <- with(shp01.sum, sprintf("%i (%.1f%%)", n, freq*100))
#reshape again
output <- dcast(variable+value~TP53_WES_class, value.var="export", data=shp01.sum, fill="missing") #use drop=F to prevent silent missings 
#'silent missings'
output$variable <- as.character(output$variable)
#make 'empty lines' 
empties <- data.frame(variable=unique(output$variable), stringsAsFactors=F)
empties[,colnames(output)[-1]] <- ""

#bind them together
output2 <- rbind(empties,output)
output2 <- output2[order(output2$variable,output2$value),]

#optional: 'remove' variable if value present
output2$variable[output2$value!=""] <- ""
write.csv(output2,"C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/Antigen analysis/BEAT_TP53MUT_WT_IP.csv")

#geom plot for ridge_score vs. survival 
TCGA_clinical_final<-readRDS('class ids/TCGA_clinical_final_class.rds')
BEAT_clinical_final<-readRDS('class ids/BEAT_clinical_final_class.rds')
BEAT_clinical_final.sub<-BEAT_clinical_final[BEAT_clinical_final$cls.GLM!='TP53_unk_WT_like',]
BEAT_clinical_final.sub$cls.GLM[BEAT_clinical_final.sub$cls.GLM=='TP53_WTunk_MT_like']<-'TP53MUT_like'
TCGA_clinical_final$cls.GLM[TCGA_clinical_final$cls.GLM=='TP53_unk_MT_like']<-'TP53MUT_like'

BEAT_clinical_final.sub$ridge.s1
BEAT_clinical_final.sub %>% ggplot(aes(x=ridge.s1*100,y=(overallSurvival/365),color=TP53_mut_stat)) +geom_point(size=3,alpha=0.65)+
  scale_color_manual(values=c("red", "grey"))
BEAT_clinical_final.sub %>% ggplot(aes(x=ridge.s1*100,y=(overallSurvival/365))) +geom_point(aes(shape=cls.GLM,color=cls.GLM),size=3.5,alpha=0.6)+
  scale_color_manual(values=c("red","grey",'blue'))+scale_shape_manual(values=c(16,16,18))

TCGA_clinical_final %>% ggplot(aes(x=ridge.s1*100,y=(OS_d/365),color=cls.GLM)) +geom_point(aes(shape=cls.GLM,color=cls.GLM),size=3,alpha=0.6)+
  scale_color_manual(values=c("red","grey",'blue'))+scale_shape_manual(values=c(16,16,18))
##
fig01<-BEAT_clinical_final.sub %>% ggplot(aes(x=(overallSurvival/365),y=ridge.s1*100,color=TP53_mut_stat)) +geom_point(size=4,alpha=0.65)+
  scale_color_manual(values=c("red", "grey"))
fig02<-BEAT_clinical_final.sub %>% ggplot(aes(x=(overallSurvival/365),y=ridge.s1*100)) +geom_point(aes(shape=cls.GLM,color=cls.GLM),size=4.5,alpha=0.6)+
  scale_color_manual(values=c("red","grey",'blue'))+scale_shape_manual(values=c(16,16,18))

fig03<-TCGA_clinical_final %>% ggplot(aes(x=(OS_d/365),y=ridge.s1*100,color=cls.GLM)) +geom_point(aes(shape=cls.GLM,color=cls.GLM),size=4,alpha=0.6)+
  scale_color_manual(values=c("red","grey",'blue'))+scale_shape_manual(values=c(16,16,18))

jpeg('BEAT_ridge_score_TP53mut.jpeg',width=576,height=480)
fig01
dev.off()

jpeg('BEAT_ridge_score_TP53mutlike.jpeg',width=576,height=480)
fig02
dev.off()

jpeg('TCGA_ridge_score_TP53mutlike.jpeg',width=576,height=480)
fig03
dev.off()

load('GLM_CV_GEP005.Rdata')
plot(TCGA_clinical.final$ridge.s1,TCGA_clinical.final$OS_d)
TCGA_clinical.final$TP53_mut_stat
TCGA_clinical.final%>% ggplot(aes(x=ridge.s1*100,y=(OS_d/365),color=TP53_mut_stat)) +geom_point(size=3,alpha=0.6)+
  scale_color_manual(values=c("red","grey"))#+scale_shape_manual(values=c(16,16,18))
fig04<-TCGA_clinical.final%>% ggplot(aes(x=(OS_d/365),y=ridge.s1*100,color=TP53_mut_stat)) +geom_point(size=4.5,alpha=0.6)+
  scale_color_manual(values=c("red","grey"))#+scale_shape_manual(values=c(16,16,18))

jpeg('TCGA_ridge_score_TP53mut.jpeg',width=576,height=480)
fig04
dev.off()
















################################################################################
################################################################################
#percentages of each antigens in TP53MUT and WT
BEAT_ct_RNA_01.WES.allmut.IP<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_AML/BEAT_AML_preprocessed/BEAT_ct_RNA_01.WES.allmut.IP.rds')

TP53MUTWT_TP<-BEAT_ct_RNA_01.WES.allmut.IP[BEAT_ct_RNA_01.WES.allmut.IP$TP53_WES_class!='unknown',c(1,161,199:274)]
# TP53_IP<-BEAT_ct_RNA_01.WES.allmut.IP[BEAT_ct_RNA_01.WES.allmut.IP$TP53_WES_class=='TP53_MUT',c(1,161,199:274)]
# TP53WT_IP<-BEAT_ct_RNA_01.WES.allmut.IP[BEAT_ct_RNA_01.WES.allmut.IP$TP53_WES_class=='WT',c(1,161,199:274)]

library(reshape2)
shp01<-melt(TP53MUTWT_TP,measure.vars = colnames(IP.df))
shp01.sum<-shp01 %>% group_by(TP53_WES_class,variable,value) %>% summarize(n=n(),na.rm=TRUE) %>% mutate(freq=n/sum(n,na.rm=T))

shp01.sum$export <- with(shp01.sum, sprintf("%i (%.1f%%)", n, freq*100))
#reshape again
output <- dcast(variable+value~TP53_WES_class, value.var="export", data=shp01.sum, fill="missing") #use drop=F to prevent silent missings 
#'silent missings'
output$variable <- as.character(output$variable)
#make 'empty lines' 
empties <- data.frame(variable=unique(output$variable), stringsAsFactors=F)
empties[,colnames(output)[-1]] <- ""

#bind them together
output2 <- rbind(empties,output)
output2 <- output2[order(output2$variable,output2$value),]

#optional: 'remove' variable if value present
output2$variable[output2$value!=""] <- ""
write.csv(output2,"C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/Antigen analysis/BEAT_TP53MUT_WT_IP.csv")

#geom plot for ridge_score vs. survival 
BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')
BEAT_clinical_final$cls.GLM
TCGA_clinical_final$cls.GLM.y[TCGA_clinical_final$cls.GLM.y=="TP53MUT_like"]<-"TP53_MUTlike"

# BEAT_clinical_final.sub<-BEAT_clinical_final[BEAT_clinical_final$cls.GLM!='TP53_unk_WT_like',]
# BEAT_clinical_final.sub$cls.GLM[BEAT_clinical_final.sub$cls.GLM=='TP53_WTunk_MT_like']<-'TP53MUT_like'
# TCGA_clinical_final$cls.GLM[TCGA_clinical_final$cls.GLM=='TP53_unk_MT_like']<-'TP53MUT_like'
BEAT_clinical_final$TP53_mut_stat<-ifelse(BEAT_clinical_final$cls.GLM=='TP53_MUT','TP53_MUT','TP53_WT')
BEAT_clinical_final$ridge.s1
BEAT_clinical_final %>% ggplot(aes(x=ridge.s1*100,y=(overallSurvival/365),color=TP53_mut_stat)) +geom_point(size=3,alpha=0.65)+
  scale_color_manual(values=c("red", "grey"))
BEAT_clinical_final %>% ggplot(aes(x=ridge.s1*100,y=(overallSurvival/365))) +geom_point(aes(shape=cls.GLM,color=cls.GLM),size=3.5,alpha=0.6)+
  scale_color_manual(values=c("red",'blue',"grey"))+scale_shape_manual(values=c(16,16,18))

TCGA_clinical_final %>% ggplot(aes(x=ridge.s1*100,y=(OS_d/365),color=cls.GLM.y)) +geom_point(aes(shape=cls.GLM.y,color=cls.GLM.y),size=3,alpha=0.6)+
  scale_color_manual(values=c("red","grey",'blue'))+scale_shape_manual(values=c(16,16,18))
##
fig01<-BEAT_clinical_final %>% ggplot(aes(y=ridge.s1*100,x=(overallSurvival/365),color=TP53_mut_stat)) +geom_point(size=3,alpha=0.65)+
  scale_color_manual(values=c("red", "grey"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                    legend.key=element_rect(fill="white"),axis.title = element_text(size=14),title = element_text(size=14),
                                                    legend.text=element_text(size=14))
fig02<-BEAT_clinical_final %>% ggplot(aes(y=ridge.s1*100,x=(overallSurvival/365))) +geom_point(aes(shape=cls.GLM,color=cls.GLM),size=3.5,alpha=0.6)+
  scale_color_manual(values=c("red",'blue',"grey"))+scale_shape_manual(values=c(16,18,16))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                 legend.key=element_rect(fill="white"),axis.title = element_text(size=14),title = element_text(size=14),
                                                                                                 legend.text=element_text(size=14))


fig03<-TCGA_clinical_final %>% ggplot(aes(y=ridge.s1*100,x=(OS_d/365),color=cls.GLM.y)) +geom_point(aes(shape=cls.GLM.y,color=cls.GLM.y),size=3,alpha=0.6)+
  scale_color_manual(values=c("red",'blue',"grey"))+scale_shape_manual(values=c(16,18,16))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                 legend.key=element_rect(fill="white"),axis.title = element_text(size=14),title = element_text(size=14),
                                                                                                 legend.text=element_text(size=14))

setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/")
jpeg('BEAT_ridge_score_TP53mut.jpeg',width=5.2,height=3.8,units='in',res=300)
fig01
dev.off()

jpeg('BEAT_ridge_score_TP53mutlike.jpeg',width=5.2,height=3.8,units='in',res=300)
fig02
dev.off()

jpeg('TCGA_ridge_score_TP53mutlike.jpeg',width=5.2,height=3.8,units='in',res=300)
fig03
dev.off()

load('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/GLM_CV_GEP005.Rdata')
plot(TCGA_clinical.final$ridge.s1,TCGA_clinical.final$OS_d)
TCGA_clinical.final$TP53_mut_stat
TCGA_clinical.final%>% ggplot(aes(x=ridge.s1*100,y=(OS_d/365),color=TP53_mut_stat)) +geom_point(size=3,alpha=0.6)+
  scale_color_manual(values=c("red","grey"))#+scale_shape_manual(values=c(16,16,18))
fig04<-TCGA_clinical.final%>% ggplot(aes(x=(OS_d/365),y=ridge.s1*100,color=TP53_mut_stat)) +geom_point(size=3,alpha=0.6)+
  scale_color_manual(values=c("red","grey"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                   legend.key=element_rect(fill="white"),axis.title = element_text(size=14),title = element_text(size=14),
                                                   legend.text=element_text(size=14))#+scale_shape_manual(values=c(16,16,18))

jpeg('TCGA_ridge_score_TP53mut.jpeg',width=5.2,height=3.8,units='in',res=300)
fig04
dev.off()