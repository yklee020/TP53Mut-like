

#2022.07.25
#Drug Sensitivity 
library(readxl)
library(dplyr)
#library(reshape2)
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version")

BEAT_cyto_clean<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/BEAT_cyto_clean.rds')
BEAT_wt_id<-BEAT_cyto_clean$LabId
BEAT_Drug_AUC<-read_excel("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/BEAT_AML_drug_sensitivity_AUC.xlsx")
#high auc resist and low auc sensitive; so transform the score to high=sensitive
BEAT_Drug_AUC.scale<-BEAT_Drug_AUC %>% group_by(inhibitor)%>%mutate(auc_scale=scale(auc)*(-1))
#test.df<-scale(BEAT_Drug_AUC.scale$auc[BEAT_Drug_AUC.scale$inhibitor=='Venetoclax'])*(-1)


patiend_id<-unique(BEAT_Drug_AUC.scale$lab_id)
BEAT_Drug_AUC.rev<-BEAT_Drug_AUC.scale %>% filter(lab_id %in% BEAT_wt_id)
BEAT_Drug_AUC.rev<-BEAT_Drug_AUC.rev %>% group_by(inhibitor)%>%mutate(auc_scale=scale(auc)*(-1))
##intersect(BEAT_wt_id, BEAT_Drug_AUC.rev$lab_id)
BEAT_Drug_AUC.rev$cls.GLM<-rep(0,dim(BEAT_Drug_AUC.rev)[1])
BEAT_Drug_AUC.rev$cls.GLM[which(BEAT_Drug_AUC.rev$lab_id %in% BEAT_cyto_clean$LabId[BEAT_cyto_clean$cls.GLM.x=='TP53_MUT'])]<-'TP53_MUT'
BEAT_Drug_AUC.rev$cls.GLM[which(BEAT_Drug_AUC.rev$lab_id %in% BEAT_cyto_clean$LabId[BEAT_cyto_clean$cls.GLM.x=='TP53_WT'])]<-'TP53_WT'
BEAT_Drug_AUC.rev$cls.GLM[which(BEAT_Drug_AUC.rev$lab_id %in% BEAT_cyto_clean$LabId[BEAT_cyto_clean$cls.GLM.x=='TP53MUT_like'])]<-'TP53MUT_like'

#!!!!watch out for dplyr and reshape2 masking group_by function it may cause scaling function to scale entire data!!!

#z-score version
#drug_res<-BEAT_Drug_AUC.rev %>% group_by(cls.GLM,inhibitor) %>% summarise(mean=mean(auc_scale),sd=sd(auc_scale))
drug_res<-BEAT_Drug_AUC.rev %>% group_by(cls.GLM,inhibitor) %>% summarise_at(vars(auc_scale),list(mean=mean,sd=sd))
drug_res.c<-BEAT_Drug_AUC.rev %>% group_by(cls.GLM,inhibitor) %>% summarise(ct=n())


#reshape again
library(reshape2)
drug_output <- dcast(cls.GLM~inhibitor, value.var="mean", data=drug_res) #use drop=F to prevent silent missings 
#auc version
drug_res.auc<-BEAT_Drug_AUC.rev %>% group_by(cls.GLM,inhibitor) %>% summarise_at(vars(auc),list(mean=mean,sd=sd))
#reshape again
drug_output.auc <- dcast(cls.GLM~inhibitor, value.var="mean", data=drug_res.auc) #use drop=F to prevent silent missings 


# BEAT_Drug_AUC.sub<-BEAT_Drug_AUC.rev[c(which(BEAT_Drug_AUC.rev$inhibitor=='Venetoclax'),which(BEAT_Drug_AUC.rev$inhibitor=='Staurosporine')),]
# BEAT_Drug_AUC.sub<-BEAT_Drug_AUC.sub %>% group_by(inhibitor) %>% mutate(ss=(scale(auc))*(-1))
# head(scale(BEAT_Drug_AUC.sub$auc[BEAT_Drug_AUC.sub$inhibitor=='Venetoclax'])*-1)
# #auc
# drug_auc<-BEAT_Drug_AUC.rev %>% group_by(cls.GLM,inhibitor) %>% summarise_at(vars(auc),list(mean=mean,sd=sd))
# drug_auc.out <- dcast(cls.GLM~inhibitor, value.var="mean", data=drug_auc) 


library(ComplexHeatmap)
a<-drug_output[,2:123]
rownames(a)<-drug_output[,1]
Heatmap(as.matrix(a), name = "TP53AUC.z", clustering_distance_columns = "pearson",
        clustering_method_columns = "average",#complete,ward.D2
        #column_title = "pre-defined distance method (1 - pearson),complete",
        show_row_names = T,show_column_names = T,column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 8),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 12),
        show_row_dend = F)
colnames(a)<-1:dim(a)[2]
inh_anno<-rowAnnotation(foo=anno_mark(at=c(28,102,117,32,76),
                                      labels=c('Elesclomol','Staurosporine','Venetoclax','Flavopiridol','Panobinostat')),
                        annotation_name_gp= gpar(fontsize = 24))
p.out<-Heatmap(t(as.matrix(a)), name = "TP53AUC.z", clustering_distance_rows = "pearson",
               clustering_method_rows = "average",#complete,ward.D2
               #column_title = "pre-defined distance method (1 - pearson),complete",
               show_row_names = F,show_column_names = T,column_names_max_height = unit(16, "cm"),
               column_names_gp = gpar(fontsize = 16),
               row_names_max_width = unit(12, "cm"),
               row_names_gp = gpar(fontsize = 4),
               show_row_dend = F,column_names_rot = 0,
               column_order = c('TP53_MUT','TP53MUT_like','TP53_WT'),right_annotation = inh_anno)
#jpeg('Analysis/Drug Sensitivity/BEA_AUC_Z.jpeg',width=5.2,height=3.8,units='in',res=300)
p.out
#dev.off()
#jpeg('Analysis/Drug Sensitivity/BEA_AUC_Z_legend.jpeg',width=640,height=896,res=300)
p.out
#dev.off()

b<-drug_output.auc[,2:123]
rownames(b)<-drug_output.auc[,1]
Heatmap((as.matrix(b)), name = "TP53AUC", #clustering_distance_columns = "pearson",
        cluster_columns = T,
        #clustering_method_columns = "average",#complete,ward.D2
        #column_title = "pre-defined distance method (1 - pearson),complete",
        show_row_names = T,show_column_names = T,column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 5),
        show_row_dend = F)

test.group<-as.symbol('cls.GLM')
num_cat<-colnames(drug_output.auc)[c(2:123)]
p.out.df<-as.data.frame(matrix(0,6,length(num_cat)))
colnames(p.out.df)<-num_cat

for(i in 1:dim(p.out.df)[2]){
  test_dat<-as.data.frame(BEAT_Drug_AUC.rev)%>%filter(inhibitor == num_cat[i])
  p<-aov(auc_scale~eval(test.group), data=as.data.frame(test_dat)) %>% TukeyHSD()
  mean_mt<-mean(test_dat$auc_scale[test_dat$cls.GLM=='TP53_MUT'])
  mean_wt<-mean(test_dat$auc_scale[test_dat$cls.GLM=='TP53_WT'])
  mean_mtl<-mean(test_dat$auc_scale[test_dat$cls.GLM=='TP53MUT_like'])
  mean_diff_mtwt<-(mean_mt-mean_wt)
  mean_diff_mtlmt<-(mean_mtl-mean_mt)
  mean_diff_mtlwt<-(mean_mtl-mean_wt)
  p.out.df[c(1,3,5),i]<-c(mean_diff_mtwt,mean_diff_mtlmt,mean_diff_mtlwt)
  p.out.df[c(2,4,6),i]<-p$`eval(test.group)`[,4]
}
rownames(p.out.df)[c(2,4,6)]<-rownames(p$`eval(test.group)`)
rownames(p.out.df)[c(1,3,5)]<-c('mean_diff_mtwt','mean_diff_mtlmt','mean_diff_mtlwt')

# test_dat<-as.data.frame(BEAT_Drug_AUC.rev)%>%filter(inhibitor == num_cat[117])
# p<-aov(auc_scale~eval(test.group), data=as.data.frame(test_dat)) %>% TukeyHSD()

p.out.df.auc<-as.data.frame(matrix(0,6,length(num_cat)))
colnames(p.out.df.auc)<-num_cat

for(i in 1:dim(p.out.df.auc)[2]){
  test_dat<-as.data.frame(BEAT_Drug_AUC.rev)%>%filter(inhibitor == num_cat[i])
  p<-aov(auc~eval(test.group), data=as.data.frame(test_dat)) %>% TukeyHSD()
  mean_mt<-mean(test_dat$auc[test_dat$cls.GLM=='TP53_MUT'])
  mean_wt<-mean(test_dat$auc[test_dat$cls.GLM=='TP53_WT'])
  mean_mtl<-mean(test_dat$auc[test_dat$cls.GLM=='TP53MUT_like'])
  mean_diff_mtwt<-(mean_mt-mean_wt)
  mean_diff_mtlmt<-(mean_mtl-mean_mt)
  mean_diff_mtlwt<-(mean_mtl-mean_wt)
  p.out.df.auc[c(1,3,5),i]<-c(mean_diff_mtwt,mean_diff_mtlmt,mean_diff_mtlwt)
  p.out.df.auc[c(2,4,6),i]<-p$`eval(test.group)`[,4]
}
rownames(p.out.df.auc)[c(2,4,6)]<-rownames(p$`eval(test.group)`)
rownames(p.out.df.auc)[c(1,3,5)]<-c('mean_diff_mtwt','mean_diff_mtlmt','mean_diff_mtlwt')

#do t.test version
p.out.df.ttest<-as.data.frame(matrix(0,12,length(num_cat)))
colnames(p.out.df.ttest)<-num_cat

for(i in 1:dim(p.out.df.ttest)[2]){
  test_dat<-as.data.frame(BEAT_Drug_AUC.rev)%>%filter(inhibitor == num_cat[i])
  #p<-aov(auc~eval(test.group), data=as.data.frame(test_dat)) %>% TukeyHSD()
  p<-pairwise.t.test(test_dat$auc_scale,test_dat$cls.GLM, p.adj = "none",na.rm=T,
                  paired=FALSE,pool.sd = F) 
  mean_mt.auc.scale<-mean(test_dat$auc_scale[test_dat$cls.GLM=='TP53_MUT'])
  mean_wt.auc.scale<-mean(test_dat$auc_scale[test_dat$cls.GLM=='TP53_WT'])
  mean_mtl.auc.scale<-mean(test_dat$auc_scale[test_dat$cls.GLM=='TP53MUT_like'])
  
  mean_mt.auc<-mean(test_dat$auc[test_dat$cls.GLM=='TP53_MUT'])
  mean_wt.auc<-mean(test_dat$auc[test_dat$cls.GLM=='TP53_WT'])
  mean_mtl.auc<-mean(test_dat$auc[test_dat$cls.GLM=='TP53MUT_like'])
  
  mean_diff_mtwt.scale<-(mean_mt.auc.scale-mean_wt.auc.scale)
  mean_diff_mtlmt.scale<-(mean_mtl.auc.scale-mean_mt.auc.scale)
  mean_diff_mtlwt.scale<-(mean_mtl.auc.scale-mean_wt.auc.scale)
  mean_diff_mtwt<-(mean_mt.auc-mean_wt.auc)
  mean_diff_mtlmt<-(mean_mtl.auc-mean_mt.auc)
  mean_diff_mtlwt<-(mean_mtl.auc-mean_wt.auc)
  p.out.df.ttest[c(1,5,9),i]<-c(mean_diff_mtwt.scale,mean_diff_mtlmt.scale,mean_diff_mtlwt.scale)
  p.out.df.ttest[c(2,6,10),i]<-c(mean_diff_mtwt,mean_diff_mtlmt,mean_diff_mtlwt)
  p.out.df.ttest[c(3,7,11),i]<-c(p$p.value[1,1],p$p.value[2,1],p$p.value[2,2])
  
}
rownames(p.out.df.ttest)[c(1,5,9)]<-c('mean_diff_mtwt.scale','mean_diff_mtlmt.scale','mean_diff_mtlwt.scale')
rownames(p.out.df.ttest)[c(2,6,10)]<-c('mean_diff_mtwt','mean_diff_mtlmt','mean_diff_mtlwt')
rownames(p.out.df.ttest)[c(3,7,11)]<-c('mtwt.p.value','mtlmt.p.value','mtlwt.p.value')
rownames(p.out.df.ttest)[c(4,8,12)]<-c('mtwt.adjp.value','mtlmt.adjp.value','mtlwt.adjp.value')
p.out.df.ttest[4,]<-p.adjust( p.out.df.ttest[3,],'BH')
p.out.df.ttest[8,]<-p.adjust( p.out.df.ttest[7,],'BH')
p.out.df.ttest[12,]<-p.adjust( p.out.df.ttest[11,],'BH')



p.out.df.t<-as.data.frame(t(p.out.df))
p.out.df.t$FDR_mtwt<-p.out.df.t$`TP53_WT-TP53_MUT`
p.out.df.t$FDR_mtlmt<-p.out.df.t$`TP53MUT_like-TP53_MUT`
p.out.df.t$FDR_mtlwt<-p.out.df.t$`TP53MUT_like-TP53_WT`

z.mtwt<-p.out.df.t[,c('mean_diff_mtwt','FDR_mtwt')]
z.mtlmt<-p.out.df.t[,c('mean_diff_mtlmt','FDR_mtlmt')]
z.mtlwt<-p.out.df.t[,c('mean_diff_mtlwt','FDR_mtlwt')]

library(ggrepel)
ggplot(data=z.mtwt, aes(x=mean_diff_mtwt, y=-log10(FDR_mtwt),
                        label=rownames(z.mtwt))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

ggplot(data=z.mtlwt, aes(x=mean_diff_mtlwt, y=-log10(FDR_mtlwt),
                         label=rownames(z.mtlwt))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")#+
  #geom_hline(yintercept=-log10(0.1), col="blue")

ggplot(data=z.mtlmt, aes(x=mean_diff_mtlmt, y=-log10(FDR_mtlmt),
                         label=rownames(z.mtlmt))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")#+
#geom_hline(yintercept=-log10(0.1), col="blue")

p.out.df.auc.t<-as.data.frame(t(p.out.df.auc))
p.out.df.auc.t$FDR_mtwt<-p.adjust(p.out.df.auc.t$`TP53_WT-TP53_MUT`,method="BH")
p.out.df.auc.t$FDR_mtlmt<-p.adjust(p.out.df.auc.t$`TP53MUT_like-TP53_MUT`,method="BH")
p.out.df.auc.t$FDR_mtlwt<-p.adjust(p.out.df.auc.t$`TP53MUT_like-TP53_WT`,method="BH")

auc.mtwt<-p.out.df.auc.t[,c('mean_diff_mtwt','FDR_mtwt')]
auc.mtlmt<-p.out.df.auc.t[,c('mean_diff_mtlmt','FDR_mtlmt')]
auc.mtlwt<-p.out.df.auc.t[,c('mean_diff_mtlwt','FDR_mtlwt')]


ggplot(data=auc.mtwt, aes(x=mean_diff_mtwt, y=-log10(FDR_mtwt),
                          label=rownames(auc.mtwt))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")#+
  #geom_hline(yintercept=-log10(0.1), col="blue")

ggplot(data=auc.mtlwt, aes(x=mean_diff_mtlwt, y=-log10(FDR_mtlwt),
                           label=rownames(auc.mtlwt))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")#+
  #geom_hline(yintercept=-log10(0.1), col="blue")


#just to do with p value from anova
auc.mtwt.p<-p.out.df.auc.t[,c('mean_diff_mtwt','TP53_WT-TP53_MUT')]
auc.mtlmt.p<-p.out.df.auc.t[,c('mean_diff_mtlmt','TP53MUT_like-TP53_MUT')]
auc.mtlwt.p<-p.out.df.auc.t[,c('mean_diff_mtlwt','TP53MUT_like-TP53_WT')]

ggplot(data=auc.mtwt.p, aes(x=mean_diff_mtwt, y=-log10(`TP53_WT-TP53_MUT`),
                            label=rownames(auc.mtwt.p))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")#+
#  geom_hline(yintercept=-log10(0.1), col="blue")

ggplot(data=auc.mtlwt.p, aes(x=mean_diff_mtlwt, y=-log10(`TP53MUT_like-TP53_WT`),
                             label=rownames(auc.mtlwt.p))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")#+
  #geom_hline(yintercept=-log10(0.1), col="blue")

#do the t.test version

p.out.df.ttest.t<-as.data.frame(t(p.out.df.ttest))
p.out.df.ttest.t$FDR_mtwt<-p.out.df.ttest.t$mtwt.adjp.value
p.out.df.ttest.t$FDR_mtlmt<-p.out.df.ttest.t$mtlmt.adjp.value
p.out.df.ttest.t$FDR_mtlwt<-p.out.df.ttest.t$mtlwt.adjp.value

z.mtwt.ttest<-p.out.df.ttest.t[,c('mean_diff_mtwt.scale','mean_diff_mtwt','FDR_mtwt')]
z.mtwt.ttest$drug<-rownames(z.mtwt.ttest)
z.mtlmt.ttest<-p.out.df.ttest.t[,c('mean_diff_mtlmt.scale','mean_diff_mtlmt','FDR_mtlmt')]
z.mtlmt.ttest$drug<-rownames(z.mtlmt.ttest)
z.mtlwt.ttest<-p.out.df.ttest.t[,c('mean_diff_mtlwt.scale','mean_diff_mtlwt','FDR_mtlwt')]
z.mtlwt.ttest$drug<-rownames(z.mtlwt.ttest)

#subset labels 
z.mtwt.ttest$drug<-rownames(z.mtwt.ttest)
z.mtlwt.ttest$drug<-rownames(z.mtlwt.ttest)
z.mtlmt.ttest$drug<-rownames(z.mtlmt.ttest)

z.mtwt.ttest$drug<-ifelse(z.mtwt.ttest$mean_diff_mtwt>-20 &z.mtwt.ttest$mean_diff_mtwt<40,"",z.mtwt.ttest$drug)
z.mtlwt.ttest$drug<-ifelse(z.mtlwt.ttest$mean_diff_mtlwt>-20 &z.mtlwt.ttest$mean_diff_mtlwt<20,"",z.mtlwt.ttest$drug)
z.mtlmt.ttest$drug<-ifelse(z.mtlmt.ttest$mean_diff_mtlmt>-20 &z.mtlmt.ttest$mean_diff_mtlmt<20,"",z.mtlmt.ttest$drug)

library(ggrepel)
ggplot(data=z.mtwt.ttest, aes(x=mean_diff_mtwt.scale, y=-log10(FDR_mtwt),
                        label=rownames(z.mtwt.ttest))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black")) 
ggplot(data=z.mtlwt.ttest, aes(x=mean_diff_mtlwt.scale, y=-log10(FDR_mtlwt),
                         label=rownames(z.mtlwt.ttest))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
#geom_hline(yintercept=-log10(0.1), col="blue")
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggplot(data=z.mtlmt.ttest, aes(x=mean_diff_mtlmt.scale, y=-log10(FDR_mtlmt),
                               label=rownames(z.mtlmt.ttest))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  #geom_hline(yintercept=-log10(0.1), col="blue")
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#ttest auc version
beat.auc.mtwt<-ggplot(data=z.mtwt.ttest, aes(x=mean_diff_mtwt, y=-log10(FDR_mtwt),
                              label=drug)) +
  geom_point(size=0.75) + 
  theme_minimal() +
  geom_text_repel(min.segment.length = Inf, seed = 42, box.padding = 0.3,size =3.5) +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=10)) 

#beat.auc.mtlwt<-
beat.auc.mtlwt<-ggplot(data=z.mtlwt.ttest, aes(x=mean_diff_mtlwt, y=-log10(FDR_mtlwt),
                               label=drug)) +
  geom_point(size=0.75) + 
  theme_minimal() +
  geom_text_repel(min.segment.length = Inf, seed = 42, box.padding = 0.2,size=3.5) +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  #geom_hline(yintercept=-log10(0.1), col="blue")
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=10))
beat.auc.mtlmt<-ggplot(data=z.mtlmt.ttest, aes(x=mean_diff_mtlmt, y=-log10(FDR_mtlmt),
                                               label=drug)) +
  geom_point(size=0.75) + 
  theme_minimal() +
  geom_text_repel(min.segment.length = Inf, seed = 42, box.padding = 0.2,size=3.5) +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")+
  #geom_hline(yintercept=-log10(0.1), col="blue")
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=10))

jpeg('Analysis/Drug Sensitivity/BEAT_auc_ttest_mtwt.jpeg',width=5.2,height=3.8,units='in',res=300)
beat.auc.mtwt
dev.off()
jpeg('Analysis/Drug Sensitivity/BEAT_auc_ttest_mtlwt.jpeg',width=5.2,height=3.8,units='in',res=300)
beat.auc.mtlwt
dev.off()
jpeg('Analysis/Drug Sensitivity/BEAT_auc_ttest_mtlmt.jpeg',width=5.2,height=3.8,units='in',res=300)
beat.auc.mtlmt
dev.off()
#generate venetoclax plot
library(rstatix)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = (sd(x[[col]], na.rm=TRUE))/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
venetoclax.df<-BEAT_Drug_AUC.rev[BEAT_Drug_AUC.rev$inhibitor=='Venetoclax',]

veneto.sum<-data_summary(venetoclax.df,varname='auc_scale',groupnames=c('cls.GLM'))

library(ggpubr)
veneto.p<-ggdotplot(venetoclax.df,x='cls.GLM',y='auc',add = 'mean_se',color='black',fill='black',size=4.5,
                    position=position_jitter(0.3),binwidth = 1.5,alpha=0.5,order=c('TP53_MUT','TP53MUT_like','TP53_WT'))+
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  
  stat_compare_means(method='t.test',p.adjust.method='BH',hide.ns=T,label ='none', #'p.signif',
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53MUT_like','TP53_WT')))+
  font("xlab", size = 12)+font("ylab", size = 14)+theme(axis.text=element_text(size=14))

# veneto.p<-ggdotplot(venetoclax.df,x='cls.GLM',y='auc',add = 'mean_se',color='black',fill='black',size=4.5,
#                     position=position_jitter(0.3),binwidth = 1.5,alpha=0.5,order=c('TP53_MUT','TP53MUT_like','TP53_WT'))+
#   stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
#   stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
#   
#   stat_compare_means(method='t.test',p.adjust.method='BH',hide.ns=T,label ='none', #'p.signif',
#                      comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53MUT_like','TP53_WT'),c('TP53MUT_like','TP53_MUT')))+
#   font("xlab", size = 12)+font("ylab", size = 14)+theme(axis.text=element_text(size=14))


stat.veneto <- venetoclax.df %>%
  #group_by(cls.GLM) %>%
  wilcox_test(auc_scale~cls.GLM) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
p.adjust(stat.veneto[c(1,3),]$p,'BH')

# jpeg('auc_veneto_scale.jpeg',width=1024,height=640)
# veneto.p.scale
# dev.off()

jpeg('Analysis/Drug Sensitivity/Veneto_auc.jpeg',width=3.3,height=3.5,units='in',res=300)
veneto.p
dev.off()

write.csv(p.out.df.ttest.t,'Analysis/Drug Sensitivity/Drug_sensitivity_final_results.csv')




###################################################################################
#replicate original results with full ids
#get TP53Mut cases
mutation_calls_wes<-read_excel('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/TableS7_variant_analysis.xlsx')
mutation_calls_wes.tp53<-mutation_calls_wes[mutation_calls_wes$symbol=='TP53',]
wes_tp53_labid<-unique(mutation_calls_wes.tp53$labId)#54
#discard that does not have exomeseq data
BEAT_count_clinical_anno_full.01.excel<-read_excel('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/BEAT AML full clinical summary.xlsx')
exomeid<-BEAT_count_clinical_anno_full.01.excel$LabId[BEAT_count_clinical_anno_full.01.excel$exomeSeq=='y']

BEAT_Drug_AUC.ori<-BEAT_Drug_AUC.scale 
BEAT_Drug_AUC.ori.filt<-BEAT_Drug_AUC.ori %>% filter(lab_id %in% exomeid)
#intersect(unique(BEAT_Drug_AUC.ori.filt$lab_id),exomeid)
BEAT_Drug_AUC.ori.filt$TP53_mut<-rep(0,dim(BEAT_Drug_AUC.ori.filt)[1])
BEAT_Drug_AUC.ori.filt$TP53_mut<-ifelse(BEAT_Drug_AUC.ori.filt$lab_id %in%wes_tp53_labid,'TP53_MUT','TP53_WT')

# unique(BEAT_Drug_AUC.ori.filt$lab_id[BEAT_Drug_AUC.ori.filt$TP53_mut=='TP53_MUT'])#MT:31
# unique(BEAT_Drug_AUC.ori.filt$lab_id[BEAT_Drug_AUC.ori.filt$TP53_mut=='TP53_WT'])#WT:339
# BEAT_Drug_AUC.ori.filt

#z-score version
drug_res.ori<-BEAT_Drug_AUC.ori.filt %>% group_by(TP53_mut,inhibitor) %>% summarise(mean=mean(auc_scale),sd=sd(auc_scale))
#reshape again
drug_output.ori <- dcast(TP53_mut~inhibitor, value.var="mean", data=drug_res.ori) #use drop=F to prevent silent missings 
#auc version
drug_res.ori.auc<-BEAT_Drug_AUC.ori.filt %>% group_by(TP53_mut,inhibitor) %>% summarise(mean=mean(auc),sd=sd(auc))
#reshape again
drug_output.ori.auc <- dcast(TP53_mut~inhibitor, value.var="mean", data=drug_res.ori.auc) #use drop=F to prevent silent missings 

library(ComplexHeatmap)
aa<-drug_output.ori[,2:123]
rownames(aa)<-drug_output.ori[,1]
Heatmap(as.matrix(aa), name = "TP53AUC.z", #clustering_distance_columns = "pearson",
        #clustering_method_columns = "average",#complete,ward.D2
        #column_title = "pre-defined distance method (1 - pearson),complete",
        show_row_names = T,show_column_names = T,column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 5),
        show_row_dend = F)
bb<-drug_output.ori.auc[,2:123]
rownames(bb)<-drug_output.ori.auc[,1]
Heatmap(as.matrix(bb), name = "TP53AUC", #clustering_distance_columns = "pearson",
        #clustering_method_columns = "average",#complete,ward.D2
        #column_title = "pre-defined distance method (1 - pearson),complete",
        show_row_names = T,show_column_names = T,column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 5),
        show_row_dend = F)


test.group<-as.symbol('TP53_mut')
num_cat.ori<-colnames(drug_output.ori)[c(2:123)]
p.out.ori.z.df<-as.data.frame(matrix(0,2,length(num_cat.ori)))
colnames(p.out.ori.z.df)<-num_cat.ori
rownames(p.out.ori.z.df)<-c('mean_diff_mtwt','p.value')

for(i in 1:dim(p.out.ori.z.df)[2]){
  test_dat<-as.data.frame(BEAT_Drug_AUC.ori.filt)%>%filter(inhibitor == num_cat.ori[i])
  p<-t.test(test_dat$auc_scale[test_dat$TP53_mut=='TP53_MUT'],test_dat$auc_scale[test_dat$TP53_mut=='TP53_WT'])
  mt_mean<-mean(test_dat$auc_scale[test_dat$TP53_mut=='TP53_MUT'])
  wt_mean<-mean(test_dat$auc_scale[test_dat$TP53_mut=='TP53_WT'])
  mean_diff<-(mt_mean-wt_mean)
  
  p.out.ori.z.df[,i]<-c(mean_diff,p$p.value)
}

p.out.ori.auc.df<-as.data.frame(matrix(0,2,length(num_cat.ori)))
colnames(p.out.ori.auc.df)<-num_cat.ori
rownames(p.out.ori.auc.df)<-c('mean_diff_mtwt','p.value')
for(i in 1:dim(p.out.ori.auc.df)[2]){
  test_dat<-as.data.frame(BEAT_Drug_AUC.ori.filt)%>%filter(inhibitor == num_cat.ori[i])
  p<-t.test(test_dat$auc[test_dat$TP53_mut=='TP53_MUT'],test_dat$auc[test_dat$TP53_mut=='TP53_WT'])
  mt_mean<-mean(test_dat$auc[test_dat$TP53_mut=='TP53_MUT'])
  wt_mean<-mean(test_dat$auc[test_dat$TP53_mut=='TP53_WT'])
  mean_diff<-(mt_mean-wt_mean)
  
  p.out.ori.auc.df[,i]<-c(mean_diff,p$p.value)
}

p.out.ori.z.df.t<-as.data.frame(t(p.out.ori.z.df))
p.out.ori.z.df.t$FDR<-p.adjust(p.out.ori.z.df.t$p.value,method="BH")

p.out.ori.auc.df.t<-as.data.frame(t(p.out.ori.auc.df))
p.out.ori.auc.df.t$FDR<-p.adjust(p.out.ori.auc.df.t$p.value,method="BH")


ggplot(data=p.out.ori.z.df.t, aes(x=mean_diff_mtwt, y=-log10(FDR),
                                  label=rownames(p.out.ori.z.df.t))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

ggplot(data=p.out.ori.auc.df.t, aes(x=mean_diff_mtwt, y=-log10(FDR),
                                    label=rownames(p.out.ori.auc.df.t))) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual(values=c("blue", "black", "red")) +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")






