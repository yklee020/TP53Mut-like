
#geom plot for ridge_score vs. survival 
BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')
BEAT_clinical_final$cls.GLM
TCGA_clinical_final$cls.GLM.y[TCGA_clinical_final$cls.GLM.y=="TP53MUT_like"]<-"TP53_MUTlike"

BEAT_clinical_final_ninit<-BEAT_clinical_final[BEAT_clinical_final$ELN2017=='NonInitial',]
BEAT_clinical_final_init<-BEAT_clinical_final[BEAT_clinical_final$ELN2017!='NonInitial',]


# BEAT_clinical_final.sub<-BEAT_clinical_final[BEAT_clinical_final$cls.GLM!='TP53_unk_WT_like',]
# BEAT_clinical_final.sub$cls.GLM[BEAT_clinical_final.sub$cls.GLM=='TP53_WTunk_MT_like']<-'TP53MUT_like'
# TCGA_clinical_final$cls.GLM[TCGA_clinical_final$cls.GLM=='TP53_unk_MT_like']<-'TP53MUT_like'
BEAT_clinical_final_init$TP53_mut_stat<-ifelse(BEAT_clinical_final_init$cls.GLM=='TP53_MUT','TP53_MUT','TP53_WT')
BEAT_clinical_final_ninit$TP53_mut_stat<-ifelse(BEAT_clinical_final_ninit$cls.GLM=='TP53_MUT','TP53_MUT','TP53_WT')

BEAT_clinical_final_init$ridge.s1
BEAT_clinical_final_ninit$ridge.s1

# BEAT_clinical_final %>% ggplot(aes(x=ridge.s1*100,y=(overallSurvival/365),color=TP53_mut_stat)) +geom_point(size=3,alpha=0.65)+
#   scale_color_manual(values=c("red", "grey"))
# BEAT_clinical_final %>% ggplot(aes(x=ridge.s1*100,y=(overallSurvival/365))) +geom_point(aes(shape=cls.GLM,color=cls.GLM),size=3.5,alpha=0.6)+
#   scale_color_manual(values=c("red",'blue',"grey"))+scale_shape_manual(values=c(16,16,18))
BEAT_clinical_final_init %>% ggplot(aes(x=ridge.s1*100,y=(overallSurvival/365),color=TP53_mut_stat)) +geom_point(size=3,alpha=0.65)+
  scale_color_manual(values=c("red", "grey"))
BEAT_clinical_final_init %>% ggplot(aes(x=ridge.s1*100,y=(overallSurvival/365))) +geom_point(aes(shape=cls.GLM,color=cls.GLM),size=3.5,alpha=0.6)+
  scale_color_manual(values=c("red",'blue',"grey"))+scale_shape_manual(values=c(16,16,18))

BEAT_clinical_final_ninit %>% ggplot(aes(x=ridge.s1*100,y=(overallSurvival/365),color=TP53_mut_stat)) +geom_point(size=3,alpha=0.65)+
  scale_color_manual(values=c("red", "grey"))
BEAT_clinical_final_ninit %>% ggplot(aes(x=ridge.s1*100,y=(overallSurvival/365))) +geom_point(aes(shape=cls.GLM,color=cls.GLM),size=3.5,alpha=0.6)+
  scale_color_manual(values=c("red",'blue',"grey"))+scale_shape_manual(values=c(16,16,18))




TCGA_clinical_final %>% ggplot(aes(x=ridge.s1*100,y=(OS_d/365),color=cls.GLM.y)) +geom_point(aes(shape=cls.GLM.y,color=cls.GLM.y),size=3,alpha=0.6)+
  scale_color_manual(values=c("red","grey",'blue'))+scale_shape_manual(values=c(16,16,18))
##
fig01.init<-BEAT_clinical_final_init %>% ggplot(aes(y=ridge.s1*100,x=(overallSurvival/365),color=TP53_mut_stat)) +geom_point(size=3,alpha=0.65)+
  scale_color_manual(values=c("red", "grey"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                    legend.key=element_rect(fill="white"),axis.title = element_text(size=14),title = element_text(size=14),
                                                    legend.text=element_text(size=14))
fig02.init<-BEAT_clinical_final_init %>% ggplot(aes(y=ridge.s1*100,x=(overallSurvival/365))) +geom_point(aes(shape=cls.GLM,color=cls.GLM),size=3.5,alpha=0.6)+
  scale_color_manual(values=c("red",'blue',"grey"))+scale_shape_manual(values=c(16,18,16))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                 legend.key=element_rect(fill="white"),axis.title = element_text(size=14),title = element_text(size=14),
                                                                                                 legend.text=element_text(size=14))
fig01.ninit<-BEAT_clinical_final_ninit %>% ggplot(aes(y=ridge.s1*100,x=(overallSurvival/365),color=TP53_mut_stat)) +geom_point(size=3,alpha=0.65)+
  scale_color_manual(values=c("red", "grey"))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                    panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                    legend.key=element_rect(fill="white"),axis.title = element_text(size=14),title = element_text(size=14),
                                                    legend.text=element_text(size=14))
fig02.ninit<-BEAT_clinical_final_ninit %>% ggplot(aes(y=ridge.s1*100,x=(overallSurvival/365))) +geom_point(aes(shape=cls.GLM,color=cls.GLM),size=3.5,alpha=0.6)+
  scale_color_manual(values=c("red",'blue',"grey"))+scale_shape_manual(values=c(16,18,16))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                 legend.key=element_rect(fill="white"),axis.title = element_text(size=14),title = element_text(size=14),
                                                                                                 legend.text=element_text(size=14))

fig03<-TCGA_clinical_final %>% ggplot(aes(y=ridge.s1*100,x=(OS_d/365),color=cls.GLM.y)) +geom_point(aes(shape=cls.GLM.y,color=cls.GLM.y),size=3,alpha=0.6)+
  scale_color_manual(values=c("red",'blue',"grey"))+scale_shape_manual(values=c(16,18,16))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                                                                 legend.key=element_rect(fill="white"),axis.title = element_text(size=14),title = element_text(size=14),
                                                                                                 legend.text=element_text(size=14))

setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis_initial AML/Ridge Regression/")
jpeg('BEAT_init_ridge_score_TP53mut.jpeg',width=5.2,height=3.8,units='in',res=300)
fig01.init
dev.off()

jpeg('BEAT_init_ridge_score_TP53mutlike.jpeg',width=5.2,height=3.8,units='in',res=300)
fig02.init
dev.off()

jpeg('BEAT_ninit_ridge_score_TP53mut.jpeg',width=5.2,height=3.8,units='in',res=300)
fig01.ninit
dev.off()

jpeg('BEAT_ninit_ridge_score_TP53mutlike.jpeg',width=5.2,height=3.8,units='in',res=300)
fig02.ninit
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