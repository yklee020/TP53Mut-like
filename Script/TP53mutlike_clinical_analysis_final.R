#2022.11.02
#Tyner group recently published paper in Cancer cell 2022 and it contains updated
#and extra data for BEAT AML. Clinical data for wave 1+2 is also updated

source("C:/Users/yklee/Desktop/Sachs Lab/R_project/TP53/TP53GEP_functions.r")
packages<-c('dplyr','tidyr','tidyverse',
            'ggplot2','survival','survminer','ComplexHeatmap')

for(p in packages){
  library(p,character.only = T)
}
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version")

BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')

#View(cbind(BEAT_clinical_final.mtwt.ELN.new$`ELN Risk Group_YL`,BEAT_clinical_final.mtwt.ELN.new$ELN2017))
fav_q_ix<-which(BEAT_clinical_final$`ELN Risk Group_YL`=='?Fav')
BEAT_clinical_final$ELN_Risk_Group_YL_new<-BEAT_clinical_final$`ELN Risk Group_YL`
#BEAT_clinical_final$ELN_Risk_Group_YL_new[fav_q_ix]<-BEAT_clinical_final$ELN2017[fav_q_ix]

for(i in (fav_q_ix)){
  cat<-BEAT_clinical_final$ELN2017[i]
  if(cat=='Favorable'){
    BEAT_clinical_final$ELN_Risk_Group_YL_new[i]<-'Favorable'
  }else if(cat=='Intermediate'){
    BEAT_clinical_final$ELN_Risk_Group_YL_new[i]<-'Intermediate'
  }else if(cat=='Adverse'){
    BEAT_clinical_final$ELN_Risk_Group_YL_new[i]<-'Adverse'
  }}
table(BEAT_clinical_final$ELN_Risk_Group_YL)
table(BEAT_clinical_final$ELN_Risk_Group_YL_new)
#use below
BEAT_clinical_final
TCGA_clinical_final

for (i in 1:dim(BEAT_clinical_final)[1]){
  cls<-names(table(BEAT_clinical_final$cls.GLM))
  if(BEAT_clinical_final$cls.GLM[i]=='TP53_WT'){#TP53wt
    BEAT_clinical_final$cls.GLM[i]<-'TP53wt'
  }
  if(BEAT_clinical_final$cls.GLM[i]=='TP53_MUT'){#TP53mut
    BEAT_clinical_final$cls.GLM[i]<-'TP53mut'
  }
  if(BEAT_clinical_final$cls.GLM[i]=='TP53_MUTlike'){#TP53mut-like
    BEAT_clinical_final$cls.GLM[i]<-'TP53mut-like'
  }
}
for (i in 1:dim(TCGA_clinical_final)[1]){
  cls<-names(table(TCGA_clinical_final$cls.GLM.y))
  if(TCGA_clinical_final$cls.GLM.y[i]=='TP53_WT'){#TP53wt
    TCGA_clinical_final$cls.GLM.y[i]<-'TP53wt'
  }
  if(TCGA_clinical_final$cls.GLM.y[i]=='TP53_MUT'){#TP53mut
    TCGA_clinical_final$cls.GLM.y[i]<-'TP53mut'
  }
  if(TCGA_clinical_final$cls.GLM.y[i]=='TP53MUT_like'){#TP53mut-like
    TCGA_clinical_final$cls.GLM.y[i]<-'TP53mut-like'
  }
}


#1. survival analysis for mut vs. wt and mutlike vs. wt

#BEAT wave 1+2
#dat.beat.mtwt<-BEAT_clinical_final[BEAT_clinical_final$cls.GLM!='TP53mut-like',]
dat.beat.mtwt<-BEAT_clinical_final
dat.beat.mtlwt<-BEAT_clinical_final[BEAT_clinical_final$cls.GLM!='TP53mut',]

surv.dat.beat.mtwt<-Surv(time = dat.beat.mtwt$overallSurvival,
                    event = dat.beat.mtwt$OS_status)
surv.dat.beat.mtlwt<-Surv(time = dat.beat.mtlwt$overallSurvival,
                         event = dat.beat.mtlwt$OS_status)
# surv.dat.beat.mtwt.fit<-survfit(surv.dat.beat.mtwt ~cls.GLM, data=dat.beat.mtwt)
surv.dat.beat.mtwt.fit<-survfit(surv.dat.beat.mtwt ~TP53_mut_stat, data=dat.beat.mtwt)
surv.dat.beat.mtlwt.fit<-survfit(surv.dat.beat.mtlwt ~cls.GLM, data=dat.beat.mtlwt)
# pdf('Analysis/Clinical/Survival analysis//BEAT_Mut_Mutlike_vs_wtothers.pdf')
# ggsurvplot(surv.dat.beat.mtwt.fit, data = dat.beat.mtwt, pval = TRUE)
# ggsurvplot(surv.dat.beat.mtlwt.fit, data = dat.beat.mtlwt, pval = TRUE)
# dev.off()

jpeg('Analysis/Clinical/Survival analysis/BEAT_Mut_vs_wtothers.jpeg',width=4.5,height=3.8,units='in',res=300)
ggsurvplot(surv.dat.beat.mtwt.fit, data = dat.beat.mtwt, pval = TRUE,palette=c('#00BFC4','#F8766D'))
dev.off()
jpeg('Analysis/Clinical/Survival analysis/BEAT_Mutlike_vs_wtothers.jpeg',width=4.5,height=3.8,units='in',res=300)
ggsurvplot(surv.dat.beat.mtlwt.fit, data = dat.beat.mtlwt, pval = TRUE)
dev.off()

jpeg('Analysis/Clinical/Survival analysis/BEAT_Mut_vs_wtothers_wopval_230118.jpeg',width=4.5,height=3.8,units='in',res=300)
ggsurvplot(surv.dat.beat.mtwt.fit, data = dat.beat.mtwt, pval = F,palette=c('#00BFC4','#F8766D'))
dev.off()
jpeg('Analysis/Clinical/Survival analysis/BEAT_Mutlike_vs_wtothers_wopval.jpeg',width=4.5,height=3.8,units='in',res=300)
ggsurvplot(surv.dat.beat.mtlwt.fit, data = dat.beat.mtlwt, pval = F)
dev.off()


#TCGA
# dat.tcga.mtwt<-TCGA_clinical_final[TCGA_clinical_final$cls.GLM.y!='TP53mut-like',]
dat.tcga.mtwt<-TCGA_clinical_final
dat.tcga.mtlwt<-TCGA_clinical_final[TCGA_clinical_final$cls.GLM.y!='TP53mut',]

surv.dat.tcga.mtwt<-Surv(time = dat.tcga.mtwt$OS_d,
                         event = dat.tcga.mtwt$OS_status)
surv.dat.tcga.mtlwt<-Surv(time = dat.tcga.mtlwt$OS_d,
                          event = dat.tcga.mtlwt$OS_status)
#surv.dat.tcga.mtwt.fit<-survfit(surv.dat.tcga.mtwt ~cls.GLM.y, data=dat.tcga.mtwt)
surv.dat.tcga.mtwt.fit<-survfit(surv.dat.tcga.mtwt ~TP53_mut_stat, data=dat.tcga.mtwt)
surv.dat.tcga.mtlwt.fit<-survfit(surv.dat.tcga.mtlwt ~cls.GLM.y, data=dat.tcga.mtlwt)

pdf('Analysis/Clinical/Survival analysis//TCGA_Mut_Mutlike_vs_wtothers.pdf')
ggsurvplot(surv.dat.tcga.mtwt.fit, data = dat.tcga.mtwt, pval = TRUE)
ggsurvplot(surv.dat.tcga.mtlwt.fit, data = dat.tcga.mtlwt, pval = TRUE)
dev.off()

jpeg('Analysis/Clinical/Survival analysis/TCGA_Mut_vs_wtothers.jpeg',width=4.5,height=3.8,units='in',res=300)
ggsurvplot(surv.dat.tcga.mtwt.fit, data = dat.tcga.mtwt, pval = TRUE)
dev.off()
jpeg('Analysis/Clinical/Survival analysis/TCGA_Mut_vs_wtothers_wopval_230118.jpeg',width=4.5,height=3.8,units='in',res=300)
ggsurvplot(surv.dat.tcga.mtwt.fit, data = dat.tcga.mtwt, pval = F)
dev.off()

jpeg('Analysis/Clinical/Survival analysis/TCGA_Mutlike_vs_wtothers.jpeg',width=4.5,height=3.8,units='in',res=300)
ggsurvplot(surv.dat.tcga.mtlwt.fit, data = dat.tcga.mtlwt, pval = F)
dev.off()

jpeg('Analysis/Clinical/Survival analysis/TCGA_Mutlike_vs_wtothers_wopavl.jpeg',width=4.5,height=3.8,units='in',res=300)
ggsurvplot(surv.dat.tcga.mtlwt.fit, data = dat.tcga.mtlwt, pval = F)
dev.off()



#generate Kaplan Meiere plots- all of them together
BEAT_clinical_final.n<-BEAT_clinical_final
TCGA_clinical_final.n<-TCGA_clinical_final
BEAT_clinical_final.n$cls.GLM<-relevel(factor(BEAT_clinical_final.n$cls.GLM),ref='TP53wt')

dat.beat.obj<-Surv(time = BEAT_clinical_final.n$overallSurvival,
                        event = BEAT_clinical_final.n$OS_status)

dat.tcga.obj<-Surv(time = TCGA_clinical_final$OS_d,
                   event = TCGA_clinical_final$OS_status)

dat.beat.fit<-survfit(dat.beat.obj ~cls.GLM, data=BEAT_clinical_final.n)
dat.tcga.fit<-survfit(dat.tcga.obj ~cls.GLM.y, data=TCGA_clinical_final)

beat.km.all<-ggsurvplot(dat.beat.fit, data = BEAT_clinical_final.n, pval = F)
tcga.km.all<-ggsurvplot(dat.tcga.fit, data = TCGA_clinical_final, pval = F)

beat.km.all$plot+scale_color_manual(values=c("#F6BA2A","#F8766D","#619CFF"))
tcga.km.all$plot+scale_color_manual(values=c("#F8766D","#619CFF","#F6BA2A"))

jpeg('Analysis/Clinical/Survival analysis/KM_BEAT_Mut_Mutl_wtothers_wopval.jpeg',width=4.5,height=3.8,units='in',res=300)
#pdf('Analysis/Clinical/Survival analysis/KM_BEAT_Mut_Mutl_wtothers_wopval.pdf',width=4.5,height=3.8)
print(beat.km.all$plot+scale_color_manual(values=c("#F6BA2A","#F8766D","#619CFF")))
dev.off()
jpeg('Analysis/Clinical/Survival analysis/KM_TCGAT_Mut_Mutl_vs_wtothers_wopval.jpeg',width=4.5,height=3.8,units='in',res=300)
#pdf('Analysis/Clinical/Survival analysis/KM_TCGAT_Mut_Mutl_vs_wtothers_wopval.pdf',width=4.5,height=3.8)
print(tcga.km.all$plot+scale_color_manual(values=c("#F8766D","#619CFF","#F6BA2A")))
dev.off()

#ELN classified survival analysis:
BEAT_clinical_final$cls.GLM<-relevel(factor(BEAT_clinical_final$cls.GLM),ref='TP53wt')
TCGA_clinical_final$cls.GLM.y<-relevel(factor(TCGA_clinical_final$cls.GLM.y),ref='TP53wt')
BEAT_clinical_final.ELN.fav<-BEAT_clinical_final[BEAT_clinical_final$`ELN Risk Group_YL`=='Favorable',]
BEAT_clinical_final.ELN.inter<-BEAT_clinical_final[BEAT_clinical_final$`ELN Risk Group_YL`=='Intermediate',]
BEAT_clinical_final.ELN.adv<-BEAT_clinical_final[BEAT_clinical_final$`ELN Risk Group_YL`=='Adverse',]
TCGA_clinical_final.ELN.fav<-TCGA_clinical_final[TCGA_clinical_final$`ELN Risk Group_YL`=='Favorable',]
TCGA_clinical_final.ELN.inter<-TCGA_clinical_final[TCGA_clinical_final$`ELN Risk Group_YL`=='Intermediate',]
TCGA_clinical_final.ELN.adv<-TCGA_clinical_final[TCGA_clinical_final$`ELN Risk Group_YL`=='Adverse',]


beat.fav.obj<-Surv(time = BEAT_clinical_final.ELN.fav$overallSurvival,
                   event = BEAT_clinical_final.ELN.fav$OS_status)
beat.inter.obj<-Surv(time = BEAT_clinical_final.ELN.inter$overallSurvival,
                     event = BEAT_clinical_final.ELN.inter$OS_status)
beat.adv.obj<-Surv(time = BEAT_clinical_final.ELN.adv$overallSurvival[BEAT_clinical_final.ELN.adv$cls.GLM!='TP53mut'],
                   event = BEAT_clinical_final.ELN.adv$OS_status[BEAT_clinical_final.ELN.adv$cls.GLM!='TP53mut'])

tcga.fav.obj<-Surv(time = TCGA_clinical_final.ELN.fav$OS_d,
                   event = TCGA_clinical_final.ELN.fav$OS_status)
tcga.inter.obj<-Surv(time = TCGA_clinical_final.ELN.inter$OS_d,
                     event = TCGA_clinical_final.ELN.inter$OS_status)
tcga.adv.obj<-Surv(time = TCGA_clinical_final.ELN.adv$OS_d[TCGA_clinical_final.ELN.adv$cls.GLM.y!='TP53mut'],
                   event = TCGA_clinical_final.ELN.adv$OS_status[TCGA_clinical_final.ELN.adv$cls.GLM.y!='TP53mut'])

beat.fav.fit<-survfit(beat.fav.obj ~cls.GLM, data=BEAT_clinical_final.ELN.fav)
beat.inter.fit<-survfit(beat.inter.obj ~cls.GLM, data=BEAT_clinical_final.ELN.inter)
beat.adv.fit<-survfit(beat.adv.obj ~cls.GLM, data=BEAT_clinical_final.ELN.adv[BEAT_clinical_final.ELN.adv$cls.GLM!='TP53mut',])

tcga.fav.fit<-survfit(tcga.fav.obj ~cls.GLM.y, data=TCGA_clinical_final.ELN.fav)
tcga.inter.fit<-survfit(tcga.inter.obj ~cls.GLM.y, data=TCGA_clinical_final.ELN.inter)
tcga.adv.fit<-survfit(tcga.adv.obj ~cls.GLM.y, data=TCGA_clinical_final.ELN.adv[TCGA_clinical_final.ELN.adv$cls.GLM.y!='TP53mut',])
# 
# pdf('Analysis/Clinical/Survival analysis/BEAT_ELN_KM.pdf')
# ggsurvplot(beat.fav.fit, data = BEAT_clinical_final.ELN.fav, pval = T,palette=c('#619CFF','#F8766D'),risk.table=F,title='BEAT_Favorable')
# ggsurvplot(beat.inter.fit, data = BEAT_clinical_final.ELN.inter, pval = T,palette=c('#619CFF','#F8766D'),risk.table=F,title='BEAT_Intermediate')
# ggsurvplot(beat.adv.fit, data = BEAT_clinical_final.ELN.adv[BEAT_clinical_final.ELN.adv$cls.GLM!='TP53_MUT',], pval = T,palette=c('#619CFF','#F8766D'),risk.table=F,title='BEAT_Adverse')
# dev.off()
# pdf('Analysis/Clinical/Survival analysis/TCGA_ELN_KM.pdf')
# ggsurvplot(tcga.fav.fit, data = TCGA_clinical_final.ELN.fav, pval = T,palette=c('#619CFF','#F8766D'),risk.table=F,title='TCGA_Favorable')
# ggsurvplot(tcga.inter.fit, data = TCGA_clinical_final.ELN.inter, pval = T,palette=c('#619CFF','#F8766D'),risk.table=F,title='TCGA_Intermediate')
# ggsurvplot(tcga.adv.fit, data = TCGA_clinical_final.ELN.adv[TCGA_clinical_final.ELN.adv$cls.GLM.y!='TP53_MUT',], pval = T,palette=c('#619CFF','#F8766D'),risk.table=F,title='TCGA_Adverse')
# dev.off()

#BEAT AML
jpeg('Analysis/Clinical/Survival analysis/ELN KM/BEAT_ELN_KM_fav.jpeg',width=5.5,height=4.5,units='in',res=300)
ggsurvplot(beat.fav.fit, data = BEAT_clinical_final.ELN.fav, pval = T,palette=c('#619CFF','#F8766D'),risk.table=F,title='BEAT_Favorable')
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN KM/BEAT_ELN_KM_inter.jpeg',width=5.5,height=4.5,units='in',res=300)
ggsurvplot(beat.inter.fit, data = BEAT_clinical_final.ELN.inter, pval = T,palette=c('#619CFF','#F8766D'),risk.table=F,title='BEAT_Intermediate')
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN KM/BEAT_ELN_KM_adv.jpeg',width=5.5,height=4.5,units='in',res=300)
ggsurvplot(beat.adv.fit, data = BEAT_clinical_final.ELN.adv[BEAT_clinical_final.ELN.adv$cls.GLM!='TP53mut',], pval = T,palette=c('#619CFF','#F8766D'),risk.table=F,title='BEAT_Adverse')
dev.off()
#wo/pval
jpeg('Analysis/Clinical/Survival analysis/ELN KM/BEAT_ELN_KM_fav_wopval.jpeg',width=5.5,height=4.5,units='in',res=300)
ggsurvplot(beat.fav.fit, data = BEAT_clinical_final.ELN.fav, pval = F,palette=c('#619CFF','#F8766D'),risk.table=F,title='BEAT_Favorable')
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN KM/BEAT_ELN_KM_inter_wopval.jpeg',width=5.5,height=4.5,units='in',res=300)
ggsurvplot(beat.inter.fit, data = BEAT_clinical_final.ELN.inter, pval = F,palette=c('#619CFF','#F8766D'),risk.table=F,title='BEAT_Intermediate')
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN KM/BEAT_ELN_KM_adv_wopval.jpeg',width=5.5,height=4.5,units='in',res=300)
ggsurvplot(beat.adv.fit, data = BEAT_clinical_final.ELN.adv[BEAT_clinical_final.ELN.adv$cls.GLM!='TP53_MUT',], pval = F,palette=c('#619CFF','#F8766D'),risk.table=F,title='BEAT_Adverse')
dev.off()


#TCGA
jpeg('Analysis/Clinical/Survival analysis/ELN KM/TCGA_ELN_KM_fav.jpeg',width=5.5,height=4.5,units='in',res=300)
ggsurvplot(tcga.fav.fit, data = TCGA_clinical_final.ELN.fav, pval = T,palette=c('#619CFF','#F8766D'),risk.table=F,title='TCGA_Favorable')
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN KM/TCGA_ELN_KM_inter.jpeg',width=5.5,height=4.5,units='in',res=300)
ggsurvplot(tcga.inter.fit, data = TCGA_clinical_final.ELN.inter, pval = T,palette=c('#619CFF','#F8766D'),risk.table=F,title='TCGA_Intermediate')
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN KM/TCGA_ELN_KM_adv.jpeg',width=5.5,height=4.5,units='in',res=300)
ggsurvplot(tcga.adv.fit, data = TCGA_clinical_final.ELN.adv[TCGA_clinical_final.ELN.adv$cls.GLM.y!='TP53mut',], pval = T,palette=c('#619CFF','#F8766D'),risk.table=F,title='TCGA_Adverse')
dev.off()
#wo/pval
jpeg('Analysis/Clinical/Survival analysis/ELN KM/TCGA_ELN_KM_fav_wopval.jpeg',width=5.5,height=4.5,units='in',res=300)
ggsurvplot(tcga.fav.fit, data = TCGA_clinical_final.ELN.fav, pval = F,palette=c('#619CFF','#F8766D'),risk.table=F,title='TCGA_Favorable')
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN KM/TCGA_ELN_KM_inter_wopval.jpeg',width=5.5,height=4.5,units='in',res=300)
ggsurvplot(tcga.inter.fit, data = TCGA_clinical_final.ELN.inter, pval = F,palette=c('#619CFF','#F8766D'),risk.table=F,title='TCGA_Intermediate')
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN KM/TCGA_ELN_KM_adv_wopval.jpeg',width=5.5,height=4.5,units='in',res=300)
ggsurvplot(tcga.adv.fit, data = TCGA_clinical_final.ELN.adv[TCGA_clinical_final.ELN.adv$cls.GLM.y!='TP53mut',], pval = F,palette=c('#619CFF','#F8766D'),risk.table=F,title='TCGA_Adverse')
dev.off()





BEAT_clinical_final.ELN.fav$cls.GLM<-relevel(droplevels(BEAT_clinical_final.ELN.fav$cls.GLM),ref='TP53wt')
BEAT_clinical_final.ELN.inter$cls.GLM<-relevel(droplevels(BEAT_clinical_final.ELN.inter$cls.GLM),ref='TP53wt')
#BEAT_clinical_final.ELN.adv$cls.GLM<-droplevels(BEAT_clinical_final.ELN.adv$cls.GLM)
BEAT_clinical_final.ELN.adv.sub<-BEAT_clinical_final.ELN.adv[BEAT_clinical_final.ELN.adv$cls.GLM!='TP53mut',]
BEAT_clinical_final.ELN.adv.sub$cls.GLM<-relevel(droplevels(BEAT_clinical_final.ELN.adv.sub$cls.GLM),ref='TP53wt')
beat.adv.sub.obj<-Surv(time = BEAT_clinical_final.ELN.adv.sub$overallSurvival,
                   event = BEAT_clinical_final.ELN.adv.sub$OS_status)


TCGA_clinical_final.ELN.fav$cls.GLM.y<-relevel(factor(TCGA_clinical_final.ELN.fav$cls.GLM.y),ref='TP53wt')
TCGA_clinical_final.ELN.inter$cls.GLM.y<-relevel(factor(TCGA_clinical_final.ELN.inter$cls.GLM.y),ref='TP53wt')
#TCGA_clinical_final.ELN.adv$cls.GLM.y<-relevel(factor(TCGA_clinical_final.ELN.adv$cls.GLM.y),ref='TP53wt')
TCGA_clinical_final.ELN.adv.sub<-TCGA_clinical_final.ELN.adv[TCGA_clinical_final.ELN.adv$cls.GLM.y!='TP53mut',]
TCGA_clinical_final.ELN.adv.sub$cls.GLM.y<-relevel(droplevels(TCGA_clinical_final.ELN.adv.sub$cls.GLM.y),ref='TP53wt')
tcga.adv.sub.obj<-Surv(time = TCGA_clinical_final.ELN.adv.sub$OS_d,
                       event = TCGA_clinical_final.ELN.adv.sub$OS_status)


beat.fav.coxfit<-coxph(beat.fav.obj ~cls.GLM, data=BEAT_clinical_final.ELN.fav)
beat.inter.coxfit<-coxph(beat.inter.obj ~cls.GLM, data=BEAT_clinical_final.ELN.inter)
beat.adv.coxfit<-coxph(beat.adv.sub.obj ~cls.GLM, data=BEAT_clinical_final.ELN.adv.sub)


tcga.fav.coxfit<-coxph(tcga.fav.obj ~cls.GLM.y, data=TCGA_clinical_final.ELN.fav)
tcga.inter.coxfit<-coxph(tcga.inter.obj ~cls.GLM.y, data=TCGA_clinical_final.ELN.inter)
tcga.adv.coxfit<-coxph(tcga.adv.sub.obj ~cls.GLM.y, data=TCGA_clinical_final.ELN.adv.sub)


ggforest3 <- function (model,
                       data = NULL,
                       main = "Hazard ratio",
                       cpositions = c(0.02,0.22, 0.4),
                       fontsize = 0.7,
                       refLabel = "reference",
                       noDigits = 2,
                       font.x.size = 20
                       )
{
  # dependencies
  require(broom)
  require(survival)
  require(grid)
  # model<-test
  # 
  # data = NULL
  # main = "Hazard ratio"
  # cpositions = c(0.02,0.22, 0.4)
  # fontsize = 0.7
  # refLabel = "reference"
  # noDigits = 2
  # font.x.size = 20
  .get_data <- function(fit, data = NULL, complain = TRUE) {
    if(is.null(data)){
      if (complain)
        warning ("The `data` argument is not provided. Data will be extracted from model fit.")
      data <- eval(fit$call$data)
      if (is.null(data))
        stop("The `data` argument should be provided either to ggsurvfit or survfit.")
    }
    data
  } # end dependencies
  
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  data <- .get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(tidy(model, conf.int = TRUE))
  gmodel <- glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if (terms[i] %in% c("factor", "character")) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    }
    else if (terms[i] == "numeric") {
      data.frame(
        var = var,
        Var1 = "",
        Freq = nrow(data),
        pos = 1
      )
    }
    else {
      vars = grep(paste0("^", var, "*."), coef$term, value = TRUE)
      data.frame(
        var = vars,
        Var1 = "",
        Freq = nrow(data),
        pos = seq_along(vars)
      )
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds,])[, c("var",
                                               "level",
                                               "N",
                                               "p.value",
                                               "estimate",
                                               "conf.low",
                                               "conf.high",
                                               "pos")]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[,
                                                              4], noDigits + 1), toShowExp)
  toShowExpClean$stars <- paste0(
    round(toShowExpClean$p.value,
          noDigits + 1),
    " ",
    ifelse(toShowExpClean$p.value <
             0.05, "*", ""),
    ifelse(toShowExpClean$p.value < 0.01,
           "*", ""),
    ifelse(toShowExpClean$p.value < 0.001, "*",
           "")
  )
  toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"],
                              " - ", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  
  #assign nmaes for the classes
  #test_cls_names
  #toShowExpClean$level<-test_cls_names
  
  
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1,]
  rangeb <-
    range(toShowExpClean$conf.low, toShowExpClean$conf.high,
          na.rm = TRUE)
  breaks <- axisTicks(rangeb / 2, log = TRUE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <-
    fontsize * as.numeric(convertX(unit(theme_get()$text$size,
                                        "pt"), "mm"))
    p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) +
      geom_rect(
        aes(
          xmin = seq_along(var) - 0.5,
          xmax = seq_along(var) +
            0.5,
          ymin = exp(rangeplot[1]),
          ymax = exp(rangeplot[2]),
          fill = ordered(seq_along(var) %% 2 + 1)
        )
      ) + scale_fill_manual(values = c("#FFFFFF33",'#FFFFFF33'),
                                     #"#00000033"),
                          guide = "none") + geom_point(pch = 15,
                                                       size = 4) + geom_errorbar(aes(ymin = exp(conf.low),
                                                                                     ymax = exp(conf.high)), width = 0.15) + geom_hline(yintercept = 1,
                                                                                                                                        linetype = 3) + coord_flip(ylim = exp(rangeplot)) +
    ggtitle(main) + scale_y_log10(
      name = "",
      labels = sprintf("%g",
                       breaks),
      expand = c(0.02, 0.02),
      breaks = breaks
    ) +
    theme_light() + theme(
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "none",
      panel.border = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      # Modified here
      axis.text.x = element_text(size=font.x.size),
      # End modification
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    # xlab("") + annotate(
    #   geom = "text",
    #   x = x_annotate,
    #   y = exp(y_variable),
    #   label = toShowExpClean$var,
    #   fontface = "bold",
    #   hjust = 0,
    #   size = annot_size_mm
    # ) + 
    annotate(
      geom = "text",
      x = x_annotate,
      y = exp(y_nlevel),
      hjust = 0,
      label = toShowExpClean$level,
      vjust = -0.1,
      size = annot_size_mm
    ) + annotate(
      geom = "text",
      x = x_annotate,
      y = exp(y_nlevel),
      label = toShowExpClean$N,
      fontface = "italic",
      hjust = 0,
      vjust = ifelse(toShowExpClean$level ==
                       "", 0.5, 1.1),
      size = annot_size_mm
    ) + annotate(
      geom = "text",
      x = x_annotate,
      y = exp(y_cistring),
      label = toShowExpClean$estimate.1,
      size = annot_size_mm,
      vjust = ifelse(toShowExpClean$estimate.1 ==
                       "reference", 0.5,-0.1)
    ) + annotate(
      geom = "text",
      x = x_annotate,
      y = exp(y_cistring),
      label = toShowExpClean$ci,
      size = annot_size_mm,
      vjust = 1.1,
      fontface = "italic"
    ) #+
      
    # annotate(
    #   geom = "text",
    #   x = x_annotate,
    #   y = exp(y_stars),
    #   label = toShowExpClean$stars,
    #   size = annot_size_mm,
    #   hjust = -0.2,
    #   fontface = "italic")
  #  + annotate(geom = "text",
  #   x = 0.5,y = exp(y_variable),
  #   label = paste0("# Events: ",gmodel$nevent,
  #     "; Global p-value (Log-Rank): ",format.pval(gmodel$p.value.log, eps = ".001"),
  #     " \nAIC: ",round(gmodel$AIC, 2),"; Concordance Index: ",
  #     round(gmodel$concordance,2)), size = annot_size_mm,
  #   hjust = 0,vjust = 1.2,fontface = "italic"
  # )
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}
# pdf('Analysis/Clinical/Survival analysis/BEAT_ELN_cox.pdf')
# ggforest3(beat.fav.coxfit,BEAT_clinical_final.ELN.fav,main='BEAT_FAV',fontsize = 1.5,font.x.size = 20)
# ggforest3(beat.inter.coxfit,BEAT_clinical_final.ELN.inter,main='BEAT_Int',fontsize = 1.5,font.x.size = 20)
# ggforest3(beat.adv.coxfit,BEAT_clinical_final.ELN.adv.sub,main='BEAT_Adv',fontsize = 1.5,font.x.size = 20)
# dev.off()
# pdf('Analysis/Clinical/Survival analysis/TCGA_ELN_cox.pdf')
# ggforest3(tcga.fav.coxfit,TCGA_clinical_final.ELN.fav,main='TCGA_FAV',fontsize = 1.5,font.x.size = 20)
# ggforest3(tcga.inter.coxfit,TCGA_clinical_final.ELN.inter,main='TCGA_Int',fontsize = 1.5,font.x.size = 20)
# ggforest3(tcga.adv.coxfit,TCGA_clinical_final.ELN.adv[TCGA_clinical_final.ELN.adv$cls.GLM.y!='TP53_MUT',],main='TCGA_Adv',fontsize = 1.5,font.x.size = 20)
# dev.off()
# 
# jpeg('Analysis/Clinical/Survival analysis/BEAT_ELN_Fav_cox.jpeg',width=5.2,height=3.8,units='in',res=300)
# ggforest3(beat.fav.coxfit,BEAT_clinical_final.ELN.fav,main='BEAT_FAV',fontsize = 0.8,font.x.size = 8)
# dev.off()
# jpeg('Analysis/Clinical/Survival analysis/BEAT_ELN_Int_cox.jpeg',width=5.2,height=3.8,units='in',res=300)
# ggforest3(beat.inter.coxfit,BEAT_clinical_final.ELN.inter,main='BEAT_Int',fontsize = 0.8,font.x.size = 8)
# dev.off()
# jpeg('Analysis/Clinical/Survival analysis/BEAT_ELN_Adv_cox.jpeg',width=5.2,height=3.8,units='in',res=300)
# ggforest3(beat.adv.coxfit,BEAT_clinical_final.ELN.adv.sub,main='BEAT_Adv',fontsize = 0.8,font.x.size = 8)
# dev.off()
# 
# jpeg('Analysis/Clinical/Survival analysis/TCGA_ELN_Fav_cox.jpeg',width=5.2,height=3.8,units='in',res=300)
# ggforest3(tcga.fav.coxfit,TCGA_clinical_final.ELN.fav,main='TCGA_FAV',fontsize = 0.8,font.x.size = 8)
# dev.off()
# jpeg('Analysis/Clinical/Survival analysis/TCGA_ELN_Int_cox.jpeg',width=5.2,height=3.8,units='in',res=300)
# ggforest3(tcga.inter.coxfit,TCGA_clinical_final.ELN.inter,main='TCGA_Int',fontsize = 0.8,font.x.size = 8)
# dev.off()
# jpeg('Analysis/Clinical/Survival analysis/TCGA_ELN_Adv_cox.jpeg',width=5.2,height=3.8,units='in',res=300)
# ggforest3(tcga.adv.coxfit,TCGA_clinical_final.ELN.adv[TCGA_clinical_final.ELN.adv$cls.GLM.y!='TP53_MUT',],main='TCGA_Adv',fontsize = 0.8,font.x.size = 8)
# dev.off()


jpeg('Analysis/Clinical/Survival analysis/ELN COX/BEAT_ELN_Fav_cox.jpeg',width=6.3,height=2,units='in',res=300)
ggforest3(beat.fav.coxfit,BEAT_clinical_final.ELN.fav,main='BEAT_FAV',fontsize = 1.75,font.x.size = 11,
         cpositions = c(0, 0, 0.35))
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN COX/BEAT_ELN_Int_cox.jpeg',width=6.3,height=2,units='in',res=300)
ggforest3(beat.inter.coxfit,BEAT_clinical_final.ELN.inter,main='BEAT_Int',fontsize = 1.75,font.x.size = 10,
          cpositions = c(0, 0, 0.35))
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN COX/BEAT_ELN_Adv_cox.jpeg',width=6.3,height=2,units='in',res=300)
ggforest3(beat.adv.coxfit,BEAT_clinical_final.ELN.adv.sub,main='BEAT_Adv',fontsize = 1.75,font.x.size = 10,
         cpositions = c(0, 0, 0.35))
dev.off()

jpeg('Analysis/Clinical/Survival analysis/ELN COX/TCGA_ELN_Fav_cox.jpeg',width=6.3,height=2,units='in',res=300)
ggforest3(tcga.fav.coxfit,TCGA_clinical_final.ELN.fav,main='TCGA_FAV',fontsize = 1.75,font.x.size = 10,
          cpositions = c(0, 0, 0.35))
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN COX/TCGA_ELN_Int_cox.jpeg',width=6.3,height=2,units='in',res=300)
ggforest3(tcga.inter.coxfit,TCGA_clinical_final.ELN.inter,main='TCGA_Int',fontsize = 1.75,font.x.size = 10,
          cpositions = c(0, 0, 0.35))
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN COX/TCGA_ELN_Adv_cox.jpeg',width=6.3,height=2,units='in',res=300)
ggforest3(tcga.adv.coxfit,TCGA_clinical_final.ELN.adv.sub,main='TCGA_Adv',fontsize = 1.75,font.x.size = 10,
          cpositions = c(0, 0, 0.35))
dev.off()

#final version: 2023.1.6

jpeg('Analysis/Clinical/Survival analysis/ELN COX/N_BEAT_ELN_Fav_cox.jpeg',width=6.3,height=2,units='in',res=300)
ggforest3(beat.fav.coxfit,BEAT_clinical_final.ELN.fav,main='BEAT_FAV',fontsize = 0,font.x.size = 10,
          cpositions = c(0, 0, 0))
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN COX/N_BEAT_ELN_Int_cox.jpeg',width=6.3,height=2,units='in',res=300)
ggforest3(beat.inter.coxfit,BEAT_clinical_final.ELN.inter,main='BEAT_Int',fontsize = 0,font.x.size = 10,
          cpositions = c(0, 0, 0))
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN COX/N_BEAT_ELN_Adv_cox.jpeg',width=6.3,height=2,units='in',res=300)
ggforest3(beat.adv.coxfit,BEAT_clinical_final.ELN.adv.sub,main='BEAT_Adv',fontsize = 0,font.x.size = 10,
          cpositions = c(0, 0, 0))
dev.off()

jpeg('Analysis/Clinical/Survival analysis/ELN COX/N_TCGA_ELN_Fav_cox.jpeg',width=6.3,height=2,units='in',res=300)
ggforest3(tcga.fav.coxfit,TCGA_clinical_final.ELN.fav,main='TCGA_FAV',fontsize = 0,font.x.size = 10,
          cpositions = c(0, 0, 0))
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN COX/N_TCGA_ELN_Int_cox.jpeg',width=6.3,height=2,units='in',res=300)
ggforest3(tcga.inter.coxfit,TCGA_clinical_final.ELN.inter,main='TCGA_Int',fontsize = 0,font.x.size = 10,
          cpositions = c(0, 0, 0))
dev.off()
jpeg('Analysis/Clinical/Survival analysis/ELN COX/N_TCGA_ELN_Adv_cox.jpeg',width=6.3,height=2,units='in',res=300)
ggforest3(tcga.adv.coxfit,TCGA_clinical_final.ELN.adv.sub,main='TCGA_Adv',fontsize = 0,font.x.size = 10,
          cpositions = c(0, 0, 0))
dev.off()

#get the pie chart, or stacked barchart or side by side pie chart
#just use excel



#check for biallelic or mono
cebpa_ix<-which(!is.na(BEAT_clinical_final.mtwt.ELN.new$CEBPA_Biallelic))
BEAT_clinical_final.mtwt.ELN.new$cls.GLM[which(BEAT_clinical_final.mtwt.ELN.new$CEBPA_Biallelic=='bi')]
BEAT_clinical_final.mtwt.ELN.new$cls.GLM[which(BEAT_clinical_final.mtwt.ELN.new$CEBPA_Biallelic=='mono')]
#check for TCGA
#####################################################################################

BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')

#View(cbind(BEAT_clinical_final.mtwt.ELN.new$`ELN Risk Group_YL`,BEAT_clinical_final.mtwt.ELN.new$ELN2017))
fav_q_ix<-which(BEAT_clinical_final$`ELN Risk Group_YL`=='?Fav')
BEAT_clinical_final$ELN_Risk_Group_YL_new<-BEAT_clinical_final$`ELN Risk Group_YL`
#BEAT_clinical_final$ELN_Risk_Group_YL_new[fav_q_ix]<-BEAT_clinical_final$ELN2017[fav_q_ix]

for(i in (fav_q_ix)){
  cat<-BEAT_clinical_final$ELN2017[i]
  if(cat=='Favorable'){
    BEAT_clinical_final$ELN_Risk_Group_YL_new[i]<-'Favorable'}
  # }else if(cat=='Intermediate'){
  #   BEAT_clinical_final$ELN_Risk_Group_YL_new[i]<-'Intermediate'
  # }else if(cat=='Adverse'){
  #   BEAT_clinical_final$ELN_Risk_Group_YL_new[i]<-'Adverse'
  # }
}
#Only for MT, WT-like_MT (no unknown) and WT
#assign
BEAT_LASSO_final_res.clinical.n<-BEAT_clinical_final

#hold on to the unknown included data on the side and work on MT and WT only
BEAT_LASSO_final_res.clinical.num<-BEAT_clinical_final[,c(1,2,4,31,71,77,80,81,84,85,102)]
#BEAT_LASSO_final_res.clinical.num.old<-BEAT_clinical_final.sub[,c(7,8,10,13,15,16,17,18,20,59,60,61)]


library(stringr)
library(dplyr)
BEAT_LASSO_final_res.clinical.num.pre<-BEAT_LASSO_final_res.clinical.num %>%
  mutate_at(vars(c("%.Blasts.in.BM","%.Blasts.in.PB" )), ~str_replace(.,"N/A",'NA')) %>% 
  mutate_at(vars(c("%.Blasts.in.BM" ,"%.Blasts.in.PB" )), ~str_replace(.,'\\.*>|\\.*<','')) %>%
  mutate_at(vars(c("%.Blasts.in.BM" ,"%.Blasts.in.PB" )), as.numeric) #%>%
#mutate_at(vars(nm), ~replace_na(., 0))

# BEAT_LASSO_final_res.clinical.num.pre.old<-BEAT_LASSO_final_res.clinical.num.old %>%mutate_at(vars('PatientId'), as.character) %>%
#   mutate_at(vars(c("%.Blasts.in.PB" )), ~str_replace(.,'rare','0.01')) %>% 
#   mutate_at(vars(c("%.Blasts.in.BM" ,"%.Blasts.in.PB" )), ~str_replace(.,'\\.*>|\\.*<','')) %>%
#   mutate_at(vars(c("%.Blasts.in.BM" ,"%.Blasts.in.PB" )), as.numeric)

#do ANOVA
library(tidyverse)
library(rstatix)
library(ggpubr)
library(dplyr)
num_res.GLM <- BEAT_LASSO_final_res.clinical.num.pre  %>% 
  group_by(cls.GLM) %>% dplyr::summarise(across(where(is.numeric),list(mean=mean,median=median),na.rm=T))
# num_res.GLM.old <- BEAT_LASSO_final_res.clinical.num.pre.old %>% 
#   group_by(cls.GLM) %>% summarise(across(where(is.numeric),list(mean=mean,median=median),na.rm=T))

colnames(BEAT_LASSO_final_res.clinical.num)

anova_num_res<-function(num_dat.pre,test_group_name){
  test.group<-as.symbol(test_group_name)
  num_cat<-colnames(num_dat.pre)[c(3:11)]
  p.out.df<-as.data.frame(matrix(0,3,length(num_cat)))
  colnames(p.out.df)<-num_cat
  beat.num.aov01<-aov(LSC_score.beat~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  beat.num.aov02<-aov(ageAtDiagnosis~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  beat.num.aov03<-aov(responseDurationToInductionTx~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  beat.num.aov04<-aov(overallSurvival~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  beat.num.aov05<-aov(`%.Blasts.in.BM`~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  beat.num.aov06<-aov(`%.Blasts.in.PB`~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  beat.num.aov07<-aov(`%.Lymphocytes.in.PB`~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  beat.num.aov08<-aov(`%.Monocytes.in.PB`~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  beat.num.aov09<-aov(`wbcCount`~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  rownames(p.out.df)<-rownames(beat.num.aov01[[1]])
  for (i in 1:length(num_cat)){
    p.out.df[,i]<-get(paste('beat.num.aov','0',i,sep=''))[[1]][,4]###########
  }
  return(p.out.df)
}

num_res.cls.GLM.aovp<-anova_num_res(BEAT_LASSO_final_res.clinical.num.pre,'cls.GLM')

#make it for pairwise wilcoxon test as well.
BEAT_LASSO_final_res.clinical.num.pre

#do pairwise-t.test for only 2 groups MUT vs. WT and MUTlike vs. WT and then p.adjust with BH

ttest_num_res<-function(num_dat.pre){
  #test.group<-as.symbol(test_group_name)
  num_cat<-colnames(num_dat.pre)[c(3:11)]
  p.out.df<-as.data.frame(matrix(0,2,length(num_cat)))
  colnames(p.out.df)<-num_cat
  beat.num.t01<-pairwise.t.test(num_dat.pre$LSC_score.beat,num_dat.pre$cls.GLM, p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t02<-pairwise.t.test(num_dat.pre$ageAtDiagnosis,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F)
  beat.num.t03<-pairwise.t.test(num_dat.pre$responseDurationToInductionTx,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t04<-pairwise.t.test(num_dat.pre$overallSurvival,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t05<-pairwise.t.test(num_dat.pre$`%.Blasts.in.BM`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t06<-pairwise.t.test(num_dat.pre$`%.Blasts.in.PB`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t07<-pairwise.t.test(num_dat.pre$`%.Lymphocytes.in.PB`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t08<-pairwise.t.test(num_dat.pre$`%.Monocytes.in.PB`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t09<-pairwise.t.test(num_dat.pre$`wbcCount`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  #c(beat.num.t01$p.value[1,1],beat.num.t01$p.value[2,2])
  rownames(p.out.df)<-c(paste(rownames(beat.num.t01$p.value)[2],colnames(beat.num.t01$p.value)[2],sep=':'),
                        paste(rownames(beat.num.t01$p.value)[2],colnames(beat.num.t01$p.value)[1],sep=':'))
  for (i in 1:length(num_cat)){
    t<-get(paste('beat.num.t','0',i,sep=''))
    p.cor<-p.adjust(c(t$p.value[2,2],t$p.value[2,1]),method='BH')
    p.out.df[,i]<-p.cor###########
  }
  return(p.out.df)
}
ttest_num_res_all_comp<-function(num_dat.pre){
  num_cat<-colnames(num_dat.pre)[c(3:11)]
  p.out.df<-as.data.frame(matrix(0,3,length(num_cat)))
  colnames(p.out.df)<-num_cat
  beat.num.t01<-pairwise.t.test(num_dat.pre$LSC_score.beat,num_dat.pre$cls.GLM, p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t02<-pairwise.t.test(num_dat.pre$ageAtDiagnosis,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F)
  beat.num.t03<-pairwise.t.test(num_dat.pre$responseDurationToInductionTx,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t04<-pairwise.t.test(num_dat.pre$overallSurvival,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t05<-pairwise.t.test(num_dat.pre$`%.Blasts.in.BM`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t06<-pairwise.t.test(num_dat.pre$`%.Blasts.in.PB`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t07<-pairwise.t.test(num_dat.pre$`%.Lymphocytes.in.PB`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t08<-pairwise.t.test(num_dat.pre$`%.Monocytes.in.PB`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  beat.num.t09<-pairwise.t.test(num_dat.pre$`wbcCount`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  #c(beat.num.t01$p.value[1,1],beat.num.t01$p.value[2,2])
  rownames(p.out.df)<-c(paste(rownames(beat.num.t01$p.value)[2],colnames(beat.num.t01$p.value)[2],sep=':'),
                        paste(rownames(beat.num.t01$p.value)[2],colnames(beat.num.t01$p.value)[1],sep=':'),
                        paste(rownames(beat.num.t01$p.value)[1],colnames(beat.num.t01$p.value)[1],sep=':'))
  for (i in 1:length(num_cat)){
    t<-get(paste('beat.num.t','0',i,sep=''))
    p.cor<-p.adjust(c(t$p.value[2,2],t$p.value[2,1],t$p.value[1,1]),method='BH')
    p.out.df[,i]<-p.cor###########
  }
  return(p.out.df)
}
num_res.cls.GLM.ttest<-ttest_num_res(BEAT_LASSO_final_res.clinical.num.pre)
num_res.cls.GLM.ttest.all<-ttest_num_res_all_comp(BEAT_LASSO_final_res.clinical.num.pre)

t<-pairwise.t.test(BEAT_LASSO_final_res.clinical.num.pre$wbcCount,BEAT_LASSO_final_res.clinical.num.pre$cls.GLM, p.adj = "none",na.rm=T,
                                       paired=FALSE,pool.sd = F)    
t<-pairwise.t.test(BEAT_LASSO_final_res.clinical.num.pre$`%.Blasts.in.BM`,BEAT_LASSO_final_res.clinical.num.pre$cls.GLM, p.adj = "none",na.rm=T,
                   paired=FALSE,pool.sd = F,)    
pair_wilcoxon_num_res_all<-function(num_dat.pre,test_group_name){
  test.group<-(test_group_name)
  num_cat<-colnames(num_dat.pre)[c(3:11)]
  p.out.df<-as.data.frame(matrix(0,3,length(num_cat)))
  colnames(p.out.df)<-num_cat
  beat.num.wilcox01<-pairwise.wilcox.test(num_dat.pre$LSC_score.beat,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  beat.num.wilcox02<-pairwise.wilcox.test(num_dat.pre$ageAtDiagnosis,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  beat.num.wilcox03<-pairwise.wilcox.test(num_dat.pre$responseDurationToInductionTx,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  beat.num.wilcox04<-pairwise.wilcox.test(num_dat.pre$overallSurvival,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  beat.num.wilcox05<-pairwise.wilcox.test(num_dat.pre$`%.Blasts.in.BM`,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  beat.num.wilcox06<-pairwise.wilcox.test(num_dat.pre$`%.Blasts.in.PB`,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  beat.num.wilcox07<-pairwise.wilcox.test(num_dat.pre$`%.Lymphocytes.in.PB`,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  beat.num.wilcox08<-pairwise.wilcox.test(num_dat.pre$`%.Monocytes.in.PB`,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  beat.num.wilcox09<-pairwise.wilcox.test(num_dat.pre$`wbcCount`,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  rownames(p.out.df)<-c(paste(colnames(beat.num.wilcox01$p.value)[2],rownames(beat.num.wilcox01$p.value)[2],sep=':'),
                        paste(colnames(beat.num.wilcox01$p.value)[1],rownames(beat.num.wilcox01$p.value)[2],sep=':'),
                        paste(rownames(beat.num.wilcox01$p.value)[1],colnames(beat.num.wilcox01$p.value)[1],sep=':'))
  for (i in 1:length(num_cat)){
    w<-get(paste('beat.num.wilcox','0',i,sep=''))
    p.cor<-p.adjust(c(w$p.value[2,2],w$p.value[2,1],w$p.value[1,1]),method='BH')
    p.out.df[,i]<-p.cor
  }
  
  return(p.out.df)
}
num_res.cls.GLM.wilcox.all<-pair_wilcoxon_num_res_all(BEAT_LASSO_final_res.clinical.num.pre,'cls.GLM')
pairwise.wilcox.test(BEAT_LASSO_final_res.clinical.num.pre$wbcCount,BEAT_LASSO_final_res.clinical.num.pre[,'cls.GLM'],p.adjust.method = 'BH')

write.csv(num_res.GLM,'Analysis/Clinical/Numerical parameter analysis/BEAT_clinical_numerical_num_res.csv')
write.csv(num_res.cls.GLM.aovp,'Analysis/Clinical/Numerical parameter analysis/BEAT_clinical_numerical_num_aov_FDR.csv')
write.csv(num_res.cls.GLM.ttest,'Analysis/Clinical/Numerical parameter analysis/BEAT_clinical_numerical_num_ttest_BH.csv')
write.csv(num_res.cls.GLM.ttest.all,'Analysis/Clinical/Numerical parameter analysis/BEAT_clinical_numerical_num_ttest_all_comp_BH_230118.csv')
write.csv(num_res.cls.GLM.wilcox.all,'Analysis/Clinical/Numerical parameter analysis/BEAT_clinical_numerical_num_wilcox_all_comp_BH_230120.csv')



####################################################
#TCGA

TCGA_clinical.final.num<-TCGA_clinical_final[,c(1,16:19,2048,2056,2066,2093,2095,2097,2101)]
#ttt<-TCGA_clinical.final.n[,c(1,30,15,16,17,18,2047,2055,2065,2095,2097:2098)]

TCGA_clinical.final.num.pre<-TCGA_clinical.final.num %>%
  mutate_at(vars('lab_procedure_abnormal_lymphocyte_result_percent_value'), as.numeric) #%>%

#double check this results one by one
num_res.GLM.tcga <- TCGA_clinical.final.num.pre %>% 
  group_by(cls.GLM.y) %>% dplyr::summarise(across(where(is.numeric),list(mean=mean,median=median),na.rm=T))

colnames(TCGA_clinical.final.num.pre)
anova_num_res_tcga<-function(num_dat.pre,test_group_name){
  test.group<-as.symbol(test_group_name)
  num_cat<-colnames(num_dat.pre)[c(2:11)]
  p.out.df<-as.data.frame(matrix(0,3,length(num_cat)))
  colnames(p.out.df)<-num_cat
  tcga.num.aov01<-aov(Age~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  tcga.num.aov02<-aov(`%BM Blast`~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  tcga.num.aov03<-aov(WBC~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  tcga.num.aov04<-aov(`%PB Blast`~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  tcga.num.aov05<-aov(lab_procedure_abnormal_lymphocyte_result_percent_value~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  tcga.num.aov06<-aov(lab_procedure_bone_marrow_lymphocyte_outcome_percent_value~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  tcga.num.aov07<-aov(lab_procedure_monocyte_result_percent_value~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  tcga.num.aov08<-aov(OS_d~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  tcga.num.aov09<-aov(EFS_d~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  tcga.num.aov10<-aov(LSC_score.tcga~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  # tcga.num.aov11<-aov(LSC_score.comb.tcga~eval(test.group), data=num_dat.pre) %>% TukeyHSD()
  num_id<-c('01','02','03','04','05','06','07','08','09','10')
  
  rownames(p.out.df)<-rownames(tcga.num.aov01[[1]])
  for (i in 1:length(num_cat)){
    p.out.df[,i]<-get(paste('tcga.num.aov',num_id[i],sep=''))[[1]][,4]###########
  }
  return(p.out.df)
}
aov(`%PB Blast`~cls.GLM.y, data=TCGA_clinical.final.num.pre) %>% TukeyHSD()
aov(`%BM Blast`~cls.GLM.y, data=TCGA_clinical.final.num.pre) %>% TukeyHSD()

num_res.cls.GLM.tcga.aovp<-anova_num_res_tcga(TCGA_clinical.final.num.pre,'cls.GLM.y')

ttest_num_res<-function(num_dat.pre){
  #test.group<-as.symbol(test_group_name)
  num_cat<-colnames(num_dat.pre)[c(2:11)]
  p.out.df<-as.data.frame(matrix(0,2,length(num_cat)))
  colnames(p.out.df)<-num_cat
  tcga.num.t01<-pairwise.t.test(num_dat.pre$Age,num_dat.pre$cls.GLM, p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t02<-pairwise.t.test(num_dat.pre$`%BM Blast`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F)
  tcga.num.t03<-pairwise.t.test(num_dat.pre$WBC,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t04<-pairwise.t.test(num_dat.pre$`%PB Blast`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t05<-pairwise.t.test(num_dat.pre$lab_procedure_abnormal_lymphocyte_result_percent_value,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t06<-pairwise.t.test(num_dat.pre$lab_procedure_bone_marrow_lymphocyte_outcome_percent_value,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t07<-pairwise.t.test(num_dat.pre$lab_procedure_monocyte_result_percent_value,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t08<-pairwise.t.test(num_dat.pre$OS_d,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t09<-pairwise.t.test(num_dat.pre$EFS_d,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t10<-pairwise.t.test(num_dat.pre$LSC_score.tcga,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  # tcga.num.t11<-pairwise.t.test(num_dat.pre$LSC_score.comb.tcga,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
  #                               paired=FALSE,pool.sd = F) 
  #c(beat.num.t01$p.value[1,1],beat.num.t01$p.value[2,2])
  rownames(p.out.df)<-c(paste(colnames(tcga.num.t01$p.value)[1],rownames(tcga.num.t01$p.value)[1],sep=':'),
                        paste(rownames(tcga.num.t01$p.value)[2],colnames(tcga.num.t01$p.value)[2],sep=':'))
  num_id<-c('01','02','03','04','05','06','07','08','09','10')
  
  #rownames(p.out.df)<-rownames(tcga.num.t01[[1]])
  for (i in 1:length(num_cat)){
    t<-get(paste('tcga.num.t',num_id[i],sep=''))
    p.cor<-p.adjust(c(t$p.value[1,1],t$p.value[2,2]),method='BH')
    p.out.df[,i]<-p.cor
  }
  return(p.out.df)
}

ttest_num_res_all_comp<-function(num_dat.pre){
  #test.group<-as.symbol(test_group_name)
  num_cat<-colnames(num_dat.pre)[c(2:11)]
  p.out.df<-as.data.frame(matrix(0,3,length(num_cat)))
  colnames(p.out.df)<-num_cat
  tcga.num.t01<-pairwise.t.test(num_dat.pre$Age,num_dat.pre$cls.GLM, p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t02<-pairwise.t.test(num_dat.pre$`%BM Blast`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F)
  tcga.num.t03<-pairwise.t.test(num_dat.pre$WBC,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t04<-pairwise.t.test(num_dat.pre$`%PB Blast`,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t05<-pairwise.t.test(num_dat.pre$lab_procedure_abnormal_lymphocyte_result_percent_value,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t06<-pairwise.t.test(num_dat.pre$lab_procedure_bone_marrow_lymphocyte_outcome_percent_value,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t07<-pairwise.t.test(num_dat.pre$lab_procedure_monocyte_result_percent_value,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t08<-pairwise.t.test(num_dat.pre$OS_d,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t09<-pairwise.t.test(num_dat.pre$EFS_d,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  tcga.num.t10<-pairwise.t.test(num_dat.pre$LSC_score.tcga,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
                                paired=FALSE,pool.sd = F) 
  # tcga.num.t11<-pairwise.t.test(num_dat.pre$LSC_score.comb.tcga,num_dat.pre$cls.GLM,p.adj = "none",na.rm=T,
  #                               paired=FALSE,pool.sd = F) 
  #c(beat.num.t01$p.value[1,1],beat.num.t01$p.value[2,2])
  rownames(p.out.df)<-c(paste(colnames(tcga.num.t01$p.value)[1],rownames(tcga.num.t01$p.value)[1],sep=':'),
                        paste(rownames(tcga.num.t01$p.value)[2],colnames(tcga.num.t01$p.value)[2],sep=':'),
                        paste(rownames(tcga.num.t01$p.value)[2],colnames(tcga.num.t01$p.value)[1],sep=':'))
  num_id<-c('01','02','03','04','05','06','07','08','09','10')
  
  #rownames(p.out.df)<-rownames(tcga.num.t01[[1]])
  for (i in 1:length(num_cat)){
    t<-get(paste('tcga.num.t',num_id[i],sep=''))
    p.cor<-p.adjust(c(t$p.value[1,1],t$p.value[2,2],t$p.value[2,1]),method='BH')
    p.out.df[,i]<-p.cor
  }
  return(p.out.df)
}
t.test(TCGA_clinical.final.num.pre$WBC[TCGA_clinical.final.num.pre$cls.GLM.y=='TP53_MUT'],
       TCGA_clinical.final.num.pre$WBC[TCGA_clinical.final.num.pre$cls.GLM.y=='TP53_WT'])

t<-pairwise.t.test(TCGA_clinical.final.num.pre$WBC,TCGA_clinical.final.num.pre$cls.GLM, p.adj = "none",na.rm=T,
                   paired=FALSE,pool.sd = F) 
p.adjust(c(t$p.value[1,1],t$p.value[2,2]),method='BH')
t<-pairwise.t.test(TCGA_clinical.final.num.pre$`%BM Blast`,TCGA_clinical.final.num.pre$cls.GLM.y, p.adj = "none",na.rm=T,
                   paired=FALSE,pool.sd = F) 
p.adjust(c(t$p.value[1,1],t$p.value[2,2]),method='BH')

num_res.cls.GLM.tcga.ttest<-ttest_num_res(TCGA_clinical.final.num.pre)
num_res.cls.GLM.tcga.ttest.all<-ttest_num_res_all_comp(TCGA_clinical.final.num.pre)



pair_wilcoxon_num_res_tcga<-function(num_dat.pre,test_group_name){
  test.group<-(test_group_name)
  num_cat<-colnames(num_dat.pre)[c(2:11)]
  p.out.df<-as.data.frame(matrix(0,3,length(num_cat)))
  colnames(p.out.df)<-num_cat
  tcga.num.wilcox01<-pairwise.wilcox.test(num_dat.pre$Age,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  tcga.num.wilcox02<-pairwise.wilcox.test(num_dat.pre$`%BM Blast`,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  tcga.num.wilcox03<-pairwise.wilcox.test(num_dat.pre$WBC,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  tcga.num.wilcox04<-pairwise.wilcox.test(num_dat.pre$`%PB Blast`,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  tcga.num.wilcox05<-pairwise.wilcox.test(num_dat.pre$lab_procedure_abnormal_lymphocyte_result_percent_value,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  tcga.num.wilcox06<-pairwise.wilcox.test(num_dat.pre$lab_procedure_bone_marrow_lymphocyte_outcome_percent_value,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  tcga.num.wilcox07<-pairwise.wilcox.test(num_dat.pre$lab_procedure_monocyte_result_percent_value,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  tcga.num.wilcox08<-pairwise.wilcox.test(num_dat.pre$OS_d,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  tcga.num.wilcox09<-pairwise.wilcox.test(num_dat.pre$EFS_d,num_dat.pre[,test.group],
                                           p.adjust.method = 'none')
  tcga.num.wilcox10<-pairwise.wilcox.test(num_dat.pre$LSC_score.tcga,num_dat.pre[,test.group],
                                          p.adjust.method = 'none')
  rownames(p.out.df)<-c(paste(colnames(beat.num.wilcox01$p.value)[2],rownames(beat.num.wilcox01$p.value)[2],sep=':'),
                        paste(colnames(beat.num.wilcox01$p.value)[1],rownames(beat.num.wilcox01$p.value)[2],sep=':'),
                        paste(rownames(beat.num.wilcox01$p.value)[1],colnames(beat.num.wilcox01$p.value)[1],sep=':'))
  num_id<-c('01','02','03','04','05','06','07','08','09','10')
  for (i in 1:length(num_cat)){
    w<-get(paste('tcga.num.wilcox',num_id[i],sep=''))
    p.cor<-p.adjust(c(w$p.value[2,2],w$p.value[1,1],w$p.value[2,1]),method='BH')
    p.out.df[,i]<-p.cor
  }
  return(p.out.df)
}
num_res.cls.GLM.tcga.wilcox.all<-pair_wilcoxon_num_res_tcga(TCGA_clinical.final.num.pre,'cls.GLM.y')
pairwise.wilcox.test(TCGA_clinical.final.num.pre$WBC,TCGA_clinical.final.num.pre[,'cls.GLM.y'],
                     p.adjust.method = 'BH')
write.csv(num_res.GLM.tcga,'Analysis/Clinical/Numerical parameter analysis/TCGA_clinical_numerical_num_res.csv')
write.csv(num_res.cls.GLM.tcga.aovp,'Analysis/Clinical/Numerical parameter analysis/TCGA_clinical_numerical_num_aov_FDR.csv')
write.csv(num_res.cls.GLM.tcga.ttest,'Analysis/Clinical/Numerical parameter analysis/TCGA_clinical_numerical_num_ttest_BH.csv')
write.csv(num_res.cls.GLM.tcga.ttest.all,'Analysis/Clinical/Numerical parameter analysis/TCGA_clinical_numerical_num_ttest_all_comp_BH_230118.csv')
write.csv(num_res.cls.GLM.tcga.wilcox.all,'Analysis/Clinical/Numerical parameter analysis/TCGA_clinical_numerical_num_wilcox_all_comp_BH_230120.csv')

#############################################################################
#get dot plots for numerical variables

#draw dot plots for BM% and WBC in BEAT and TCGA:2022.08.26
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Clinical/Numerical parameter analysis/Figures/")

# TCGA_clinical_final<-readRDS('class ids/TCGA_clinical_final_class.rds')
# BEAT_clinical_final<-readRDS('class ids/BEAT_clinical_final_class.rds')
# BEAT_clinical_final.sub<-BEAT_clinical_final[BEAT_clinical_final$cls.GLM!='TP53_unk_WT_like',]
# BEAT_clinical_final.sub$cls.GLM[BEAT_clinical_final.sub$cls.GLM=='TP53_WTunk_MT_like']<-'TP53MUT_like'
# TCGA_clinical_final$cls.GLM[TCGA_clinical_final$cls.GLM=='TP53_unk_MT_like']<-'TP53MUT_like'

BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')

library(ggpubr)
library(rstatix)
library(tibble)
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
BEAT_clinical_final$cls.GLM<-factor(BEAT_clinical_final$cls.GLM,levels=c('TP53_MUT','TP53_MUTlike','TP53_WT'))
BEAT_clinical_final$`%.Blasts.in.BM`<-as.numeric(BEAT_clinical_final$`%.Blasts.in.BM`)
TCGA_clinical_final$cls.GLM.y[TCGA_clinical_final$cls.GLM.y=='TP53MUT_like']<-'TP53_MUTlike'
TCGA_clinical_final$cls.GLM.y<-factor(TCGA_clinical_final$cls.GLM.y,levels=c('TP53_MUT','TP53_MUTlike','TP53_WT'))

df.d<-BEAT_clinical_final %>% filter(!is.na(`%.Blasts.in.BM`))
df.d$BM<-as.numeric(df.d$`%.Blasts.in.BM`)
df.wbc<-BEAT_clinical_final %>% filter(!is.na(wbcCount))
df.wbc$wbcCount<-as.numeric(df.wbc$wbcCount)
df.d$LSC_score.beat
df.d$ageAtDiagnosis
df.BM.tcga<-TCGA_clinical_final
df.BM.tcga$BM<-df.BM.tcga$`%BM Blast`
df.wbc.tcga<-TCGA_clinical_final
df.BM.tcga$LSC_score.tcga
df.BM.tcga$age_at_initial_pathologic_diagnosis
# log transformation is not a good idea - make a multiple lines figure
# df.d$BM.log<-log(df.d$BM)
# df.wbc$wbcCount.log<-log(df.wbc$wbcCount)
# df.BM.tcga$BM.log<-log(df.BM.tcga$BM)
# df.wbc.tcga$WBC.log<-log(df.wbc.tcga$WBC)

stat.test <- df.d %>%
  #group_by(cls.GLM) %>%
  t_test(BM~cls.GLM) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
p.adjust(stat.test[c(2,3),]$p,'BH')
stat.wbc <- df.wbc %>%
  #group_by(cls.GLM) %>%
  t_test(wbcCount~cls.GLM) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
p.adjust(stat.wbc[c(2,3),]$p,'BH')

stat.BM.tcga <- df.BM.tcga %>%
  #group_by(cls.GLM) %>%
  t_test(BM~cls.GLM.y) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
p.adjust(stat.BM.tcga[c(2,3),]$p,'BH')

stat.wbc.tcga <- df.wbc.tcga %>%
  #group_by(cls.GLM) %>%
  t_test(WBC~cls.GLM.y) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
p.adjust(stat.wbc.tcga[c(2,3),]$p,'BH')

df.d
df2.bm<-data_summary(df.d,varname='BM',groupnames=c('cls.GLM'))
df2.wbc<-data_summary(df.wbc,varname='wbcCount',groupnames=c('cls.GLM'))
df2.lsc<-data_summary(df.d,varname='LSC_score.beat',groupnames=c('cls.GLM'))
df2.age<-data_summary(df.d,varname='ageAtDiagnosis',groupnames=c('cls.GLM'))
df2.bm.tcga<-data_summary(df.BM.tcga,varname='BM',groupnames=c('cls.GLM.y'))
df2.wbc.tcga<-data_summary(df.wbc.tcga,varname='WBC',groupnames=c('cls.GLM.y'))
df2.lsc.tcga<-data_summary(df.BM.tcga,varname='LSC_score.tcga',groupnames=c('cls.GLM.y'))
df2.age.tcga<-data_summary(df.BM.tcga,varname='age_at_initial_pathologic_diagnosis',groupnames=c('cls.GLM.y'))
#######################################
beat.bm<-ggdotplot(df.d,x='cls.GLM',y='BM',add = 'mean_se',color='black',fill='black',
                   position=position_jitter(0.225),binwidth = 2,alpha=0.5)+ 
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+
  theme(axis.text=element_text(size=16))+scale_y_continuous(breaks=seq(0,119,20))#ylim(0,100)#expand_limits(y=c(0,100))
beat.bm.wo.ast<-ggdotplot(df.d,x='cls.GLM',y='BM',add = 'mean_se',color='black',fill='black',
                   position=position_jitter(0.225),binwidth = 2,alpha=0.5)+ 
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("", "", "", "", "")),
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+
  theme(axis.text=element_text(size=16))+scale_y_continuous(breaks=seq(0,119,20))
#label ='p.signif',
# beat.bm<-ggplot(df2.bm,aes(x=cls.GLM,y=BM,group=cls.GLM)) +
#   geom_bar(stat='identity',fill='steelblue',width=0.5)+theme_minimal()+
#   geom_errorbar(aes(x=cls.GLM, ymin=BM-se, ymax=BM+se), width=0.1, colour="orange", alpha=0.9, size=0.3) + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))  
df.wbc%>%filter(wbcCount<400)
beat.wbc<-ggdotplot(df.wbc%>%filter(wbcCount<400),x='cls.GLM',y='wbcCount',add = 'mean_se',color='black',fill='black',size=3,
                    position=position_jitter(0.2),binwidth = 2,alpha=0.5)+
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  #geom_text(x=3.25, y=310, label="WBC>400 is removed(n=1,wt)",size=3)+
  font("xlab", size = 14)+font("ylab", size = 18)+theme(axis.text=element_text(size=16))

beat.wbc.wo.ast<-ggdotplot(df.wbc%>%filter(wbcCount<400),x='cls.GLM',y='wbcCount',add = 'mean_se',color='black',fill='black',size=3,
                    position=position_jitter(0.2),binwidth = 2,alpha=0.5)+
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("", "", "", "", "")),
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  #geom_text(x=3.25, y=310, label="WBC>400 is removed(n=1,wt)",size=3)+
  font("xlab", size = 14)+font("ylab", size = 18)+theme(axis.text=element_text(size=16))

# # beat.bm<-ggplot(df2.bm,aes(x=cls.GLM,y=BM,group=cls.GLM)) +
# #   geom_bar(stat='identity',fill='steelblue',width=0.5)+theme_minimal()+
# #   geom_errorbar(aes(x=cls.GLM, ymin=BM-se, ymax=BM+se), width=0.1, colour="orange", alpha=0.9, size=0.3) + 
# #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
# #         panel.background = element_blank(), axis.line = element_line(colour = "black"))  
# 
# beat.wbc<-ggplot(df2.wbc,aes(x=cls.GLM,y=WBC.Count,group=cls.GLM)) +
#   geom_bar(stat='identity',fill='steelblue',width=0.5)+theme_minimal()+
#   geom_errorbar( aes(x=cls.GLM, ymin=WBC.Count-se, ymax=WBC.Count+se), width=0.1, colour="orange", alpha=0.9, size=0.3)

beat.lsc<-ggdotplot(df.d,x='cls.GLM',y='LSC_score.beat',add = 'mean_se',color='black',fill='black',
                     position=position_jitter(0.25),alpha=0.5,size=1, binwidth = 1/35)+ 
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+
  theme(axis.text=element_text(size=16))

beat.lsc.wo.ast<-ggdotplot(df.d,x='cls.GLM',y='LSC_score.beat',add = 'mean_se',color='black',fill='black',
                           position=position_jitter(0.25),alpha=0.5,size=1, binwidth = 1/35)+ 
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("", "", "", "", "")),
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+
  theme(axis.text=element_text(size=16))
beat.age<-ggdotplot(df.d,x='cls.GLM',y='ageAtDiagnosis',add = 'mean_se',color='black',fill='black',
                    position=position_jitter(0.25),alpha=0.5,size=1, binwidth = 1/35)+ 
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+
  theme(axis.text=element_text(size=16))

beat.age.wo.ast<-ggdotplot(df.d,x='cls.GLM',y='ageAtDiagnosis',add = 'mean_se',color='black',fill='black',
                           position=position_jitter(0.25),alpha=0.5,size=1, binwidth = 1/35)+ 
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("", "", "", "", "")),
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+
  theme(axis.text=element_text(size=16))
##########
tcga.bm<-ggdotplot(df.BM.tcga,x='cls.GLM.y',y='BM',add = 'mean_se',color='black',fill='black',
                   position=position_jitter(0.225),binwidth = 2,alpha=0.5)+ 
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,hide.ns = T,
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+theme(axis.text=element_text(size=16))+expand_limits(y=0)+
  scale_y_continuous(breaks=seq(0,119,20))

tcga.bm.wo.ast<-ggdotplot(df.BM.tcga,x='cls.GLM.y',y='BM',add = 'mean_se',color='black',fill='black',
                   position=position_jitter(0.225),binwidth = 2,alpha=0.5)+ 
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("", "", "", "", "")),
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+theme(axis.text=element_text(size=16))+expand_limits(y=0)+
  scale_y_continuous(breaks=seq(0,119,20))


tcga.wbc<-ggdotplot(df.wbc.tcga,x='cls.GLM.y',y='WBC',add = 'mean_se',color='black',fill='black',size=3,
                    position=position_jitter(0.2),binwidth = 2,alpha=0.5)+
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+theme(axis.text=element_text(size=16))

tcga.wbc.wo.ast<-ggdotplot(df.wbc.tcga,x='cls.GLM.y',y='WBC',add = 'mean_se',color='black',fill='black',size=3,
                    position=position_jitter(0.2),binwidth = 2,alpha=0.5)+
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("", "", "", "", "")),
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+theme(axis.text=element_text(size=16))

tcga.lsc<-ggdotplot(df.BM.tcga,x='cls.GLM.y',y='LSC_score.tcga',add = 'mean_se',color='black',fill='black',
                    position=position_jitter(0.25),alpha=0.5,size=1, binwidth = 1/35)+ 
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,hide.ns = T,
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+theme(axis.text=element_text(size=16))

tcga.lsc.wo.ast<-ggdotplot(df.BM.tcga,x='cls.GLM.y',y='LSC_score.tcga',add = 'mean_se',color='black',fill='black',
                           position=position_jitter(0.25),alpha=0.5,size=1, binwidth = 1/35)+ 
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,hide.ns = T,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("", "", "", "", "")),
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+theme(axis.text=element_text(size=16))
#  geom_text(x=3.25, y=310, label="WBC>400 is removed(n=1,wt)",size=3)
tcga.age<-ggdotplot(df.BM.tcga,x='cls.GLM.y',y='age_at_initial_pathologic_diagnosis',add = 'mean_se',color='black',fill='black',
                    position=position_jitter(0.25),alpha=0.5,size=1, binwidth = 1/35)+ 
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,hide.ns = T,
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+theme(axis.text=element_text(size=16))

tcga.age.wo.ast<-ggdotplot(df.BM.tcga,x='cls.GLM.y',y='age_at_initial_pathologic_diagnosis',add = 'mean_se',color='black',fill='black',
                           position=position_jitter(0.25),alpha=0.5,size=1, binwidth = 1/35)+ 
  stat_summary(shape=95,fun='mean',geom='point',color='red',size=18)+
  stat_summary(fun.data='mean_se',geom='errorbar',color='red',width=.05)+
  stat_compare_means(method='t.test',p.adjust.method='BH',label ='p.signif',bracket.size=0.5,hide.ns = T,
                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("", "", "", "", "")),
                     comparisons = list(c('TP53_MUT','TP53_WT'),c('TP53_MUTlike','TP53_WT')))+
  font("xlab", size = 14)+font("ylab", size = 18)+theme(axis.text=element_text(size=16))


age_at_initial_pathologic_diagnosis
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Clinical/Numerical parameter analysis/")

jpeg('BEAT_BM_dot.jpeg',width=5.2,height=3.8,units='in',res=300)
beat.bm
dev.off()
jpeg('BEAT_wbc_dot.jpeg',width=5.2,height=3.8,units='in',res=300)
beat.wbc
dev.off()
jpeg('TCGA_BM_dot.jpeg',width=5.2,height=3.8,units='in',res=300)
tcga.bm
dev.off()
jpeg('TCGA_WBC_dot.jpeg',width=5.2,height=3.8,units='in',res=300)
tcga.wbc
dev.off()


jpeg('BEAT_BM_dot_wo_ast.jpeg',width=5.2,height=3.8,units='in',res=300)
beat.bm.wo.ast
dev.off()
jpeg('BEAT_wbc_dot_wo_ast.jpeg',width=5.2,height=3.8,units='in',res=300)
beat.wbc.wo.ast
dev.off()
jpeg('TCGA_BM_dot_wo_ast.jpeg',width=5.2,height=3.8,units='in',res=300)
tcga.bm.wo.ast
dev.off()
jpeg('TCGA_WBC_dot_wo_ast.jpeg',width=5.2,height=3.8,units='in',res=300)
tcga.wbc.wo.ast
dev.off()



jpeg('BEAT_LSC_dot.jpeg',width=5.2,height=3.8,units='in',res=300)
beat.lsc
dev.off()
jpeg('BEAT_LSC_dot_wo_ast.jpeg',width=5.2,height=3.8,units='in',res=300)
beat.lsc.wo.ast
dev.off()

jpeg('TCGA_LSC_dot.jpeg',width=5.2,height=3.8,units='in',res=300)
tcga.lsc
dev.off()
jpeg('TCGA_LSC_dot_wo_ast.jpeg',width=5.2,height=3.8,units='in',res=300)
tcga.lsc.wo.ast
dev.off()

# ggpar(beat.bm, ylim = c(0, 100))
# ggpar(tcga.bm, ylim = c(0, 100))

jpeg('BEAT_BM_dot.jpeg',width=1024,height=640)
ggpar(beat.bm, ylim = c(0, 100))
dev.off()
jpeg('BEAT_wbc_dot.jpeg',width=1024,height=640)
beat.wbc
dev.off()
jpeg('TCGA_BM_dot.jpeg',width=1024,height=640)
ggpar(tcga.bm, ylim = c(0, 100))
dev.off()
jpeg('TCGA_WBC_dot.jpeg',width=1024,height=640)
tcga.wbc
dev.off()


setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Clinical/Numerical parameter analysis/")


jpeg('BEAT_age_dot.jpeg',width=5.2,height=3.8,units='in',res=300)
beat.age
dev.off()
jpeg('BEAT_age_dot_wo_aster.jpeg',width=5.2,height=3.8,units='in',res=300)
beat.age.wo.ast
dev.off()
jpeg('TCGA_age_dot.jpeg',width=5.2,height=3.8,units='in',res=300)
tcga.age
dev.off()
jpeg('TCGA_age_dot_wo_aster.jpeg',width=5.2,height=3.8,units='in',res=300)
tcga.age.wo.ast
dev.off()
t.test(df.d$ageAtDiagnosis[df.d$cls.GLM=='TP53_MUT'],df.d$ageAtDiagnosis[df.d$cls.GLM=='TP53_WT'])

################################
#Categorical parameters
#1.Mutational status
#2. Cytogenetics : ELN and other delq and all
#3. other categorical data
BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")

library(dplyr)
#TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/TCGA_clinical_final_class.rds')
BEAT_clinical_final.old<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/BEAT_clinical_final_class.rds')
BEAT_clinical_final.old_403<-BEAT_clinical_final.old[match(BEAT_clinical_final$LabId,BEAT_clinical_final.old$LabId),]
#mutational status
#need to preprocess mutational subtypes FLT3,PTPN11,SF3B1(21:23) for 1,0 and NA
table(BEAT_clinical_final.old_403$FLT3)
table(BEAT_clinical_final.old_403$WES_PTPN11)
table(BEAT_clinical_final.old_403$WES_SF3B1)

restind<-BEAT_clinical_final.old_403$responseToInductionTx
restind[which(restind=='UNKNOWN')]<-'Unknown'
restind[which(restind=='Complete Response i')]<-'Complete Response'
flt3_itd<-ifelse(BEAT_clinical_final.old_403$`FLT3-ITD`=='positive',1,
                 ifelse(BEAT_clinical_final.old_403$`FLT3-ITD`=='negative',0,NA))
npm1<-ifelse(BEAT_clinical_final.old_403$NPM1=='positive',1,
             ifelse(BEAT_clinical_final.old_403$NPM1=='negative',0,NA))
# flt3<-ifelse(BEAT_LASSO_final_res.clinical.n$FLT3=='negative',0,
#              ifelse(BEAT_LASSO_final_res.clinical.n$FLT3!='negative',1,NA))
# ptpn11<-ifelse(is.na(BEAT_LASSO_final_res.clinical.n$PTPN11),NA,1)
# sf3b1<-ifelse(is.na(BEAT_LASSO_final_res.clinical.n$SF3B1),NA,1)
BEAT_clinical_final.old_403$responseToInductionTx<-restind
BEAT_clinical_final.old_403$`FLT3-ITD`<-flt3_itd
BEAT_clinical_final.old_403$NPM1<-npm1
# BEAT_LASSO_final_res.clinical.n$FLT3<-flt3
# BEAT_LASSO_final_res.clinical.n$PTPN11<-ptpn11
# BEAT_LASSO_final_res.clinical.n$SF3B1<-sf3b1

BEAT_LASSO_final_res_cat<-BEAT_clinical_final.old_403[,c(4,6,9,11,12,14,19,21:58)]
BEAT_LASSO_final_res_cat<-right_join(BEAT_LASSO_final_res_cat,BEAT_clinical_final[,c(1,2,111)],by='LabId')
BEAT_LASSO_final_res_cat[,c(8:44)]<-data.frame(lapply(BEAT_LASSO_final_res_cat[,c(8:44)],as.character))
BEAT_LASSO_final_cat<-BEAT_LASSO_final_res_cat[,-c(8:44)]
# for(i in 8:44){
#   BEAT_LASSO_final_res_cat[,i]<-ifelse(BEAT_LASSO_final_res_cat[,i]==1,'MUT','WT')
# }
# cat_res_long<-BEAT_LASSO_final_res_cat %>% group_by(new_group) %>%
#   summarize(across(c(2:11),list(table),na.rm=T))
Mut_list<-c('DNMT3A','TET2','IDH1','IDH2','WT1','EZH2','ASXL1','SRSF2','STAG2','RAD21',
            'FLT3','KIT','KRAS','NRAS','PTPN11','NF1','JAK2','BCOR','GATA2','ETV6',
            'CEBPA','RUNX1','CDKN2A','NPM1','MYC','U2AF1','PHF6','SMC3','SMC1A','CSF3R',
            'PDS5B','SF3B1','BRINP3','CROCC','PLCE1','SPEN')
Mut_list.n<-str_c('WES',Mut_list,sep='_')
#use clinical summary data for NPM1 and FLT3-ITD - these two are complete
Mut_list.n[Mut_list.n =='WES_NPM1']<-'NPM1'
# Mut_list.n[Mut_list.n =='TP53_mut_stat']<-'TP53'
Mut_list.n<-c(Mut_list.n,'FLT3-ITD')

library(reshape2)
# shp<-melt(BEAT_LASSO_final_res_cat,measure.vars = c("isRelapse","responseToInductionTx","typeInductionTx",'vitalStatus',
#                                                     "FAB/Blast.Morphology",Mut_list.n,'ZS_anno_group'))#"TP53_mut_stat"
BEAT_LASSO_final_res_cat$TP53<-ifelse(BEAT_LASSO_final_res_cat$TP53=='TP53_MUT',1,0)
shp<-melt(BEAT_LASSO_final_res_cat,measure.vars = c(Mut_list.n,"TP53"))#"TP53_mut_stat"

cat_res<-shp %>% group_by(cls.GLM,variable,value) %>% dplyr::summarize(n=n(),na.rm=TRUE) %>% mutate(freq=n/sum(n,na.rm=T))

#make an 'export' variable
cat_res$export <- with(cat_res, sprintf("%i (%.1f%%)", n, freq*100))
#reshape again
output <- dcast(variable+value~cls.GLM, value.var="export", data=cat_res, fill="missing") #use drop=F to prevent silent missings 
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


write.csv(output2,'Analysis/Clinical/Categorical parameter analysis/BEAT_AML403_Mutations.csv')
#save.image('GLM_CV_GEP005.Rdata')#one that before using cv$fit.preval for train data score
#load('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GEP/GLM_CV_GEP005.Rdata')

#do chi-square/pairwise-fisher exact testing for one vairable for now
#how to extract count data?make a separate data to do counting data

cat_res_t<-shp %>% group_by(cls.GLM,variable,value) %>% dplyr::summarize(n=n(),na.rm=TRUE) 
#cat_res_t$export <- with(cat_res_t, sprintf("%i (%.1f%%)", n, freq*100))
#use NA then 0 - for the testing purpose
output_t <- dcast(variable+value~cls.GLM, value.var="n", data=cat_res_t, fill=0) #use drop=F to prevent silent missings 

#output_t.n<-output_t[-c(6,12,15,27,76,106),]
output_t.n<-output_t
#Mut_list.n
output_t.n.mut<-output_t.n[output_t.n$variable %in% c(Mut_list.n,'TP53'),]#or 24:95

# output_t.n.cat<-output_t.n[output_t.n$variable %in% c("isRelapse","responseToInductionTx","typeInductionTx",'vitalStatus',
#                                                       "FAB/Blast.Morphology",'ZS_anno_group'),]
#Mut first
Cat_group<-c('TP53_MUT','TP53_MUTlike','TP53_WT')
Cat_var<-unique(output_t.n.mut$variable)
var_num<-seq(1:length(Cat_var))
xmut<-vector(mode='list',length=length(Cat_var))
names(xmut)<-Cat_var
for(i in 1:length(Cat_var)){
  ix<-which(output_t.n.mut$variable%in% Cat_var[i])
  xcat.sub<-as.table(as.matrix(output_t.n.mut[ix,3:5]))
  dimnames(xcat.sub)<-list(output_t.n.mut[ix,2],
                           Group=Cat_group)
  xmut[[i]]<-xcat.sub
}
library(rstatix)
mut_ref<-pairwise_fisher_test(xmut[[1]],p.adjust.method = "none")
mut_ref$p


BEAT_mut_stat.mat<-matrix(0,2,length(xmut))
rownames(BEAT_mut_stat.mat)<-str_c(mut_ref$group1,mut_ref$group2,sep = ':')[c(2:3)]
colnames(BEAT_mut_stat.mat)<-names(xmut)
for (i in 1:length(xmut)){
  tryCatch({
    print(i)
    BEAT_mut_stat.mat[,i]<-pairwise_fisher_test(xmut[[i]])$p[c(2:3)] %>% p.adjust(method='BH')  
  },error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
  
  #eval(as.symbol(paste('xcat',num_list[i],'_res',sep='')))$p.adj
}

output_t.n.mut$variable <- as.character(output_t.n.mut$variable)
#make 'empty lines' 
empties.mut <- data.frame(variable=unique(output_t.n.mut$variable), stringsAsFactors=F)
empties.mut[,colnames(output_t.n.mut)[-1]] <- ""
output2.mut <- rbind(empties.mut,output_t.n.mut)
output2.mut <- output2.mut[order(output2.mut$variable,output2.mut$value),]
output2.mut$variable[output2.mut$value!=""] <- ""

write.csv(output2.mut,'Analysis/Clinical/Categorical parameter analysis/BEAT_AML403_Mutations_no_freq.csv')
write.csv(BEAT_mut_stat.mat,'Analysis/Clinical/Categorical parameter analysis/BEAT_AML403_Mutations_fisher_tests.csv')

#3 combination comparisons

BEAT_mut_stat.mat<-matrix(0,3,length(xmut))
rownames(BEAT_mut_stat.mat)<-str_c(mut_ref$group1,mut_ref$group2,sep = ':')[c(2:3,1)]
colnames(BEAT_mut_stat.mat)<-names(xmut)
for (i in 1:length(xmut)){
  tryCatch({
    print(i)
    BEAT_mut_stat.mat[,i]<-pairwise_fisher_test(xmut[[i]])$p[c(2:3,1)] %>% p.adjust(method='BH')  
  },error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
  
  #eval(as.symbol(paste('xcat',num_list[i],'_res',sep='')))$p.adj
}

output_t.n.mut$variable <- as.character(output_t.n.mut$variable)
#make 'empty lines' 
empties.mut <- data.frame(variable=unique(output_t.n.mut$variable), stringsAsFactors=F)
empties.mut[,colnames(output_t.n.mut)[-1]] <- ""
output2.mut <- rbind(empties.mut,output_t.n.mut)
output2.mut <- output2.mut[order(output2.mut$variable,output2.mut$value),]
output2.mut$variable[output2.mut$value!=""] <- ""

#write.csv(output2.mut,'Analysis/Clinical/Categorical parameter analysis/BEAT_AML403_Mutations_no_freq.csv')
write.csv(BEAT_mut_stat.mat,'Analysis/Clinical/Categorical parameter analysis/BEAT_AML403_Mutations_fisher_tests_3comb.csv')

#################
#get TCGA clinical data - not many data to look at
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')
TCGA_clinical.final.n<-TCGA_clinical_final


#categorical results
#make it long form
Mut_list<-c('DNMT3A','TET2','IDH1','IDH2','WT1','EZH2','ASXL1','SRSF2','STAG2','RAD21',
                 'FLT3','KIT','KRAS','NRAS','PTPN11','NF1','JAK2','BCOR','GATA2','ETV6',
                 'CEBPA','RUNX1','CDKN2A','NPM1','MYC','U2AF1','PHF6','SMC3','SMC1A','CSF3R',
                 'PDS5B','SF3B1','BRINP3','CROCC','PLCE1','SPEN')
#clinical.rna.GDC.178.final.mut<-TCGA_clinical.final.n[,which(colnames(TCGA_clinical.final.n)%in%Mut_list)]

dim(TCGA_clinical.final.n)

wes.ix01<-which(c(Mut_list,'TP53')%in%colnames(TCGA_clinical.final.n) )
wes.ix01.id<-c(Mut_list,'TP53')[wes.ix01]
df01<-as.data.frame(matrix(0,dim(TCGA_clinical.final.n)[1],length(wes.ix01.id)))

rownames(df01)<-TCGA_clinical.final.n$patient_RNA_id
colnames(df01)<-wes.ix01.id

for(i in 1:length((wes.ix01.id))){
  ix<-which(colnames(TCGA_clinical.final.n) %in% wes.ix01.id[i])
  mut_stat<-TCGA_clinical.final.n[,ix]
  df01[,i]<-ifelse(is.na(mut_stat),0,1)
}


#BEAT_LASSO_final_res_cat<-BEAT_LASSO_final_res.clinical.n[,c(1,5,8,10,11,18,20:24,26)]
# TCGA_clinical.final.cat<-TCGA_clinical.final.n[,c(1,14,25,27,35,2067,2098,2099:2102)]
TCGA_clinical.final.cat<-TCGA_clinical.final.n[,c(1,15,26,27,36,2014,2068,2088,2101,2102)]
TCGA_clinical.final.cat.mut<-right_join(TCGA_clinical.final.cat,
                                        rownames_to_column(df01,var='patient_RNA_id'),by='patient_RNA_id')
TCGA_clinical.final.mut<-TCGA_clinical.final.cat.mut[,-c(2:8,10)]
TCGA_clinical.final.cat<-TCGA_clinical.final.cat.mut[,c(2:8,10)]

shp.tcga.mut<-melt(TCGA_clinical.final.mut,measure.vars = colnames(TCGA_clinical.final.mut)[c(3:37)])
output_mut.tcga<-shp.tcga.mut %>% group_by(cls.GLM.y,variable,value) %>% dplyr::summarize(n=n(),na.rm=TRUE) %>% mutate(freq=n/sum(n,na.rm=T))

#use NA then 0 - for the testing purpose
output_mut.tcga$export <- with(output_mut.tcga, sprintf("%i (%.1f%%)", n, freq*100))

#output_mut <- dcast(variable+value~cls.GLM.y, value.var="n", data=mut_res.tcga, fill=0) #use drop=F to prevent silent missings 
output_mut.tcga <- dcast(variable+value~cls.GLM.y, value.var="export", data=output_mut.tcga, fill=0)

output_mut.tcga$variable <- as.character(output_mut.tcga$variable)
#make 'empty lines' 
empties <- data.frame(variable=unique(output_mut.tcga$variable), stringsAsFactors=F)
empties[,colnames(output_mut.tcga)[-1]] <- ""

#bind them together
output2.tcga <- rbind(empties,output_mut.tcga)
output2.tcga <- output2.tcga[order(output2.tcga$variable,output2.tcga$value),]

#optional: 'remove' variable if value present
output2.tcga$variable[output2.tcga$value!=""] <- ""
write.csv(output2.tcga,'Analysis/Clinical/Categorical parameter analysis/TCGA_Mutations.csv')


# shp.tcga.mut<-melt(TCGA_clinical.final.mut,measure.vars = colnames(TCGA_clinical.final.mut)[c(3:37)])
output_mut.tcga<-shp.tcga.mut %>% group_by(cls.GLM.y,variable,value) %>% dplyr::summarize(n=n(),na.rm=TRUE)# %>% mutate(freq=n/sum(n,na.rm=T))
output_mut.tcga <- dcast(variable+value~cls.GLM.y, value.var="n", data=output_mut.tcga, fill=0)


#Mut first
Cat_group<-c('TP53_MUT','TP53_WT','TP53MUT_like')
Cat_var<-unique(output_mut.tcga$variable)
var_num<-seq(1:length(Cat_var))
xmut<-vector(mode='list',length=length(Cat_var))
names(xmut)<-Cat_var
for(i in 1:length(Cat_var)){
  ix<-which(output_mut.tcga$variable%in% Cat_var[i])
  xcat.sub<-as.table(as.matrix(output_mut.tcga[ix,3:5]))
  dimnames(xcat.sub)<-list(output_mut.tcga[ix,2],
                           Group=Cat_group)
  xmut[[i]]<-xcat.sub
}

mut_ref<-pairwise_fisher_test(xmut[[1]])

#2 group comparisons
TCGA_mut_stat.mat<-matrix(0,2,length(xmut))
rownames(TCGA_mut_stat.mat)<-str_c(mut_ref$group1,mut_ref$group2,sep = ':')[c(1,3)]
colnames(TCGA_mut_stat.mat)<-names(xmut)
for (i in 1:length(xmut)){
  tryCatch({
    TCGA_mut_stat.mat[,i]<-pairwise_fisher_test(xmut[[i]])$p[c(1,3)] %>% p.adjust(method='BH')  
  },error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
  
  #eval(as.symbol(paste('xcat',num_list[i],'_res',sep='')))$p.adj
}

output_mut.tcga$variable <- as.character(output_mut.tcga$variable)
#make 'empty lines' 
empties.mut <- data.frame(variable=unique(output_mut.tcga$variable), stringsAsFactors=F)
empties.mut[,colnames(output_mut.tcga)[-1]] <- ""
output2.mut <- rbind(empties.mut,output_mut.tcga)
output2.mut <- output2.mut[order(output2.mut$variable,output2.mut$value),]
output2.mut$variable[output2.mut$value!=""] <- ""

write.csv(output2.mut,'Analysis/Clinical/Categorical parameter analysis/TCGA_Mutations_no_freq.csv')
write.csv(TCGA_mut_stat.mat,'Analysis/Clinical/Categorical parameter analysis/TCGA_Mutations_fisher_tests.csv')

#3 group comparisons
TCGA_mut_stat.mat<-matrix(0,3,length(xmut))
rownames(TCGA_mut_stat.mat)<-str_c(mut_ref$group1,mut_ref$group2,sep = ':')[c(1,3,2)]
colnames(TCGA_mut_stat.mat)<-names(xmut)
for (i in 1:length(xmut)){
  tryCatch({
    TCGA_mut_stat.mat[,i]<-pairwise_fisher_test(xmut[[i]])$p[c(1,3,2)] %>% p.adjust(method='BH')  
  },error=function(e){cat("ERROR:",conditionMessage(e),"\n")})
  
  #eval(as.symbol(paste('xcat',num_list[i],'_res',sep='')))$p.adj
}

output_mut.tcga$variable <- as.character(output_mut.tcga$variable)
#make 'empty lines' 
empties.mut <- data.frame(variable=unique(output_mut.tcga$variable), stringsAsFactors=F)
empties.mut[,colnames(output_mut.tcga)[-1]] <- ""
output2.mut <- rbind(empties.mut,output_mut.tcga)
output2.mut <- output2.mut[order(output2.mut$variable,output2.mut$value),]
output2.mut$variable[output2.mut$value!=""] <- ""

#write.csv(output2.mut,'Analysis/Clinical/Categorical parameter analysis/TCGA_Mutations_no_freq.csv')
write.csv(TCGA_mut_stat.mat,'Analysis/Clinical/Categorical parameter analysis/TCGA_Mutations_fisher_tests_3comb.csv')


#save.image('Analysis/Clinical/Categorical parameter analysis/TP53mutl_clinical_mut.RData')
load('Analysis/Clinical/Categorical parameter analysis/TP53mutl_clinical_mut.RData')

#####################################################################################
#other categorical results - do it later

# BEAT_LASSO_final_res_cat<-BEAT_clinical_final.old_403[,c(4,6,9,11,12,14,19,21:58)]
# BEAT_LASSO_final_res_cat<-right_join(BEAT_LASSO_final_res_cat,BEAT_clinical_final[,c(1,2,113)],by='LabId')
# BEAT_LASSO_final_res_cat[,c(8:44)]<-data.frame(lapply(BEAT_LASSO_final_res_cat[,c(8:44)],as.character))
# BEAT_LASSO_final_cat<-BEAT_LASSO_final_res_cat[,-c(8:44)]
# 
# TCGA_clinical.final.cat<-TCGA_clinical.final.n[,c(1,15,26,27,36,2014,2068,2088,2101,2102)]
# TCGA_clinical.final.cat.mut<-right_join(TCGA_clinical.final.cat,
#                                         rownames_to_column(df01,var='patient_RNA_id'),by='patient_RNA_id')
# TCGA_clinical.final.mut<-TCGA_clinical.final.cat.mut[,-c(2:8,10)]
# TCGA_clinical.final.cat<-TCGA_clinical.final.cat.mut[,c(2:8,10)]

BEAT_LASSO_final_cat
TCGA_clinical.final.cat

###################################
#ELN and Other cytogenetics(LB annotated)
library(readxl)
library(stringr)
library(dplyr)

BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')

#analyze cytogenetics that Dr. Linda Baughn and ZS re-annotation
BEAT_final_cyto.fin01<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/BEAT_final_cyto.fin01.rds')
BEAT_final_cyto.fin02<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/BEAT_final_cyto.fin02.rds')

BEAT_final_cyto.fin01<-BEAT_final_cyto.fin01[which(BEAT_final_cyto.fin01$LabId%in% BEAT_clinical_final$LabId),]
BEAT_final_cyto.fin02<-BEAT_final_cyto.fin02[which(BEAT_final_cyto.fin02$LabId%in%BEAT_clinical_final$LabId),]

TCGA_final_cyto.fin01<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/TCGA_final_cyto.fin01.rds')
TCGA_final_cyto.fin02<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/TCGA_final_cyto.fin02.rds')

# c("cls.GLM.x","ageAtDiagnosis.x","isRelapse.x","isDenovo","FAB/Blast.Morphology.x","overallSurvival.y","t(8;21)",
#   "t(16;16) or inv(16)","t(6;9)","t(v;11) (except t(9;11))","t(9;22)","inv(3) or t(3;3)", "-5","17p alteration/deletion","t(8;21) or inv(16) or t(16;16)", "tri 8",                         
# "tri 21","tri 8 and tri 21","FLT3-ITD","NPM1.y","WES_ASXL1","WES_CEBPA", "WES_RUNX1")   
colnames(BEAT_final_cyto.fin01)[5:dim(BEAT_final_cyto.fin01)[2]]
library(reshape2)
shp01<-melt(BEAT_final_cyto.fin01,measure.vars = colnames(BEAT_final_cyto.fin01)[5:dim(BEAT_final_cyto.fin01)[2]])#"TP53_mut_stat"
cat_res01<-shp01 %>% group_by(cls.GLM.x,variable,value) %>% dplyr::summarize(n=n(),na.rm=TRUE) %>% mutate(freq=n/sum(n,na.rm=T))

#make an 'export' variable
cat_res01$export <- with(cat_res01, sprintf("%i (%.1f%%)", n, freq*100))
#reshape again
output01 <- dcast(variable+value~cls.GLM.x, value.var="export", data=cat_res01, fill="0") #use drop=F to prevent silent missings 
#'silent missings'
output01$variable <- as.character(output01$variable)
#make 'empty lines' 
empties <- data.frame(variable=unique(output01$variable), stringsAsFactors=F)
empties[,colnames(output01)[-1]] <- ""

#bind them together
output01.n <- rbind(empties,output01)
output01.n <- output01.n[order(output01.n$variable,output01.n$value),]

#optional: 'remove' variable if value present
output01.n$variable[output01.n$value!=""] <- ""

write.csv(output01.n,'Analysis/Clinical/Categorical parameter analysis/BEAT_AML_BL_ZS_cytogenetics_v01.csv')

cat_res_t01<-shp01 %>% group_by(cls.GLM.x,variable,value) %>% dplyr::summarize(n=n(),na.rm=F) 
#cat_res_t$export <- with(cat_res_t, sprintf("%i (%.1f%%)", n, freq*100))
#use NA then 0 - for the testing purpose
output_t01 <- dcast(variable+value~cls.GLM.x, value.var="n", data=cat_res_t01, fill=0) #use drop=F to prevent silent missings 
#-c(13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,60,63,66,69)
output_t01.n<-output_t01[-c(13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,63,66,69),]
Cat_group<-c('TP53_MUT','TP53_WT','TP53MUT_like')
Cat_var<-unique(output_t01.n$variable)
var_num<-seq(1:length(Cat_var))
xcat<-vector(mode='list',length=length(Cat_var))
names(xcat)<-Cat_var
for(i in 1:length(Cat_var)){
  ix<-which(output_t01.n$variable%in% Cat_var[i])
  xcat.sub<-as.table(as.matrix(output_t01.n[ix,3:5]))
  dimnames(xcat.sub)<-list(output_t01.n[ix,2],
                           Group=Cat_group)
  xcat[[i]]<-xcat.sub
}
xcat[[1]]
xcat[[2]]

#fisher:1,4,6:42
#chi-square:2,3,5,43
library(rstatix)
xcat01_res<-pairwise_fisher_test(xcat[[1]],p.adjust.method = 'BH')
xcat02_res<-pairwise_fisher_test(xcat[[2]])
xcat03_res<-pairwise_fisher_test(xcat[[3]])
xcat04_res<-pairwise_fisher_test(xcat[[4]])
xcat05_res<-pairwise_fisher_test(xcat[[5]])#-5
#pairwise_fisher_test(xcat[[5]][c(2,1),])
xcat06_res<-pairwise_fisher_test(xcat[[6]])#-7:TP53MUT
xcat07_res<-pairwise_fisher_test(xcat[[7]])#comp or mono
xcat08_res<-pairwise_fisher_test(xcat[[8]])#17p alt
xcat09_res<-pairwise_fisher_test(xcat[[9]])
xcat10_res<-pairwise_fisher_test(xcat[[10]])
xcat11_res<-pairwise_fisher_test(xcat[[11]])#FLT3-ITD: TP53MUT
xcat12_res<-pairwise_fisher_test(xcat[[12]])#NPM1
#pairwise_fisher_test(xcat[[12]][c(2,1),])
xcat13_res<-pairwise_fisher_test(xcat[[13]])
xcat14_res<-pairwise_fisher_test(xcat[[14]])
xcat15_res<-pairwise_fisher_test(xcat[[15]])
xcat16_res<-pairwise_fisher_test(xcat[[16]])
xcat17_res<-pairwise_fisher_test(xcat[[17]])
xcat18_res<-pairwise_fisher_test(xcat[[18]])
xcat19_res<-pairwise_fisher_test(xcat[[19]])
xcat20_res<-pairwise_fisher_test(xcat[[20]])
xcat21_res<-pairwise_fisher_test(xcat[[21]])
xcat22_res<-pairwise_fisher_test(xcat[[22]])
xcat23_res<-pairwise_fisher_test(xcat[[23]])
xcat24_res<-pairwise_fisher_test(xcat[[24]])
xcat25_res<-pairwise_fisher_test(xcat[[25]])


#BEAT_cyto_stat.mat<-matrix(0,3,length(xcat))
BEAT_cyto_stat.mat<-matrix(0,2,length(xcat))
rownames(BEAT_cyto_stat.mat)<-str_c(xcat01_res$group1,xcat01_res$group2,sep = ':')[c(1,3)]
colnames(BEAT_cyto_stat.mat)<-names(xcat)
num_list<-c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',
            '16','17','18','19','20','21','22','23','24','25')
for (i in 1:length(xcat)){
  BEAT_cyto_stat.mat[,i]<-eval(as.symbol(paste('xcat',num_list[i],'_res',sep='')))$p[c(1,3)] %>%p.adjust(method = 'BH')
}
write.csv(BEAT_cyto_stat.mat,'Analysis/Clinical/Categorical parameter analysis/BEAT_AML_BL_ZS_cytogenetics_fisher_v01.csv')
#3 comparisons
BEAT_cyto_stat.mat<-matrix(0,3,length(xcat))
#BEAT_cyto_stat.mat<-matrix(0,2,length(xcat))
rownames(BEAT_cyto_stat.mat)<-str_c(xcat01_res$group1,xcat01_res$group2,sep = ':')[c(1,3,2)]
colnames(BEAT_cyto_stat.mat)<-names(xcat)
num_list<-c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',
            '16','17','18','19','20','21','22','23','24','25')
for (i in 1:length(xcat)){
  BEAT_cyto_stat.mat[,i]<-eval(as.symbol(paste('xcat',num_list[i],'_res',sep='')))$p[c(1,3,2)] %>%p.adjust(method = 'BH')
}
write.csv(BEAT_cyto_stat.mat,'Analysis/Clinical/Categorical parameter analysis/BEAT_AML_BL_ZS_cytogenetics_fisher_v01_3comb.csv')

#version02
shp02<-melt(BEAT_final_cyto.fin02,measure.vars = colnames(BEAT_final_cyto.fin02)[10:dim(BEAT_final_cyto.fin02)[2]])#"TP53_mut_stat"
cat_res02<-shp02 %>% group_by(cls.GLM.x,variable,value) %>% dplyr::summarize(n=n(),na.rm=TRUE) %>% mutate(freq=n/sum(n,na.rm=T))

#make an 'export' variable
cat_res02$export <- with(cat_res02, sprintf("%i (%.1f%%)", n, freq*100))
#reshape again
output02 <- dcast(variable+value~cls.GLM.x, value.var="export", data=cat_res02, fill="0") #use drop=F to prevent silent missings 
#'silent missings'
output02$variable <- as.character(output02$variable)
#make 'empty lines' 
empties02 <- data.frame(variable=unique(output02$variable), stringsAsFactors=F)
empties02[,colnames(output02)[-1]] <- ""

#bind them together
output02.n <- rbind(empties02,output02)
output02.n <- output02.n[order(output02.n$variable,output02.n$value),]

#optional: 'remove' variable if value present
output02.n$variable[output02.n$value!=""] <- ""

write.csv(output02.n,'Analysis/Clinical/Categorical parameter analysis/BEAT_AML_BL_ZS_cytogenetics_v02.csv')

cat_res_t02<-shp02 %>% group_by(cls.GLM.x,variable,value) %>% dplyr::summarize(n=n(),na.rm=F) 
#cat_res_t$export <- with(cat_res_t, sprintf("%i (%.1f%%)", n, freq*100))
#use NA then 0 - for the testing purpose
output_t02 <- dcast(variable+value~cls.GLM.x, value.var="n", data=cat_res_t02, fill=0) #use drop=F to prevent silent missings 
#-c(13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,60,63,66,69)
output_t02.n<-output_t02[-c(3,6,9,12,15,18,21,24,29),]
Cat_group<-c('TP53_MUT','TP53_WT','TP53MUT_like')
Cat_var02<-unique(output_t02.n$variable)
var_num02<-seq(1:length(Cat_var02))
xcat02<-vector(mode='list',length=length(Cat_var02))
names(xcat02)<-Cat_var02
for(i in 1:length(Cat_var02)){
  ix<-which(output_t02.n$variable%in% Cat_var02[i])
  xcat.sub<-as.table(as.matrix(output_t02.n[ix,3:5]))
  dimnames(xcat.sub)<-list(output_t02.n[ix,2],
                           Group=Cat_group)
  xcat02[[i]]<-xcat.sub
}
xcat02[[1]]
xcat02[[2]]

#fisher:1,4,6:42
#chi-square:2,3,5,43
library(rstatix)
xcat02_res01<-pairwise_fisher_test(xcat02[[1]])
xcat02_res02<-pairwise_fisher_test(xcat02[[2]])
xcat02_res03<-pairwise_fisher_test(xcat02[[3]])
xcat02_res04<-pairwise_fisher_test(xcat02[[4]])
xcat02_res05<-pairwise_fisher_test(xcat02[[5]])#-5
#pairwise_fisher_test(xcat[[5]][c(2,1),])
xcat02_res06<-pairwise_fisher_test(xcat02[[6]])#-7:TP53MUT
xcat02_res07<-pairwise_fisher_test(xcat02[[7]])#comp or mono
xcat02_res08<-pairwise_fisher_test(xcat02[[8]])#17p alt
xcat02_res09<-pairwise_fisher_test(xcat02[[9]])
xcat02_res10<-pairwise_fisher_test(xcat02[[10]])


#BEAT_cyto_stat02.mat<-matrix(0,3,length(xcat02))
BEAT_cyto_stat02.mat<-matrix(0,2,length(xcat02))
rownames(BEAT_cyto_stat02.mat)<-str_c(xcat02_res01$group1,xcat02_res01$group2,sep = ':')[c(1,3)]
colnames(BEAT_cyto_stat02.mat)<-names(xcat02)
num_list<-c('01','02','03','04','05','06','07','08','09','10')
for (i in 1:length(xcat02)){
  BEAT_cyto_stat02.mat[,i]<-eval(as.symbol(paste('xcat02','_res',num_list[i],sep='')))$p[c(1,3)] %>%p.adjust(method = 'BH')
}
write.csv(BEAT_cyto_stat02.mat,'Analysis/Clinical/Categorical parameter analysis/BEAT_AML_BL_ZS_cytogenetics_fisher_v02.csv')

#3combinations

BEAT_cyto_stat02.mat<-matrix(0,3,length(xcat02))
#BEAT_cyto_stat02.mat<-matrix(0,2,length(xcat02))
rownames(BEAT_cyto_stat02.mat)<-str_c(xcat02_res01$group1,xcat02_res01$group2,sep = ':')[c(1,3,2)]
colnames(BEAT_cyto_stat02.mat)<-names(xcat02)
num_list<-c('01','02','03','04','05','06','07','08','09','10')
for (i in 1:length(xcat02)){
  BEAT_cyto_stat02.mat[,i]<-eval(as.symbol(paste('xcat02','_res',num_list[i],sep='')))$p[c(1,3,2)] %>%p.adjust(method = 'BH')
}
write.csv(BEAT_cyto_stat02.mat,'Analysis/Clinical/Categorical parameter analysis/BEAT_AML_BL_ZS_cytogenetics_fisher_v02_3comb.csv')

###########################################################################################
#TCGA_cyto


# saveRDS(TCGA_final_cyto.fin01,'class ids/TCGA_final_cyto.fin01.rds')
# saveRDS(TCGA_final_cyto.fin02,'class ids/TCGA_final_cyto.fin02.rds')

# c("cls.GLM.x","ageAtDiagnosis.x","isRelapse.x","isDenovo","FAB/Blast.Morphology.x","overallSurvival.y","t(8;21)",
#   "t(16;16) or inv(16)","t(6;9)","t(v;11) (except t(9;11))","t(9;22)","inv(3) or t(3;3)", "-5","17p alteration/deletion","t(8;21) or inv(16) or t(16;16)", "tri 8",                         
# "tri 21","tri 8 and tri 21","FLT3-ITD","NPM1.y","WES_ASXL1","WES_CEBPA", "WES_RUNX1")   
colnames(TCGA_final_cyto.fin01)[5:dim(TCGA_final_cyto.fin01)[2]]
library(reshape2)
shp.tcga01<-melt(TCGA_final_cyto.fin01,measure.vars = colnames(TCGA_final_cyto.fin01)[5:dim(TCGA_final_cyto.fin01)[2]])#"TP53_mut_stat"
cat_res.tcga01<-shp.tcga01 %>% group_by(cls.GLM.x,variable,value) %>% dplyr::summarize(n=n(),na.rm=TRUE) %>% mutate(freq=n/sum(n,na.rm=T))

#make an 'export' variable
cat_res.tcga01$export <- with(cat_res.tcga01, sprintf("%i (%.1f%%)", n, freq*100))
#reshape again
output.tcga01 <- dcast(variable+value~cls.GLM.x, value.var="export", data=cat_res.tcga01, fill="0") #use drop=F to prevent silent missings 
#'silent missings'
output.tcga01$variable <- as.character(output.tcga01$variable)
#make 'empty lines' 
empties <- data.frame(variable=unique(output.tcga01$variable), stringsAsFactors=F)
empties[,colnames(output.tcga01)[-1]] <- ""

#bind them together
output.tcga01.n <- rbind(empties,output.tcga01)
output.tcga01.n <- output.tcga01.n[order(output.tcga01.n$variable,output.tcga01.n$value),]

#optional: 'remove' variable if value present
output.tcga01.n$variable[output.tcga01.n$value!=""] <- ""

write.csv(output.tcga01.n,'Analysis/Clinical/Categorical parameter analysis/TCGA_AML_BL_ZS_cytogenetics_v01.csv')

cat_res_t.tcga01<-shp.tcga01 %>% group_by(cls.GLM.x,variable,value) %>% dplyr::summarize(n=n(),na.rm=F) 
#cat_res_t$export <- with(cat_res_t, sprintf("%i (%.1f%%)", n, freq*100))
#use NA then 0 - for the testing purpose
output_t.tcga01 <- dcast(variable+value~cls.GLM.x, value.var="n", data=cat_res_t.tcga01, fill=0) #use drop=F to prevent silent missings 
#-c(13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,60,63,66,69)
output_t.tcga01<-output_t.tcga01[-c(11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,61,64),]
Cat_group<-c('TP53_MUT','TP53_WT','TP53MUT_like')
Cat_var<-unique(output_t.tcga01$variable)


var_num<-seq(1:length(Cat_var))
xcat01.tcga<-vector(mode='list',length=length(Cat_var))
names(xcat01.tcga)<-Cat_var
for(i in 1:length(Cat_var)){
  ix<-which(output_t.tcga01$variable%in% Cat_var[i])
  xcat.sub<-as.table(as.matrix(output_t.tcga01[ix,3:5]))
  dimnames(xcat.sub)<-list(output_t.tcga01[ix,2],
                           Group=Cat_group)
  xcat01.tcga[[i]]<-xcat.sub
}
xcat01.tcga[[1]]
xcat01.tcga[[2]]

#fisher:1,4,6:42
#chi-square:2,3,5,43
library(rstatix)
xcat01.tcga01_res<-pairwise_fisher_test(xcat01.tcga[[1]])
xcat01.tcga02_res<-pairwise_fisher_test(xcat01.tcga[[2]])
xcat01.tcga03_res<-pairwise_fisher_test(xcat01.tcga[[3]])
xcat01.tcga04_res<-pairwise_fisher_test(xcat01.tcga[[4]])
xcat01.tcga05_res<-pairwise_fisher_test(xcat01.tcga[[5]])#-5
#pairwise_fisher_test(xcat[[5]][c(2,1),])
xcat01.tcga06_res<-pairwise_fisher_test(xcat01.tcga[[6]])#-7:TP53MUT
xcat01.tcga07_res<-pairwise_fisher_test(xcat01.tcga[[7]])#comp or mono
xcat01.tcga08_res<-pairwise_fisher_test(xcat01.tcga[[8]])#17p alt
xcat01.tcga09_res<-pairwise_fisher_test(xcat01.tcga[[9]])
xcat01.tcga10_res<-pairwise_fisher_test(xcat01.tcga[[10]])
xcat01.tcga11_res<-pairwise_fisher_test(xcat01.tcga[[11]])#FLT3-ITD: TP53MUT
xcat01.tcga12_res<-pairwise_fisher_test(xcat01.tcga[[12]])#NPM1
#pairwise_fisher_test(xcat[[12]][c(2,1),])
xcat01.tcga13_res<-pairwise_fisher_test(xcat01.tcga[[13]])
xcat01.tcga14_res<-pairwise_fisher_test(xcat01.tcga[[14]])
xcat01.tcga15_res<-pairwise_fisher_test(xcat01.tcga[[15]])
xcat01.tcga16_res<-pairwise_fisher_test(xcat01.tcga[[16]])
xcat01.tcga17_res<-pairwise_fisher_test(xcat01.tcga[[17]])
xcat01.tcga18_res<-pairwise_fisher_test(xcat01.tcga[[18]])
xcat01.tcga19_res<-pairwise_fisher_test(xcat01.tcga[[19]])
xcat01.tcga20_res<-pairwise_fisher_test(xcat01.tcga[[20]])
xcat01.tcga21_res<-pairwise_fisher_test(xcat01.tcga[[21]])
xcat01.tcga22_res<-pairwise_fisher_test(xcat01.tcga[[22]])
xcat01.tcga23_res<-pairwise_fisher_test(xcat01.tcga[[23]])


#TCGA_cyto_stat01.mat<-matrix(0,3,length(xcat01.tcga))
TCGA_cyto_stat01.mat<-matrix(0,2,length(xcat01.tcga))
rownames(TCGA_cyto_stat01.mat)<-str_c(xcat01.tcga01_res$group1,xcat01.tcga01_res$group2,sep = ':')[c(1,3)]
colnames(TCGA_cyto_stat01.mat)<-names(xcat01.tcga)
num_list<-c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',
            '16','17','18','19','20','21','22','23')
for (i in 1:length(xcat01.tcga)){
  TCGA_cyto_stat01.mat[,i]<-eval(as.symbol(paste('xcat01.tcga',num_list[i],'_res',sep='')))$p[c(1,3)]%>% p.adjust(method = 'BH')
}
write.csv(TCGA_cyto_stat01.mat,'Analysis/Clinical/Categorical parameter analysis/TCGA_AML_BL_ZS_cytogenetics_fisher_v01.csv')

#3group comparisons
TCGA_cyto_stat01.mat<-matrix(0,3,length(xcat01.tcga))
#TCGA_cyto_stat01.mat<-matrix(0,2,length(xcat01.tcga))
rownames(TCGA_cyto_stat01.mat)<-str_c(xcat01.tcga01_res$group1,xcat01.tcga01_res$group2,sep = ':')[c(1,3,2)]
colnames(TCGA_cyto_stat01.mat)<-names(xcat01.tcga)
num_list<-c('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15',
            '16','17','18','19','20','21','22','23')
for (i in 1:length(xcat01.tcga)){
  TCGA_cyto_stat01.mat[,i]<-eval(as.symbol(paste('xcat01.tcga',num_list[i],'_res',sep='')))$p[c(1,3,2)]%>% p.adjust(method = 'BH')
}
write.csv(TCGA_cyto_stat01.mat,'Analysis/Clinical/Categorical parameter analysis/TCGA_AML_BL_ZS_cytogenetics_fisher_v01_3comb.csv')


#version02
shp.tcga02<-melt(TCGA_final_cyto.fin02,measure.vars = colnames(TCGA_final_cyto.fin02)[9:dim(TCGA_final_cyto.fin02)[2]])#"TP53_mut_stat"
cat_res.tcga02<-shp.tcga02 %>% group_by(cls.GLM.x,variable,value) %>% dplyr::summarize(n=n(),na.rm=TRUE) %>% mutate(freq=n/sum(n,na.rm=T))

#make an 'export' variable
cat_res.tcga02$export <- with(cat_res.tcga02, sprintf("%i (%.1f%%)", n, freq*100))
#reshape again
output.tcga02 <- dcast(variable+value~cls.GLM.x, value.var="export", data=cat_res.tcga02, fill="0") #use drop=F to prevent silent missings 
#'silent missings'
output.tcga02$variable <- as.character(output.tcga02$variable)
#make 'empty lines' 
empties02 <- data.frame(variable=unique(output.tcga02$variable), stringsAsFactors=F)
empties02[,colnames(output.tcga02)[-1]] <- ""

#bind them together
output.tcga02.n <- rbind(empties02,output.tcga02)
output.tcga02.n <- output.tcga02.n[order(output.tcga02.n$variable,output.tcga02.n$value),]

#optional: 'remove' variable if value present
output.tcga02.n$variable[output.tcga02.n$value!=""] <- ""

write.csv(output.tcga02.n,'Analysis/Clinical/Categorical parameter analysis//TCGA_AML_BL_ZS_cytogenetics_v02.csv')

cat_res_t.tcga02<-shp.tcga02 %>% group_by(cls.GLM.x,variable,value) %>% dplyr::summarize(n=n(),na.rm=F) 
#cat_res_t$export <- with(cat_res_t, sprintf("%i (%.1f%%)", n, freq*100))
#use NA then 0 - for the testing purpose
output_t.tcga02 <- dcast(variable+value~cls.GLM.x, value.var="n", data=cat_res_t.tcga02, fill=0) #use drop=F to prevent silent missings 
#-c(13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,60,63,66,69)
output_t.tcga02<-output_t.tcga02[-c(3,6,9,12,15,18,21,24,27,30),]
Cat_group<-c('TP53_MUT','TP53_WT','TP53MUT_like')
Cat_var02<-unique(output_t.tcga02$variable)
var_num02<-seq(1:length(Cat_var02))
xcat02.tcga<-vector(mode='list',length=length(Cat_var02))
names(xcat02.tcga)<-Cat_var02
for(i in 1:length(Cat_var02)){
  ix<-which(output_t.tcga02$variable%in% Cat_var02[i])
  xcat.sub<-as.table(as.matrix(output_t.tcga02[ix,3:5]))
  dimnames(xcat.sub)<-list(output_t.tcga02[ix,2],
                           Group=Cat_group)
  xcat02.tcga[[i]]<-xcat.sub
}
xcat02.tcga[[1]]
xcat02.tcga[[2]]

#fisher:1,4,6:42
#chi-square:2,3,5,43
library(rstatix)
xcat02.tcga01<-pairwise_fisher_test(xcat02.tcga[[1]])
xcat02.tcga02<-pairwise_fisher_test(xcat02.tcga[[2]])
xcat02.tcga03<-pairwise_fisher_test(xcat02.tcga[[3]])
xcat02.tcga04<-pairwise_fisher_test(xcat02.tcga[[4]])
xcat02.tcga05<-pairwise_fisher_test(xcat02.tcga[[5]])#-5
#pairwise_fisher_test(xcat[[5]][c(2,1),])
xcat02.tcga06<-pairwise_fisher_test(xcat02.tcga[[6]])#-7:TP53MUT
xcat02.tcga07<-pairwise_fisher_test(xcat02.tcga[[7]])#comp or mono
xcat02.tcga08<-pairwise_fisher_test(xcat02.tcga[[8]])#17p alt
xcat02.tcga09<-pairwise_fisher_test(xcat02.tcga[[9]])
xcat02.tcga10<-pairwise_fisher_test(xcat02.tcga[[10]])
xcat02.tcga11<-pairwise_fisher_test(xcat02.tcga[[11]])


#TCGA_cyto_stat02.mat<-matrix(0,3,length(xcat02.tcga))
TCGA_cyto_stat02.mat<-matrix(0,2,length(xcat02.tcga))
rownames(TCGA_cyto_stat02.mat)<-str_c(xcat02.tcga01$group1,xcat02.tcga01$group2,sep = ':')[c(1,3)]
colnames(TCGA_cyto_stat02.mat)<-names(xcat02.tcga)
num_list<-c('01','02','03','04','05','06','07','08','09','10','11')
for (i in 1:length(xcat02.tcga)){
  TCGA_cyto_stat02.mat[,i]<-eval(as.symbol(paste('xcat02.tcga',num_list[i],sep='')))$p[c(1,3)]%>%p.adjust(method = 'BH')
}
write.csv(TCGA_cyto_stat02.mat,'Analysis/Clinical/Categorical parameter analysis/TCGA_AML_BL_ZS_cytogenetics_fisher_v02.csv')
save.image('Analysis/Clinical/Categorical parameter analysis/TP53mutl_cyto.RData')
load('Analysis/Clinical/Categorical parameter analysis/TP53mutl_cyto.RData')

#3 groups combinations
TCGA_cyto_stat02.mat<-matrix(0,3,length(xcat02.tcga))
#TCGA_cyto_stat02.mat<-matrix(0,2,length(xcat02.tcga))
rownames(TCGA_cyto_stat02.mat)<-str_c(xcat02.tcga01$group1,xcat02.tcga01$group2,sep = ':')[c(1,3,2)]
colnames(TCGA_cyto_stat02.mat)<-names(xcat02.tcga)
num_list<-c('01','02','03','04','05','06','07','08','09','10','11')
for (i in 1:length(xcat02.tcga)){
  TCGA_cyto_stat02.mat[,i]<-eval(as.symbol(paste('xcat02.tcga',num_list[i],sep='')))$p[c(1,3,2)]%>%p.adjust(method = 'BH')
}
write.csv(TCGA_cyto_stat02.mat,'Analysis/Clinical/Categorical parameter analysis/TCGA_AML_BL_ZS_cytogenetics_fisher_v02_3comb.csv')



#Visualization: binarized map(mutation and cytogenetics)

binary_sorting<-function(mat){
  #mdf.n<-Final.DEG.mat.n.t[,which(Final.DEG.mat.n.anno$cls=='TP53MUT_like')]
  mdf.n<-mat
  mdf.n = mdf.n[order(mdf.n[, ncol(mdf.n)], decreasing = TRUE), ]
  #colnames(mdf.n) = gsub(pattern = "^X", replacement = "", colnames(mdf.n))
  nMut = mdf.n[, ncol(mdf.n)]
  
  #mdf.n = mdf.n[, -ncol(mdf.n)]
  
  mdf.n.temp.copy = mdf.n #temp copy of original unsorted numeric coded matrix
  
  mdf.n[mdf.n != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
  tmdf.n = t(mdf.n) #transposematrix
  mdf.n = t(tmdf.n[do.call(order, c(as.list(as.data.frame(tmdf.n)), decreasing = TRUE)), ]) #sort
  
  mdf.n.temp.copy = mdf.n.temp.copy[rownames(mdf.n),] #organise original matrix into sorted matrix
  mdf.n.temp.copy = mdf.n.temp.copy[,colnames(mdf.n)]
  mdf.n = mdf.n.temp.copy
  return(mdf.n)
}

#####################################################################
BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/TCGA_clinical_final_class.rds')
# BEAT_clinical_final.old<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/BEAT_clinical_final_class.rds')
# BEAT_clinical_final.old_403<-BEAT_clinical_final.old[match(BEAT_clinical_final$LabId,BEAT_clinical_final.old$LabId),]

BEAT_clinical_cyto_fin_v01<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/BEAT_final_cyto.fin01.rds')
BEAT_clinical_cyto_fin_v02<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/BEAT_final_cyto.fin02.rds')

BEAT_clinical_cyto_fin_v01<-BEAT_clinical_cyto_fin_v01[which(BEAT_clinical_cyto_fin_v01$LabId%in% BEAT_clinical_final$LabId),]
BEAT_clinical_cyto_fin_v02<-BEAT_clinical_cyto_fin_v02[which(BEAT_clinical_cyto_fin_v02$LabId%in%BEAT_clinical_final$LabId),]

TCGA_clinical_cyto_fin_v01<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/TCGA_final_cyto.fin01.rds')
TCGA_clinical_cyto_fin_v02<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/TCGA_final_cyto.fin02.rds')


for (i in 1:dim(BEAT_clinical_cyto_fin_v02)[1]){
  cls<-names(table(BEAT_clinical_cyto_fin_v02$cls.GLM))
  if(BEAT_clinical_cyto_fin_v02$cls.GLM[i]=='TP53_WT'){#TP53wt
    BEAT_clinical_cyto_fin_v02$cls.GLM[i]<-'TP53wt'
  }
  if(BEAT_clinical_cyto_fin_v02$cls.GLM[i]=='TP53_MUT'){#TP53mut
    BEAT_clinical_cyto_fin_v02$cls.GLM[i]<-'TP53mut'
  }
  if(BEAT_clinical_cyto_fin_v02$cls.GLM[i]=='TP53MUT_like'){#TP53mut-like
    BEAT_clinical_cyto_fin_v02$cls.GLM[i]<-'TP53mut-like'
  }
}
for (i in 1:dim(TCGA_clinical_cyto_fin_v02)[1]){
  cls<-names(table(TCGA_clinical_cyto_fin_v02$cls.GLM.x))
  if(TCGA_clinical_cyto_fin_v02$cls.GLM.x[i]=='TP53_WT'){#TP53wt
    TCGA_clinical_cyto_fin_v02$cls.GLM.x[i]<-'TP53wt'
  }
  if(TCGA_clinical_cyto_fin_v02$cls.GLM.x[i]=='TP53_MUT'){#TP53mut
    TCGA_clinical_cyto_fin_v02$cls.GLM.x[i]<-'TP53mut'
  }
  if(TCGA_clinical_cyto_fin_v02$cls.GLM.x[i]=='TP53MUT_like'){#TP53mut-like
    TCGA_clinical_cyto_fin_v02$cls.GLM.x[i]<-'TP53mut-like'
  }
}






library(ComplexHeatmap)
# Comb_clinical_cyto_fin_v02<-data.frame('ID'=c(BEAT_clinical_cyto_fin_v02[,1],TCGA_clinical_cyto_fin_v02[,1]),
#                                        'cls'=c(BEAT_clinical_cyto_fin_v02[,2],TCGA_clinical_cyto_fin_v02[,2]),
#                                        'FLT3-ITD'=c(BEAT_clinical_cyto_fin_v02$`FLT3-ITD`,rep(NA,dim(TCGA_clinical_cyto_fin_v02)[1])),###
#                                        'NPM1'=c(BEAT_clinical_cyto_fin_v02$NPM1,TCGA_clinical_cyto_fin_v02$NPM1),
#                                        'ASXL1'=c(BEAT_clinical_cyto_fin_v02$WES_ASXL1,TCGA_clinical_cyto_fin_v02$ASXL1),
#                                        'CEBPA'=c(BEAT_clinical_cyto_fin_v02$WES_CEBPA,TCGA_clinical_cyto_fin_v02$CEBPA),
#                                        'RUNX1'=c(BEAT_clinical_cyto_fin_v02$WES_RUNX1,TCGA_clinical_cyto_fin_v02$RUNX1),
#                                        'Favorable'=c(BEAT_clinical_cyto_fin_v02$Favorable,TCGA_clinical_cyto_fin_v02$Favorable),
#                                        'del5'=c(BEAT_clinical_cyto_fin_v02$del5,TCGA_clinical_cyto_fin_v02$del5),
#                                        '17pAlt'=c(BEAT_clinical_cyto_fin_v02$`17pAlt`,TCGA_clinical_cyto_fin_v02$`17pAlt`),
#                                        'del7'=c(BEAT_clinical_cyto_fin_v02$del7,TCGA_clinical_cyto_fin_v02$del7),
#                                        'complex'=c(BEAT_clinical_cyto_fin_v02$complex,TCGA_clinical_cyto_fin_v02$complex),
#                                        'Intermediate'=c(BEAT_clinical_cyto_fin_v02$Intermediate,TCGA_clinical_cyto_fin_v02$Intermediate),
#                                        'OtherAdverse'=c(BEAT_clinical_cyto_fin_v02$OtherAdverse,TCGA_clinical_cyto_fin_v02$OtherAdverse),
#                                        'monosomal'=c(BEAT_clinical_cyto_fin_v02$monosomal,TCGA_clinical_cyto_fin_v02$monosomal))
# rownames(Comb_clinical_cyto_fin_v02)<-Comb_clinical_cyto_fin_v02$ID
# Comb_clinical_cyto_fin_v02.anno<-Comb_clinical_cyto_fin_v02[,c(1,2)]

BEAT_clinical_cyto_fin_v02.df<-data.frame('ID'=c(BEAT_clinical_cyto_fin_v02[,1]),
                                       'cls'=c(BEAT_clinical_cyto_fin_v02[,20]),
                                       'FLT3-ITD'=c(BEAT_clinical_cyto_fin_v02$`FLT3-ITD`),###
                                       'NPM1'=c(BEAT_clinical_cyto_fin_v02$NPM1),
                                       'ASXL1'=c(BEAT_clinical_cyto_fin_v02$WES_ASXL1),
                                       'CEBPA'=c(BEAT_clinical_cyto_fin_v02$WES_CEBPA),
                                       'RUNX1'=c(BEAT_clinical_cyto_fin_v02$WES_RUNX1),
                                       'Favorable'=c(BEAT_clinical_cyto_fin_v02$Favorable),
                                       'del5'=c(BEAT_clinical_cyto_fin_v02$del5),
                                       '17pAlt'=c(BEAT_clinical_cyto_fin_v02$`17pAlt`),
                                       'del7'=c(BEAT_clinical_cyto_fin_v02$del7),
                                       'complex'=c(BEAT_clinical_cyto_fin_v02$complex),
                                       'Intermediate'=c(BEAT_clinical_cyto_fin_v02$Intermediate),
                                       'OtherAdverse'=c(BEAT_clinical_cyto_fin_v02$OtherAdverse),
                                       'monosomal'=c(BEAT_clinical_cyto_fin_v02$monosomal))
rownames(BEAT_clinical_cyto_fin_v02.df)<-BEAT_clinical_cyto_fin_v02.df$ID
BEAT_clinical_cyto_fin_v02.anno<-BEAT_clinical_cyto_fin_v02.df[,c(1,2)]

TCGA_clinical_cyto_fin_v02.df<-data.frame('ID'=c(TCGA_clinical_cyto_fin_v02[,1]),
                                       'cls'=c(TCGA_clinical_cyto_fin_v02[,2]),
                                       'FLT3-ITD'=c(rep(NA,dim(TCGA_clinical_cyto_fin_v02)[1])),###
                                       'NPM1'=c(TCGA_clinical_cyto_fin_v02$NPM1),
                                       'ASXL1'=c(TCGA_clinical_cyto_fin_v02$ASXL1),
                                       'CEBPA'=c(TCGA_clinical_cyto_fin_v02$CEBPA),
                                       'RUNX1'=c(TCGA_clinical_cyto_fin_v02$RUNX1),
                                       'Favorable'=c(TCGA_clinical_cyto_fin_v02$Favorable),
                                       'del5'=c(TCGA_clinical_cyto_fin_v02$del5),
                                       '17pAlt'=c(TCGA_clinical_cyto_fin_v02$`17pAlt`),
                                       'del7'=c(TCGA_clinical_cyto_fin_v02$del7),
                                       'complex'=c(TCGA_clinical_cyto_fin_v02$complex),
                                       'Intermediate'=c(TCGA_clinical_cyto_fin_v02$Intermediate),
                                       'OtherAdverse'=c(TCGA_clinical_cyto_fin_v02$OtherAdverse),
                                       'monosomal'=c(TCGA_clinical_cyto_fin_v02$monosomal))
rownames(TCGA_clinical_cyto_fin_v02.df)<-TCGA_clinical_cyto_fin_v02.df$ID
TCGA_clinical_cyto_fin_v02.anno<-TCGA_clinical_cyto_fin_v02.df[,c(1,2)]

# comb.ref<-as.data.frame(Final.DEG.mat[,1])
# colnames(comb.ref)<-'ID'
# comb.ref.n<-right_join(comb.ref,Final.DEG.mat,by='ID')
# comb.ref.n<-right_join(comb.ref.n,Comb_clinical_cyto_fin_v02,by='ID')
#Final.DEG.mat.n<-Comb_clinical_cyto_fin_v02#change this all the time the DEG chagned
#rownames(Final.DEG.mat.n)<-Comb_clinical_cyto_fin_v02$ID
#Final.DEG.mat.n.t<-as.data.frame(t(Final.DEG.mat.n[,3:dim(Final.DEG.mat.n)[2]]))
BEAT.anno<-data.frame('cls'=BEAT_clinical_cyto_fin_v02.df$cls,
                                 row.names = BEAT_clinical_cyto_fin_v02.df$ID)
BEAT_HM_anno<-HeatmapAnnotation(df=BEAT.anno,col = list(cls=c('TP53mut'='red','TP53wt'='blue','TP53mut-like'='purple')))

TCGA.anno<-data.frame('cls'=TCGA_clinical_cyto_fin_v02.df$cls,
                      row.names = TCGA_clinical_cyto_fin_v02.df$ID)
TCGA_HM_anno<-HeatmapAnnotation(df=TCGA.anno,col = list(cls=c('TP53mut'='red','TP53wt'='blue','TP53mut-like'='purple')))


# mtl_row_sort_ix<-order(rowSums(Final.DEG.mat.n.t[,Final.DEG.mat.n.anno$cls=='TP53MUT_like'],na.rm = T),decreasing = T)
# mtl_col_sort_ix<-order(colSums(Final.DEG.mat.n.t[,Final.DEG.mat.n.anno$cls=='TP53MUT_like'],na.rm = T),decreasing = T)
# 
# #wt_row_sort_ix<-order(rowSums(Final.DEG.mat.n.t[,Final.DEG.mat.n.anno$cls=='TP53_WT'],na.rm = T),decreasing = T)
# wt_col_sort_ix<-order(colSums(Final.DEG.mat.n.t[,Final.DEG.mat.n.anno$cls=='TP53_WT'],na.rm = T),decreasing = T)

#rowSums(Final.DEG.mat.n.t[,Final.DEG.mat.n.anno$cls=='TP53_WT'],na.rm = T)
##################################################################################
##########################################################################
binary_sorting<-function(mat){
  #mdf.n<-Final.DEG.mat.n.t[,which(Final.DEG.mat.n.anno$cls=='TP53MUT_like')]
  mdf.n<-mat
  mdf.n = mdf.n[order(mdf.n[, ncol(mdf.n)], decreasing = TRUE), ]
  #colnames(mdf.n) = gsub(pattern = "^X", replacement = "", colnames(mdf.n))
  nMut = mdf.n[, ncol(mdf.n)]
  
  #mdf.n = mdf.n[, -ncol(mdf.n)]
  
  mdf.n.temp.copy = mdf.n #temp copy of original unsorted numeric coded matrix
  
  mdf.n[mdf.n != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
  tmdf.n = t(mdf.n) #transposematrix
  mdf.n = t(tmdf.n[do.call(order, c(as.list(as.data.frame(tmdf.n)), decreasing = TRUE)), ]) #sort
  
  mdf.n.temp.copy = mdf.n.temp.copy[rownames(mdf.n),] #organise original matrix into sorted matrix
  mdf.n.temp.copy = mdf.n.temp.copy[,colnames(mdf.n)]
  mdf.n = mdf.n.temp.copy
  return(mdf.n)
}
BEAT.anno[(BEAT.anno$cls=='TP53mut-like'),]
TCGA.anno[(TCGA.anno$cls=='TP53mut-like'),]
BEAT_clinical_cyto_fin_v02.df.t<-as.data.frame(t(BEAT_clinical_cyto_fin_v02.df[,3:dim(BEAT_clinical_cyto_fin_v02.df)[2]]))
TCGA_clinical_cyto_fin_v02.df.t<-as.data.frame(t(TCGA_clinical_cyto_fin_v02.df[,3:dim(TCGA_clinical_cyto_fin_v02.df)[2]]))
#sort only using 4 marker genes
beat.mtl.mat.sort<-binary_sorting(BEAT_clinical_cyto_fin_v02.df.t[,which(BEAT.anno$cls=='TP53mut-like')])
beat.mt.mat.sort<-binary_sorting(BEAT_clinical_cyto_fin_v02.df.t[,which(BEAT.anno$cls=='TP53mut')])
beat.wt.mat.sort<-binary_sorting(BEAT_clinical_cyto_fin_v02.df.t[,which(BEAT.anno$cls=='TP53wt')])

tcga.mtl.mat.sort<-binary_sorting(TCGA_clinical_cyto_fin_v02.df.t[,which(TCGA.anno$cls=='TP53mut-like')])
tcga.mt.mat.sort<-binary_sorting(TCGA_clinical_cyto_fin_v02.df.t[,which(TCGA.anno$cls=='TP53mut')])
tcga.wt.mat.sort<-binary_sorting(TCGA_clinical_cyto_fin_v02.df.t[,which(TCGA.anno$cls=='TP53wt')])


wri.beat<-c("FLT3.ITD","NPM1", "ASXL1","CEBPA",
       "RUNX1","del5","X17pAlt","del7","monosomal","Favorable","Intermediate","complex" ,"OtherAdverse") 
wri.tcga<-c("FLT3.ITD","NPM1", "ASXL1","CEBPA",
            "RUNX1","del5","X17pAlt","del7","monosomal","Favorable","Intermediate","complex" ,"OtherAdverse") 
####################################

library(circlize)
col_fun = colorRamp2(c(0,1), c("grey", "black"),transparency = 0.35)



#get Zohar's ELN classification (2022.08.18)

BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')

#View(cbind(BEAT_clinical_final.mtwt.ELN.new$`ELN Risk Group_YL`,BEAT_clinical_final.mtwt.ELN.new$ELN2017))
fav_q_ix<-which(BEAT_clinical_final$`ELN Risk Group_YL`=='?Fav')
BEAT_clinical_final$ELN_Risk_Group_YL_new<-BEAT_clinical_final$`ELN Risk Group_YL`
#BEAT_clinical_final$ELN_Risk_Group_YL_new[fav_q_ix]<-BEAT_clinical_final$ELN2017[fav_q_ix]

for(i in (fav_q_ix)){
  cat<-BEAT_clinical_final$ELN2017[i]
  if(cat=='Favorable'){
    BEAT_clinical_final$ELN_Risk_Group_YL_new[i]<-'Favorable'}
  # }else if(cat=='Intermediate'){
  #   BEAT_clinical_final$ELN_Risk_Group_YL_new[i]<-'Intermediate'
  # }else if(cat=='Adverse'){
  #   BEAT_clinical_final$ELN_Risk_Group_YL_new[i]<-'Adverse'
  # }
}
BEAT_clinical_final$`ELN_Risk_Group_YL_new`[BEAT_clinical_final$`ELN_Risk_Group_YL_new`=='?Fav']<-'unknown'
TCGA_clinical_final$`ELN Risk Group_YL`[TCGA_clinical_final$`ELN Risk Group_YL`=='?Fav']<-'unknown'
table(BEAT_clinical_final$ELN_Risk_Group_YL_new)
library(readxl)
# ELN_ZS<-read_excel("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/p53 mutant Like Project/Version for final/BEAT_TCGA_inte.join_re.ZSeds.8.18.22.xlsx")
# ELN_ZS.n<-ELN_ZS[,c(1,15,17)]
# colnames(ELN_ZS.n)[1]<-'rownames'
# ELN_ZS.n$cohort<-ifelse(substr(ELN_ZS.n$rownames,1,4)=='TCGA','TCGA','BEAT')
# ELN_ZS.n.beat<-ELN_ZS.n[ELN_ZS.n$cohort=='BEAT',]
# tt<-as.data.frame(ELN_ZS.n.beat %>% filter(rownames %in% BEAT_clinical_final$LabId))
# rownames(tt)<-tt$rownames
# View(cbind(tt[BEAT_clinical_final$LabId,],BEAT_clinical_final[,c(1,113)]))
# ELN_ZS.n.tcga<-ELN_ZS.n[ELN_ZS.n$cohort=='TCGA',]
#create reference for final order
BEAT.ref<-data.frame('rownames'=colnames(BEAT_clinical_cyto_fin_v02.df.t[,c(colnames(beat.mt.mat.sort),colnames(beat.mtl.mat.sort),colnames(beat.wt.mat.sort))]))
TCGA.ref<-data.frame('rownames'=colnames(TCGA_clinical_cyto_fin_v02.df.t[,c(colnames(tcga.mt.mat.sort),colnames(tcga.mtl.mat.sort),colnames(tcga.wt.mat.sort))]))
BEAT.ref.anno<-data.frame('cls'=BEAT_clinical_cyto_fin_v02.df$cls,
                                 row.names = BEAT_clinical_cyto_fin_v02.df$ID)%>%rownames_to_column(var='rownames')
TCGA.ref.anno<-data.frame('cls'=TCGA_clinical_cyto_fin_v02.df$cls,
                                 row.names = TCGA_clinical_cyto_fin_v02.df$ID)%>%rownames_to_column(var='rownames')
BEAT_ELN<-BEAT_clinical_final[,c(1,113)]
colnames(BEAT_ELN)[1]<-'rownames'
BEAT.ref.anno.n<-right_join(BEAT.ref,BEAT.ref.anno,by='rownames')
BEAT.ref.anno.n<-inner_join(BEAT.ref.anno.n,BEAT_ELN,by='rownames')
colnames(BEAT.ref.anno.n)[3]<-'ELN2017 Risk Group'

TCGA_ELN<-TCGA_clinical_final[,c(1,2102)]
colnames(TCGA_ELN)[1]<-'rownames'
TCGA.ref.anno.n<-right_join(TCGA.ref,TCGA.ref.anno,by='rownames')
TCGA.ref.anno.n<-inner_join(TCGA.ref.anno.n,TCGA_ELN,by='rownames')
colnames(TCGA.ref.anno.n)[3]<-'ELN2017 Risk Group'

BEAT.ref.anno.n.clean<-BEAT.ref.anno.n %>% column_to_rownames('rownames')
TCGA.ref.anno.n.clean<-TCGA.ref.anno.n %>% column_to_rownames('rownames')

BEAT_anno.all<-HeatmapAnnotation(df=BEAT.ref.anno.n.clean,
                                  col = list(cls=c('TP53mut'='red','TP53wt'='blue1','TP53mut-like'='purple'),
                                             `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                 simple_anno_size = unit(0.3, "cm"))
BEAT_anno.mtlwt<-HeatmapAnnotation(df=BEAT.ref.anno.n.clean[c(BEAT.ref.anno.n.clean$cls=='TP53mut-like'|BEAT.ref.anno.n.clean$cls=='TP53wt'),],
                                 col = list(cls=c('TP53wt'='blue1','TP53mut-like'='purple'),
                                            `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                 simple_anno_size = unit(0.3, "cm"))
TCGA_anno.all<-HeatmapAnnotation(df=TCGA.ref.anno.n.clean,
                                 col = list(cls=c('TP53mut'='red','TP53wt'='blue1','TP53mut-like'='purple'),
                                            `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                 simple_anno_size = unit(0.3, "cm"))
TCGA_anno.mtlwt<-HeatmapAnnotation(df=TCGA.ref.anno.n.clean[c(TCGA.ref.anno.n.clean$cls=='TP53mut-like'|TCGA.ref.anno.n.clean$cls=='TP53wt'),],
                                   col = list(cls=c('TP53wt'='blue1','TP53mut-like'='purple'),
                                              `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                   simple_anno_size = unit(0.3, "cm"))

library(circlize)
col_fun = colorRamp2(c(0,1), c("grey", "black"),transparency = 0.35)

BEAT.cyto.all.binary<-Heatmap(as.matrix(BEAT_clinical_cyto_fin_v02.df.t)[,c(colnames(beat.mt.mat.sort),colnames(beat.mtl.mat.sort),colnames(beat.wt.mat.sort))],
        show_row_names = T,show_column_names = F,top_annotation = BEAT_anno.all,
        column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 12),
        show_row_dend = F,cluster_columns =F,
        cluster_rows = F,
        col=col_fun,
        na_col = 'white',
        rect_gp = gpar(col = "white", lwd = 0.4),show_heatmap_legend = F)
TCGA.cyto.all.binary<-Heatmap(as.matrix(TCGA_clinical_cyto_fin_v02.df.t)[,c(colnames(tcga.mt.mat.sort),colnames(tcga.mtl.mat.sort),colnames(tcga.wt.mat.sort))],
        show_row_names = T,show_column_names = F,top_annotation = TCGA_anno.all,
        column_names_max_height = unit(4, "cm"),
        column_names_gp = gpar(fontsize = 6),
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 12),
        show_row_dend = F,cluster_columns =F,
        cluster_rows = F,
        col=col_fun,
        na_col = 'white',
        rect_gp = gpar(col = "white", lwd = 0.4))

BEAT.cyto.mtlwt.binary<-Heatmap(as.matrix(BEAT_clinical_cyto_fin_v02.df.t)[,c(colnames(beat.mtl.mat.sort),colnames(beat.wt.mat.sort))],
                              show_row_names = T,show_column_names = F,top_annotation = BEAT_anno.mtlwt,
                              column_names_max_height = unit(4, "cm"),
                              column_names_gp = gpar(fontsize = 6),
                              row_names_max_width = unit(6, "cm"),
                              row_names_gp = gpar(fontsize = 12),
                              show_row_dend = F,cluster_columns =F,
                              cluster_rows = F,
                              col=col_fun,
                              na_col = 'white',
                              rect_gp = gpar(col = "white", lwd = 0.4))
TCGA.cyto.mtlwt.binary<-Heatmap(as.matrix(TCGA_clinical_cyto_fin_v02.df.t)[,c(colnames(tcga.mtl.mat.sort),colnames(tcga.wt.mat.sort))],
                              show_row_names = T,show_column_names = F,top_annotation = TCGA_anno.mtlwt,
                              column_names_max_height = unit(4, "cm"),
                              column_names_gp = gpar(fontsize = 6),
                              row_names_max_width = unit(6, "cm"),
                              row_names_gp = gpar(fontsize = 12),
                              show_row_dend = F,cluster_columns =F,
                              cluster_rows = F,
                              col=col_fun,
                              na_col = 'white',
                              rect_gp = gpar(col = "white", lwd = 0.4))
#map only: for publication figure generation purposes only
BEAT_anno.all<-HeatmapAnnotation(df=BEAT.ref.anno.n.clean,
                                 col = list(cls=c('TP53mut'='red','TP53wt'='blue1','TP53mut-like'='purple'),
                                            `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                 show_legend = c(F,F),show_annotation_name = F,simple_anno_size = unit(0.3, "cm"))
TCGA_anno.all<-HeatmapAnnotation(df=TCGA.ref.anno.n.clean,
                                 col = list(cls=c('TP53mut'='red','TP53wt'='blue1','TP53mut-like'='purple'),
                                            `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                 show_legend = c(F,F),show_annotation_name = F,simple_anno_size = unit(0.3, "cm"))


BEAT_anno.mtlwt<-HeatmapAnnotation(df=BEAT.ref.anno.n.clean[c(BEAT.ref.anno.n.clean$cls=='TP53mut-like'|BEAT.ref.anno.n.clean$cls=='TP53wt'),],
                                   col = list(cls=c('TP53wt'='blue1','TP53mut-like'='purple'),
                                              `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                   show_legend = c(F,F),show_annotation_name = F,simple_anno_size = unit(0.3, "cm"))
TCGA_anno.mtlwt<-HeatmapAnnotation(df=TCGA.ref.anno.n.clean[c(TCGA.ref.anno.n.clean$cls=='TP53mut-like'|TCGA.ref.anno.n.clean$cls=='TP53wt'),],
                                   col = list(cls=c('TP53wt'='blue1','TP53mut-like'='purple'),
                                              `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                   show_legend = c(F,F),show_annotation_name = F,simple_anno_size = unit(0.3, "cm"))


BEAT.cyto.all.binary.maponly<-Heatmap(as.matrix(BEAT_clinical_cyto_fin_v02.df.t)[,c(colnames(beat.mt.mat.sort),colnames(beat.mtl.mat.sort),colnames(beat.wt.mat.sort))],
                                      show_row_names = F,show_column_names = F,top_annotation = BEAT_anno.all,
                                      column_names_max_height = unit(4, "cm"),
                                      column_names_gp = gpar(fontsize = 6),
                                      row_names_max_width = unit(6, "cm"),
                                      row_names_gp = gpar(fontsize = 12),
                                      show_row_dend = F,cluster_columns =F,
                                      cluster_rows = F,
                                      col=col_fun,
                                      na_col = 'white',
                                      rect_gp = gpar(col = "white", lwd = 0.4),show_heatmap_legend = F)

TCGA.cyto.all.binary.maponly<-Heatmap(as.matrix(TCGA_clinical_cyto_fin_v02.df.t)[,c(colnames(tcga.mt.mat.sort),colnames(tcga.mtl.mat.sort),colnames(tcga.wt.mat.sort))],
                                      show_row_names = F,show_column_names = F,top_annotation = TCGA_anno.all,
                                      column_names_max_height = unit(4, "cm"),
                                      column_names_gp = gpar(fontsize = 6),
                                      row_names_max_width = unit(6, "cm"),
                                      row_names_gp = gpar(fontsize = 12),
                                      show_row_dend = F,cluster_columns =F,
                                      cluster_rows = F,
                                      col=col_fun,
                                      na_col = 'white',
                                      rect_gp = gpar(col = "white", lwd = 0.4),show_heatmap_legend = F)


BEAT.cyto.mtlwt.binary.maponly<-Heatmap(as.matrix(BEAT_clinical_cyto_fin_v02.df.t)[,c(colnames(beat.mtl.mat.sort),colnames(beat.wt.mat.sort))],
                                show_row_names = F,show_column_names = F,top_annotation = BEAT_anno.mtlwt,
                                column_names_max_height = unit(4, "cm"),
                                column_names_gp = gpar(fontsize = 6),
                                row_names_max_width = unit(6, "cm"),
                                row_names_gp = gpar(fontsize = 12),
                                show_row_dend = F,cluster_columns =F,
                                cluster_rows = F,
                                col=col_fun,
                                na_col = 'white',
                                rect_gp = gpar(col = "white", lwd = 0.4),show_heatmap_legend = F)
TCGA.cyto.mtlwt.binary.maponly<-Heatmap(as.matrix(TCGA_clinical_cyto_fin_v02.df.t)[,c(colnames(tcga.mtl.mat.sort),colnames(tcga.wt.mat.sort))],
                                show_row_names = F,show_column_names = F,top_annotation = TCGA_anno.mtlwt,
                                column_names_max_height = unit(4, "cm"),
                                column_names_gp = gpar(fontsize = 6),
                                row_names_max_width = unit(6, "cm"),
                                row_names_gp = gpar(fontsize = 12),
                                show_row_dend = F,cluster_columns =F,
                                cluster_rows = F,
                                col=col_fun,
                                na_col = 'white',
                                rect_gp = gpar(col = "white", lwd = 0.4),show_heatmap_legend = F)

setwd('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/')
jpeg("Clinical/Categorical parameter analysis/Binary plot/BEAT_cyto_ver02_binary_all.jpeg",width=4.5,height=3.3,units='in',res=300)
BEAT.cyto.all.binary
dev.off()
jpeg('Clinical/Categorical parameter analysis/Binary plot/TCGA_cyto_ver02_binary_all.jpeg',width=4.5,height=3.8,units='in',res=300)
TCGA.cyto.all.binary
dev.off()

jpeg('Clinical/Categorical parameter analysis/Binary plot/BEAT_cyto_ver02_binary_mtlwt.jpeg',width=4.5,height=3.8,units='in',res=300)
BEAT.cyto.mtlwt.binary
dev.off()
jpeg('Clinical/Categorical parameter analysis/Binary plot/TCGA_cyto_ver02_binary_mtlwt.jpeg.jpeg',width=4.5,height=3.8,units='in',res=300)
TCGA.cyto.mtlwt.binary
dev.off()
#heatmap only
jpeg("Clinical/Categorical parameter analysis/Binary plot/BEAT_cyto_ver02_binary_all_heatmap.jpeg",width=4.5,height=3.3,units='in',res=300)
BEAT.cyto.all.binary.maponly
dev.off()
jpeg('Clinical/Categorical parameter analysis/Binary plot/TCGA_cyto_ver02_binary_all_heatmap.jpeg',width=4.5,height=3.8,units='in',res=300)
TCGA.cyto.all.binary.maponly
dev.off()

jpeg('Clinical/Categorical parameter analysis/Binary plot/BEAT_cyto_ver02_binary_mtlwt_heatmap.jpeg',width=4.5,height=3.8,units='in',res=300)
BEAT.cyto.mtlwt.binary.maponly
dev.off()
jpeg('Clinical/Categorical parameter analysis/Binary plot/TCGA_cyto_ver02_binary_mtlwt_heatmap.jpeg.jpeg',width=4.5,height=3.8,units='in',res=300)
TCGA.cyto.mtlwt.binary.maponly
dev.off()
#################################################################
#get mutational profile

library(dplyr)
#TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/TCGA_clinical_final_class.rds')
BEAT_clinical_final.old<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/BEAT_TCGA_analysis_WES/GLM/GEP/class ids/BEAT_clinical_final_class.rds')
BEAT_clinical_final.old_403<-BEAT_clinical_final.old[match(BEAT_clinical_final$LabId,BEAT_clinical_final.old$LabId),]
###
BEAT_clinical_final<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')

fav_q_ix<-which(BEAT_clinical_final$`ELN Risk Group_YL`=='?Fav')
BEAT_clinical_final$ELN_Risk_Group_YL_new<-BEAT_clinical_final$`ELN Risk Group_YL`

for(i in (fav_q_ix)){
  cat<-BEAT_clinical_final$ELN2017[i]
  if(cat=='Favorable'){
    BEAT_clinical_final$ELN_Risk_Group_YL_new[i]<-'Favorable'}
  # }else if(cat=='Intermediate'){
  #   BEAT_clinical_final$ELN_Risk_Group_YL_new[i]<-'Intermediate'
  # }else if(cat=='Adverse'){
  #   BEAT_clinical_final$ELN_Risk_Group_YL_new[i]<-'Adverse'
  # }
}
BEAT_clinical_final$`ELN_Risk_Group_YL_new`[BEAT_clinical_final$`ELN_Risk_Group_YL_new`=='?Fav']<-'unknown'
TCGA_clinical_final$`ELN Risk Group_YL`[TCGA_clinical_final$`ELN Risk Group_YL`=='?Fav']<-'unknown'




#mutational status
#need to preprocess mutational subtypes FLT3,PTPN11,SF3B1(21:23) for 1,0 and NA
table(BEAT_clinical_final.old_403$FLT3)
table(BEAT_clinical_final.old_403$WES_PTPN11)
table(BEAT_clinical_final.old_403$WES_SF3B1)

flt3_itd<-ifelse(BEAT_clinical_final.old_403$`FLT3-ITD`=='positive',1,
                 ifelse(BEAT_clinical_final.old_403$`FLT3-ITD`=='negative',0,NA))
npm1<-ifelse(BEAT_clinical_final.old_403$NPM1=='positive',1,
             ifelse(BEAT_clinical_final.old_403$NPM1=='negative',0,NA))
# flt3<-ifelse(BEAT_LASSO_final_res.clinical.n$FLT3=='negative',0,
#              ifelse(BEAT_LASSO_final_res.clinical.n$FLT3!='negative',1,NA))
# ptpn11<-ifelse(is.na(BEAT_LASSO_final_res.clinical.n$PTPN11),NA,1)
# sf3b1<-ifelse(is.na(BEAT_LASSO_final_res.clinical.n$SF3B1),NA,1)
BEAT_clinical_final.old_403$`FLT3-ITD`<-flt3_itd
BEAT_clinical_final.old_403$NPM1<-npm1
# BEAT_LASSO_final_res.clinical.n$FLT3<-flt3
# BEAT_LASSO_final_res.clinical.n$PTPN11<-ptpn11
# BEAT_LASSO_final_res.clinical.n$SF3B1<-sf3b1

BEAT_LASSO_final_res_cat<-BEAT_clinical_final.old_403[,c(4,6,9,11,12,14,19,21:58)]
BEAT_LASSO_final_res_cat<-right_join(BEAT_LASSO_final_res_cat,BEAT_clinical_final[,c(1,2,113)],by='LabId')
#BEAT_LASSO_final_res_cat[,c(8:44)]<-data.frame(lapply(BEAT_LASSO_final_res_cat[,c(8:44)],as.character))
#BEAT_LASSO_final_cat<-BEAT_LASSO_final_res_cat[,-c(8:44)]
BEAT_LASSO_final_mut<-BEAT_LASSO_final_res_cat[,c(1,46,8:44)]
Mut_list<-c('DNMT3A','TET2','IDH1','IDH2','WT1','EZH2','ASXL1','SRSF2','STAG2','RAD21',
            'FLT3','KIT','KRAS','NRAS','PTPN11','NF1','JAK2','BCOR','GATA2','ETV6',
            'CEBPA','RUNX1','CDKN2A','NPM1','MYC','U2AF1','PHF6','SMC3','SMC1A','CSF3R',
            'PDS5B','SF3B1','BRINP3','CROCC','PLCE1','SPEN')
Mut_list.n<-str_c('WES',Mut_list,sep='_')
#use clinical summary data for NPM1 and FLT3-ITD - these two are complete
Mut_list.n[Mut_list.n =='WES_NPM1']<-'NPM1'
# Mut_list.n[Mut_list.n =='TP53_mut_stat']<-'TP53'
Mut_list.n<-c(Mut_list.n,'FLT3-ITD')
#################
#get TCGA clinical data - not many data to look at
#TCGA_clinical_final<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')
TCGA_clinical.final.n<-TCGA_clinical_final


#categorical results
#make it long form
Mut_list<-c('DNMT3A','TET2','IDH1','IDH2','WT1','EZH2','ASXL1','SRSF2','STAG2','RAD21',
            'FLT3','KIT','KRAS','NRAS','PTPN11','NF1','JAK2','BCOR','GATA2','ETV6',
            'CEBPA','RUNX1','CDKN2A','NPM1','MYC','U2AF1','PHF6','SMC3','SMC1A','CSF3R',
            'PDS5B','SF3B1','BRINP3','CROCC','PLCE1','SPEN',)
#clinical.rna.GDC.178.final.mut<-TCGA_clinical.final.n[,which(colnames(TCGA_clinical.final.n)%in%Mut_list)]

wes.ix01<-which(c(Mut_list)%in%colnames(TCGA_clinical.final.n) )
wes.ix01.id<-c(Mut_list,'TP53')[wes.ix01]
df01<-as.data.frame(matrix(0,dim(TCGA_clinical.final.n)[1],length(wes.ix01.id)))

rownames(df01)<-TCGA_clinical.final.n$patient_RNA_id
colnames(df01)<-wes.ix01.id

for(i in 1:length((wes.ix01.id))){
  ix<-which(colnames(TCGA_clinical.final.n) %in% wes.ix01.id[i])
  mut_stat<-TCGA_clinical.final.n[,ix]
  df01[,i]<-ifelse(is.na(mut_stat),0,1)
}


#BEAT_LASSO_final_res_cat<-BEAT_LASSO_final_res.clinical.n[,c(1,5,8,10,11,18,20:24,26)]
# TCGA_clinical.final.cat<-TCGA_clinical.final.n[,c(1,14,25,27,35,2067,2098,2099:2102)]
TCGA_clinical.final.cat<-TCGA_clinical.final.n[,c(1,15,26,27,36,2014,2068,2088,2101,2102)]
TCGA_clinical.final.cat.mut<-right_join(TCGA_clinical.final.cat,
                                        rownames_to_column(df01,var='patient_RNA_id'),by='patient_RNA_id')
TCGA_clinical.final.mut<-TCGA_clinical.final.cat.mut[,-c(2:8,10)]
#TCGA_clinical.final.cat<-TCGA_clinical.final.cat.mut[,c(2:8,10)]
#TCGA_clinical.final.mut[,c(3:36)]<-data.frame(lapply(TCGA_clinical.final.mut[,c(3:36)],as.character))
#make heatmap
BEAT_mut<-BEAT_LASSO_final_mut %>% column_to_rownames(var='LabId')
TCGA_mut<-TCGA_clinical.final.mut%>% column_to_rownames(var='patient_RNA_id')

for (i in 1:dim(BEAT_mut)[1]){
  cls<-names(table(BEAT_mut$cls.GLM))
  if(BEAT_mut$cls.GLM[i]=='TP53_WT'){#TP53wt
    BEAT_mut$cls.GLM[i]<-'TP53wt'
  }
  if(BEAT_mut$cls.GLM[i]=='TP53_MUT'){#TP53mut
    BEAT_mut$cls.GLM[i]<-'TP53mut'
  }
  if(BEAT_mut$cls.GLM[i]=='TP53_MUTlike'){#TP53mut-like
    BEAT_mut$cls.GLM[i]<-'TP53mut-like'
  }
}
for (i in 1:dim(TCGA_mut)[1]){
  cls<-names(table(TCGA_mut$cls.GLM.x))
  if(TCGA_mut$cls.GLM.y[i]=='TP53_WT'){#TP53wt
    TCGA_mut$cls.GLM.y[i]<-'TP53wt'
  }
  if(TCGA_mut$cls.GLM.y[i]=='TP53_MUT'){#TP53mut
    TCGA_mut$cls.GLM.y[i]<-'TP53mut'
  }
  if(TCGA_mut$cls.GLM.y[i]=='TP53MUT_like'){#TP53mut-like
    TCGA_mut$cls.GLM.y[i]<-'TP53mut-like'
  }
}

binary_sorting<-function(mat){
  #mdf.n<-Final.DEG.mat.n.t[,which(Final.DEG.mat.n.anno$cls=='TP53MUT_like')]
  mdf.n<-mat
  mdf.n = mdf.n[order(mdf.n[, ncol(mdf.n)], decreasing = TRUE), ]
  #colnames(mdf.n) = gsub(pattern = "^X", replacement = "", colnames(mdf.n))
  nMut = mdf.n[, ncol(mdf.n)]
  
  #mdf.n = mdf.n[, -ncol(mdf.n)]
  
  mdf.n.temp.copy = mdf.n #temp copy of original unsorted numeric coded matrix
  
  mdf.n[mdf.n != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
  tmdf.n = t(mdf.n) #transposematrix
  mdf.n = t(tmdf.n[do.call(order, c(as.list(as.data.frame(tmdf.n)), decreasing = TRUE)), ]) #sort
  
  mdf.n.temp.copy = mdf.n.temp.copy[rownames(mdf.n),] #organise original matrix into sorted matrix
  mdf.n.temp.copy = mdf.n.temp.copy[,colnames(mdf.n)]
  mdf.n = mdf.n.temp.copy
  return(mdf.n)
}
BEAT_mut.t<-as.data.frame(t(BEAT_mut[,2:dim(BEAT_mut)[2]]))
Mut_names<-c(unlist(str_split(rownames(BEAT_mut.t),'_')[c(1,2)]),str_split(rownames(BEAT_mut.t),'_')[3:dim(BEAT_mut.t)[1]]%>%map_chr(2))
rownames(BEAT_mut.t)<-Mut_names
TCGA_mut.t<-as.data.frame(t(TCGA_mut[,2:dim(TCGA_mut)[2]]))

beat.mtl.mut.sort<-binary_sorting(BEAT_mut.t[,which(BEAT_mut$cls.GLM=='TP53mut-like')])
beat.mt.mut.sort<-binary_sorting(BEAT_mut.t[,which(BEAT_mut$cls.GLM=='TP53mut')])
beat.wt.mut.sort<-binary_sorting(BEAT_mut.t[,which(BEAT_mut$cls.GLM=='TP53wt')])

tcga.mtl.mut.sort<-binary_sorting(TCGA_mut.t[,which(TCGA_mut$cls.GLM.y=='TP53mut-like')])
tcga.mt.mut.sort<-binary_sorting(TCGA_mut.t[,which(TCGA_mut$cls.GLM.y=='TP53mut')])
tcga.wt.mut.sort<-binary_sorting(TCGA_mut.t[,which(TCGA_mut$cls.GLM.y=='TP53wt')])

library(circlize)
col_fun = colorRamp2(c(0,1), c("grey", "black"),transparency = 0.35)
BEAT.mut.ref<-data.frame('rownames'=colnames(BEAT_mut.t[,c(colnames(beat.mt.mut.sort),colnames(beat.mtl.mut.sort),colnames(beat.wt.mut.sort))]))
TCGA.mut.ref<-data.frame('rownames'=colnames(TCGA_mut.t[,c(colnames(tcga.mt.mut.sort),colnames(tcga.mtl.mut.sort),colnames(tcga.wt.mut.sort))]))


BEAT.mut.ref.anno<-data.frame('cls'=BEAT_mut$cls.GLM,
                          row.names = rownames(BEAT_mut))%>%rownames_to_column(var='rownames')
TCGA.mut.ref.anno<-data.frame('cls'=TCGA_mut$cls.GLM.y,
                          row.names = rownames(TCGA_mut))%>%rownames_to_column(var='rownames')
BEAT_ELN<-BEAT_clinical_final[,c(1,113)]
colnames(BEAT_ELN)[1]<-'rownames'
BEAT.ref.anno.n<-right_join(BEAT.mut.ref,BEAT.mut.ref.anno,by='rownames')
BEAT.ref.anno.n<-inner_join(BEAT.ref.anno.n,BEAT_ELN,by='rownames')
colnames(BEAT.ref.anno.n)[3]<-'ELN2017 Risk Group'

TCGA_ELN<-TCGA_clinical_final[,c(1,2102)]
colnames(TCGA_ELN)[1]<-'rownames'
TCGA.ref.anno.n<-right_join(TCGA.mut.ref,TCGA.mut.ref.anno,by='rownames')
TCGA.ref.anno.n<-inner_join(TCGA.ref.anno.n,TCGA_ELN,by='rownames')
colnames(TCGA.ref.anno.n)[3]<-'ELN2017 Risk Group'

BEAT.ref.anno.n.clean<-BEAT.ref.anno.n %>% column_to_rownames('rownames')
TCGA.ref.anno.n.clean<-TCGA.ref.anno.n %>% column_to_rownames('rownames')

BEAT_anno.all<-HeatmapAnnotation(df=BEAT.ref.anno.n.clean,
                                 col = list(cls=c('TP53mut'='red','TP53wt'='blue1','TP53mut-like'='purple'),
                                            `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                 simple_anno_size = unit(0.3, "cm"))
BEAT_anno.mtlwt<-HeatmapAnnotation(df=BEAT.ref.anno.n.clean[c(BEAT.ref.anno.n.clean$cls=='TP53mut-like'|BEAT.ref.anno.n.clean$cls=='TP53wt'),],
                                   col = list(cls=c('TP53wt'='blue1','TP53mut-like'='purple'),
                                              `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                   simple_anno_size = unit(0.3, "cm"))
TCGA_anno.all<-HeatmapAnnotation(df=TCGA.ref.anno.n.clean,
                                 col = list(cls=c('TP53mut'='red','TP53wt'='blue1','TP53mut-like'='purple'),
                                            `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                 simple_anno_size = unit(0.3, "cm"))
TCGA_anno.mtlwt<-HeatmapAnnotation(df=TCGA.ref.anno.n.clean[c(TCGA.ref.anno.n.clean$cls=='TP53mut-like'|TCGA.ref.anno.n.clean$cls=='TP53wt'),],
                                   col = list(cls=c('TP53wt'='blue1','TP53mut-like'='purple'),
                                              `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                   simple_anno_size = unit(0.3, "cm"))

library(circlize)
col_fun = colorRamp2(c(0,1), c("grey", "black"),transparency = 0.35)

BEAT.mut.all.binary<-Heatmap(as.matrix(BEAT_mut.t)[,c(colnames(beat.mt.mut.sort),colnames(beat.mtl.mut.sort),colnames(beat.wt.mut.sort))],
                              show_row_names = T,show_column_names = F,top_annotation = BEAT_anno.all,
                              column_names_max_height = unit(4, "cm"),
                              column_names_gp = gpar(fontsize = 6),
                              row_names_max_width = unit(6, "cm"),
                              row_names_gp = gpar(fontsize = 8),
                              show_row_dend = F,cluster_columns =F,
                              cluster_rows = F,
                              col=col_fun,
                              na_col = 'white',
                              rect_gp = gpar(col = "white", lwd = 0.4),show_heatmap_legend = F)
TCGA.mut.all.binary<-Heatmap(as.matrix(TCGA_mut.t)[,c(colnames(tcga.mt.mut.sort),colnames(tcga.mtl.mut.sort),colnames(tcga.wt.mut.sort))],
                              show_row_names = T,show_column_names = F,top_annotation = TCGA_anno.all,
                              column_names_max_height = unit(4, "cm"),
                              column_names_gp = gpar(fontsize = 6),
                              row_names_max_width = unit(6, "cm"),
                              row_names_gp = gpar(fontsize = 8),
                              show_row_dend = F,cluster_columns =F,
                              cluster_rows = F,
                              col=col_fun,
                              na_col = 'white',
                              rect_gp = gpar(col = "white", lwd = 0.4))

BEAT.mut.mtlwt.binary<-Heatmap(as.matrix(BEAT_mut.t)[,c(colnames(beat.mtl.mut.sort),colnames(beat.wt.mut.sort))],
                                show_row_names = T,show_column_names = F,top_annotation = BEAT_anno.mtlwt,
                                column_names_max_height = unit(4, "cm"),
                                column_names_gp = gpar(fontsize = 6),
                                row_names_max_width = unit(6, "cm"),
                                row_names_gp = gpar(fontsize = 8),
                                show_row_dend = F,cluster_columns =F,
                                cluster_rows = F,
                                col=col_fun,
                                na_col = 'white',
                                rect_gp = gpar(col = "white", lwd = 0.4))
TCGA.mut.mtlwt.binary<-Heatmap(as.matrix(TCGA_mut.t)[,c(colnames(tcga.mtl.mut.sort),colnames(tcga.wt.mut.sort))],
                                show_row_names = T,show_column_names = F,top_annotation = TCGA_anno.mtlwt,
                                column_names_max_height = unit(4, "cm"),
                                column_names_gp = gpar(fontsize = 6),
                                row_names_max_width = unit(6, "cm"),
                                row_names_gp = gpar(fontsize = 8),
                                show_row_dend = F,cluster_columns =F,
                                cluster_rows = F,
                                col=col_fun,
                                na_col = 'white',
                                rect_gp = gpar(col = "white", lwd = 0.4))
#map only: for publication figure generation purposes only
BEAT_anno.all<-HeatmapAnnotation(df=BEAT.ref.anno.n.clean,
                                 col = list(cls=c('TP53mut'='red','TP53wt'='blue1','TP53mut-like'='purple'),
                                            `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                 show_legend = c(F,F),show_annotation_name = F,simple_anno_size = unit(0.3, "cm"))
TCGA_anno.all<-HeatmapAnnotation(df=TCGA.ref.anno.n.clean,
                                 col = list(cls=c('TP53mut'='red','TP53wt'='blue1','TP53mut-like'='purple'),
                                            `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                 show_legend = c(F,F),show_annotation_name = F,simple_anno_size = unit(0.3, "cm"))


BEAT_anno.mtlwt<-HeatmapAnnotation(df=BEAT.ref.anno.n.clean[c(BEAT.ref.anno.n.clean$cls=='TP53mut-like'|BEAT.ref.anno.n.clean$cls=='TP53wt'),],
                                   col = list(cls=c('TP53wt'='blue1','TP53mut-like'='purple'),
                                              `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                   show_legend = c(F,F),show_annotation_name = F,simple_anno_size = unit(0.3, "cm"))
TCGA_anno.mtlwt<-HeatmapAnnotation(df=TCGA.ref.anno.n.clean[c(TCGA.ref.anno.n.clean$cls=='TP53mut-like'|TCGA.ref.anno.n.clean$cls=='TP53wt'),],
                                   col = list(cls=c('TP53wt'='blue1','TP53mut-like'='purple'),
                                              `ELN2017 Risk Group`=c('Adverse'='orange2','Favorable'='green3','Intermediate'='yellow2','unknown'='white')),
                                   show_legend = c(F,F),show_annotation_name = F,simple_anno_size = unit(0.3, "cm"))


BEAT.mut.all.binary.maponly<-Heatmap(as.matrix(BEAT_mut.t)[,c(colnames(beat.mt.mut.sort),colnames(beat.mtl.mut.sort),colnames(beat.wt.mut.sort))],
                                      show_row_names = F,show_column_names = F,top_annotation = BEAT_anno.all,
                                      column_names_max_height = unit(4, "cm"),
                                      column_names_gp = gpar(fontsize = 6),
                                      row_names_max_width = unit(6, "cm"),
                                      row_names_gp = gpar(fontsize = 12),
                                      show_row_dend = F,cluster_columns =F,
                                      cluster_rows = F,
                                      col=col_fun,
                                      na_col = 'white',
                                      rect_gp = gpar(col = "white", lwd = 0.4),show_heatmap_legend = F)

TCGA.mut.all.binary.maponly<-Heatmap(as.matrix(TCGA_mut.t)[,c(colnames(tcga.mt.mut.sort),colnames(tcga.mtl.mut.sort),colnames(tcga.wt.mut.sort))],
                                      show_row_names = F,show_column_names = F,top_annotation = TCGA_anno.all,
                                      column_names_max_height = unit(4, "cm"),
                                      column_names_gp = gpar(fontsize = 6),
                                      row_names_max_width = unit(6, "cm"),
                                      row_names_gp = gpar(fontsize = 12),
                                      show_row_dend = F,cluster_columns =F,
                                      cluster_rows = F,
                                      col=col_fun,
                                      na_col = 'white',
                                      rect_gp = gpar(col = "white", lwd = 0.4),show_heatmap_legend = F)


BEAT.mut.mtlwt.binary.maponly<-Heatmap(as.matrix(BEAT_mut.t)[,c(colnames(beat.mtl.mut.sort),colnames(beat.wt.mut.sort))],
                                        show_row_names = F,show_column_names = F,top_annotation = BEAT_anno.mtlwt,
                                        column_names_max_height = unit(4, "cm"),
                                        column_names_gp = gpar(fontsize = 6),
                                        row_names_max_width = unit(6, "cm"),
                                        row_names_gp = gpar(fontsize = 12),
                                        show_row_dend = F,cluster_columns =F,
                                        cluster_rows = F,
                                        col=col_fun,
                                        na_col = 'white',
                                        rect_gp = gpar(col = "white", lwd = 0.4),show_heatmap_legend = F)
TCGA.mut.mtlwt.binary.maponly<-Heatmap(as.matrix(TCGA_mut.t)[,c(colnames(tcga.mtl.mut.sort),colnames(tcga.wt.mut.sort))],
                                        show_row_names = F,show_column_names = F,top_annotation = TCGA_anno.mtlwt,
                                        column_names_max_height = unit(4, "cm"),
                                        column_names_gp = gpar(fontsize = 6),
                                        row_names_max_width = unit(6, "cm"),
                                        row_names_gp = gpar(fontsize = 12),
                                        show_row_dend = F,cluster_columns =F,
                                        cluster_rows = F,
                                        col=col_fun,
                                        na_col = 'white',
                                        rect_gp = gpar(col = "white", lwd = 0.4),show_heatmap_legend = F)

setwd('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/')
jpeg("Clinical/Categorical parameter analysis/Binary plot/Mutataion/BEAT_mut_binary_all.jpeg",width=4.5,height=3.3,units='in',res=300)
BEAT.mut.all.binary
dev.off()
jpeg('Clinical/Categorical parameter analysis/Binary plot/Mutataion/TCGA_mut_binary_all.jpeg',width=4.5,height=3.8,units='in',res=300)
TCGA.mut.all.binary
dev.off()

jpeg('Clinical/Categorical parameter analysis/Binary plot/Mutataion/BEAT_mut_binary_mtlwt.jpeg',width=4.5,height=3.8,units='in',res=300)
BEAT.mut.mtlwt.binary
dev.off()
jpeg('Clinical/Categorical parameter analysis/Binary plot/Mutataion/TCGA_mut_binary_mtlwt.jpeg.jpeg',width=4.5,height=3.8,units='in',res=300)
TCGA.mut.mtlwt.binary
dev.off()
#heatmap only
jpeg("Clinical/Categorical parameter analysis/Binary plot/Mutataion/BEAT_mut_binary_all_heatmap.jpeg",width=4.5,height=3.3,units='in',res=300)
BEAT.mut.all.binary.maponly
dev.off()
jpeg('Clinical/Categorical parameter analysis/Binary plot/Mutataion/TCGA_mut_binary_all_heatmap.jpeg',width=4.5,height=3.8,units='in',res=300)
TCGA.mut.all.binary.maponly
dev.off()

jpeg('Clinical/Categorical parameter analysis/Binary plot/Mutataion/BEAT_mut_binary_mtlwt_heatmap.jpeg',width=4.5,height=3.8,units='in',res=300)
BEAT.mut.mtlwt.binary.maponly
dev.off()
jpeg('Clinical/Categorical parameter analysis/Binary plot/Mutataion/TCGA_mut_binary_mtlwt_heatmap.jpeg.jpeg',width=4.5,height=3.8,units='in',res=300)
TCGA.mut.mtlwt.binary.maponly
dev.off()



