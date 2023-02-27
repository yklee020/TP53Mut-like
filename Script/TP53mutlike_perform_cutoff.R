#get the precision and accuracy in different cutoffs
#1. Ridge model01 (performance TP53Mut classification): AUROC and AUPRC
#2. Ridge model02 (performance TP53Mut classification): AUROC and AUPRC
#3. 25 signature genes()
#get the function to draw AUC 

source("C:/Users/yklee/Desktop/Sachs Lab/R_project/TP53/TP53GEP_functions.r")
packages<-c('edgeR','limma','dplyr','tidyr','sva','tidyverse',
            'ggplot2','glmnet','caret','pROC','survival','survminer','HGNChelper')

for(p in packages){
  library(p,character.only = T)
}
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/")
# model01.trn.y<-trn.y
# model01.trn.prob<-prob.ridge.res
# model01.test.y<-test.y
# model01.test.prob<-prob.ridge.res.test
# model01.tcga.y<-tcga.y
# model01.tcga.prob<-prob.ridge.res.tcga
# 
# saveRDS(model01.trn.y,'model01_trn_y.rds')
# saveRDS(model01.trn.prob,'model01_trn_prob.rds')
# saveRDS(model01.test.y,'model01_test_y.rds')
# saveRDS(model01.test.prob,'model01_test_prob.rds')
# saveRDS(model01.tcga.y,'model01_tcga_y.rds')
# saveRDS(model01.tcga.prob,'model01_tcga_prob.rds')
# 
# model02.trn.y<-trn.y
# model02.trn.prob<-prob.ridge.res
# model02.test.y<-test.y
# model02.test.prob<-prob.ridge.res.test
# model02.tcga.y<-tcga.y
# model02.tcga.prob<-prob.ridge.res.tcga
# setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/")
# 
# saveRDS(model02.trn.y,'model02_trn_y.rds')
# saveRDS(model02.trn.prob,'model02_trn_prob.rds')
# saveRDS(model02.test.y,'model02_test_y.rds')
# saveRDS(model02.test.prob,'model02_test_prob.rds')
# saveRDS(model02.tcga.y,'model02_tcga_y.rds')
# saveRDS(model02.tcga.prob,'model02_tcga_prob.rds')


# sig25_model.trn.y<-trn.y
# sig25_model.trn.prob<-prob.ridge.res
# sig25_model.test.y<-test.y
# sig25_model.test.prob<-prob.ridge.res.test
# sig25_model.tcga.y<-tcga.y
# sig25_model.tcga.prob<-prob.ridge.res.tcga
# 
# setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/")
# saveRDS(sig25_model.trn.y,'sig25_model_trn_y.rds')
# saveRDS(sig25_model.trn.prob,'sig25_model_trn_prob.rds')
# saveRDS(sig25_model.test.y,'sig25_model_test_y.rds')
# saveRDS(sig25_model.test.prob,'sig25_model_test_prob.rds')
# saveRDS(sig25_model.tcga.y,'sig25_model_tcga_y.rds')
# saveRDS(sig25_model.tcga.prob,'sig25_model_tcga_prob.rds')


setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/")

#model01.trn.y<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/model01_trn_y.rds")
#model01.trn.prob<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/model01_trn_prob.rds")
model01.test.y<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/model01_test_y.rds")
model01.test.prob<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/model01_test_prob.rds")
model01.tcga.y<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/model01_tcga_y.rds")
model01.tcga.prob<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/model01_tcga_prob.rds")


#model02.trn.y<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/model02_trn_y.rds")
#model02.trn.prob<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/model02_trn_prob.rds")
model02.test.y<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/model02_test_y.rds")
model02.test.prob<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/model02_test_prob.rds")
model02.tcga.y<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/model02_tcga_y.rds")
model02.tcga.prob<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/model02_tcga_prob.rds")


#sig25_model.trn.y<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/sig25_model_trn_y.rds")
#sig25_model.trn.prob<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/sig25_model_trn_prob.rds")
sig25_model.test.y<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/sig25_model_test_y.rds")
sig25_model.test.prob<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/sig25_model_test_prob.rds")
sig25_model.tcga.y<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/sig25_model_tcga_y.rds")
sig25_model.tcga.prob<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/Analysis/Ridge Regression/Model Performance/sig25_model_tcga_prob.rds")

pROC_summary(model01.test.y,model01.test.prob)
Prec_Recall_summary(model01.test.y,model01.test.prob)
pROC_summary(model01.tcga.y,model01.tcga.prob)
Prec_Recall_summary(model01.tcga.y,model01.tcga.prob)

pROC_summary(model02.test.y,model02.test.prob)
Prec_Recall_summary(model02.test.y,model02.test.prob)
pROC_summary(model02.tcga.y,model02.tcga.prob)
Prec_Recall_summary(model02.tcga.y,model02.tcga.prob)

pROC_summary(sig25_model.test.y,sig25_model.test.prob)
Prec_Recall_summary(sig25_model.test.y,sig25_model.test.prob)
pROC_summary(sig25_model.tcga.y,sig25_model.tcga.prob)
Prec_Recall_summary(sig25_model.tcga.y,sig25_model.tcga.prob)


pROC_summary<-function(y.response,probabilities_res){
  roc(y.response, as.vector(probabilities_res), plot=TRUE, legacy.axes=TRUE, percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", col="#377eb8", lwd=4, print.auc=TRUE)
  roc.info<-roc(y.response, as.vector(probabilities_res), legacy.axes=TRUE)
  roc.df<-data.frame(tpp=roc.info$sensitivities*100, ## tpp = true positive percentage
                     fpp=(1 - roc.info$specificities)*100, ## fpp = false positive precentage
                     thresholds=roc.info$thresholds)
  return(roc.df)
  #return(roc.info)
}

Prec_Recall_summary<-function(y.response,probabilities_res){
  a<-roc(y.response, as.vector(probabilities_res),percent=T,legacy.axes=T)
  PR.df<-coords(a, "all", ret = c("recall", "precision",'threshold'), transpose = F)
  #plot(precision ~ recall,PR.df,type="l", ylim = c(0, 100))
  library(PRROC)
  
  sc1=probabilities_res[y.response==1]
  sc0=probabilities_res[y.response==0]
  
  #b<-pr.curve(scores.class0 = sc1, scores.class1 = sc0,curve=T )
  # tryCatch( { b<-pr.curve(scores.class0 = sc1, scores.class1 = sc0,curve=T ); print('noerror') }
  #           , error = function(e) {b<<-data.frame('auc.integral'=NA)})
  tryCatch( { b<-pr.curve(scores.class0 = sc1, scores.class1 = sc0,curve=T ) }
            , error = function(e) {b<<-data.frame('auc.integral'=NA)})
  plot(pr.curve(scores.class0 = sc1, scores.class1 = sc0,curve=T ))
  return(PR.df)
  #return(b)
}


pROC_summary(model01.test.y,model01.test.prob)
Prec_Recall_summary(model01.test.y,model01.test.prob)
pROC_summary(model01.tcga.y,model01.tcga.prob)
Prec_Recall_summary(model01.tcga.y,model01.tcga.prob)

pROC_summary(model02.test.y,model02.test.prob)
Prec_Recall_summary(model02.test.y,model02.test.prob)
pROC_summary(model02.tcga.y,model02.tcga.prob)
Prec_Recall_summary(model02.tcga.y,model02.tcga.prob)

pROC_summary(sig25_model.test.y,sig25_model.test.prob)
Prec_Recall_summary(sig25_model.test.y,sig25_model.test.prob)
pROC_summary(sig25_model.tcga.y,sig25_model.tcga.prob)
Prec_Recall_summary(sig25_model.tcga.y,sig25_model.tcga.prob)


#get plot and manipulate it
roc_trn.ridge.0<-pROC_summary.0(trn.y,prob.ridge.res)
PR_trn.ridge.0<-Prec_Recall_summary.0(trn.y,prob.ridge.res)

roc_test.ridge.0<-pROC_summary.0(test.y,prob.ridge.res.test)
PR_test.ridge.0<-Prec_Recall_summary.0(test.y,prob.ridge.res.test)



#2022.12.12
#generate figures for the publication
# ridge.trn.roc.obj<-roc(trn.y, as.vector(prob.ridge.res), plot=TRUE, 
#                        legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)
# ridge.trn.roc.auc<-round(auc(trn.y,as.numeric(prob.ridge.res)),3)
# 
# ridge.test.roc.obj<-roc(test.y, as.vector(prob.ridge.res.test), plot=TRUE, 
#                         legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)
# ridge.test.roc.auc<-round(auc(test.y,as.numeric(prob.ridge.res.test)),3)
# 

m01.test.roc.obj<-roc(model01.test.y, as.vector(model01.test.prob), plot=TRUE, 
                        legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)
m01.tcga.roc.obj<-roc(model01.tcga.y, as.vector(model01.tcga.prob), plot=TRUE, 
                      legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)

m02.test.roc.obj<-roc(model02.test.y, as.vector(model02.test.prob), plot=TRUE, 
                      legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)
m02.tcga.roc.obj<-roc(model02.tcga.y, as.vector(model02.tcga.prob), plot=TRUE, 
                      legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)

sig25.test.roc.obj<-roc(sig25_model.test.y, as.vector(sig25_model.test.prob), plot=TRUE, 
                      legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)
sig25.tcga.roc.obj<-roc(sig25_model.tcga.y, as.vector(sig25_model.tcga.prob), plot=TRUE, 
                      legacy.axes=TRUE, percent=F, col="#377eb8", lwd=4, print.auc=T)

pROC_summary

m01.test.roc.0<-pROC_summary.0(model01.test.y,model01.test.prob)
m01.tcga.roc.0<-pROC_summary.0(model01.tcga.y,model01.tcga.prob)
m02.test.roc.0<-pROC_summary.0(model02.test.y,model02.test.prob)
m02.tcga.roc.0<-pROC_summary.0(model02.tcga.y,model02.tcga.prob)
sig25.test.roc.0<-pROC_summary.0(sig25_model.test.y,sig25_model.test.prob)
sig25.tcga.roc.0<-pROC_summary.0(sig25_model.tcga.y,sig25_model.tcga.prob)


m01.test.pr.0<-Prec_Recall_summary.0(model01.test.y,model01.test.prob)
m01.tcga.pr.0<-Prec_Recall_summary.0(model01.tcga.y,model01.tcga.prob)
m02.test.pr.0<-Prec_Recall_summary.0(model02.test.y,model02.test.prob)
m02.tcga.pr.0<-Prec_Recall_summary.0(model02.tcga.y,model02.tcga.prob)
sig25.test.pr.0<-Prec_Recall_summary.0(sig25_model.test.y,sig25_model.test.prob)
sig25.tcga.pr.0<-Prec_Recall_summary.0(sig25_model.tcga.y,sig25_model.tcga.prob)

#precision recall checking
# m01.test.roc.obj
# a<-sapply(seq(0, 1, by=0.005), function(x) coords(m01.test.roc.obj, x, ret=c("recall", "precision")))
# aa<-do.call(rbind.data.frame, t(a))
# aa<-as.data.frame(t(a))[1:176,]
# fscore<-c()
# for(i in 1:dim(aa)[1]){
#   fscore <-c(fscore,(2 * aa$precision[i]* aa$recall[i]) / (aa$precision[i] + aa$recall[i]))
#   
# }
tt<-m01.test.pr.0<-coords(m01.test.roc.obj, "all", ret = c('specificity',"recall", "precision",'threshold'), transpose = F)

m01.test.pr.0<-coords(m01.test.roc.obj, "all", ret = c('specificity',"recall", "precision",'threshold'), transpose = F)
m01.tcga.pr.0<-coords(m01.tcga.roc.obj, "all", ret = c('specificity',"recall", "precision",'threshold'), transpose = F)
m02.test.pr.0<-coords(m02.test.roc.obj, "all", ret = c('specificity',"recall", "precision",'threshold'), transpose = F)
m02.tcga.pr.0<-coords(m02.tcga.roc.obj, "all", ret = c('specificity',"recall", "precision",'threshold'), transpose = F)
sig25.test.pr.0<-coords(sig25.test.roc.obj, "all", ret = c('specificity',"recall", "precision",'threshold'), transpose = F)
sig25.tcga.pr.0<-coords(sig25.tcga.roc.obj, "all", ret = c('specificity',"recall", "precision",'threshold'), transpose = F)


#calculate threshold between sensitivity(TPR)==recall and specificity(1-FPR) using geometric mean
Gmeans.m01.test.ix<-sqrt(m01.test.roc.obj$sensitivities*m01.test.roc.obj$specificities) %>% which.max()
Gmeans.m01.tcga.ix<-sqrt(m01.tcga.roc.obj$sensitivities*m01.tcga.roc.obj$specificities) %>% which.max()

Gmeans.m02.test.ix<-sqrt(m02.test.roc.obj$sensitivities*m02.test.roc.obj$specificities) %>% which.max()
Gmeans.m02.tcga.ix<-sqrt(m02.tcga.roc.obj$sensitivities*m02.tcga.roc.obj$specificities) %>% which.max()

Gmeans.sig25.test.ix<-sqrt(sig25.test.roc.obj$sensitivities*sig25.test.roc.obj$specificities) %>% which.max()
Gmeans.sig25.tcga.ix<-sqrt(sig25.tcga.roc.obj$sensitivities*sig25.tcga.roc.obj$specificities) %>% which.max()

print(paste('model01 test sens:',m01.test.roc.obj$sensitivities[Gmeans.m01.test.ix],'/','sepc:',m01.test.roc.obj$specificities[Gmeans.m01.test.ix]))
print(paste('model01 tcga sens:',m01.tcga.roc.obj$sensitivities[Gmeans.m01.tcga.ix],'/','sepc:',m01.tcga.roc.obj$specificities[Gmeans.m01.tcga.ix]))

print(paste('model02 test sens:',m02.test.roc.obj$sensitivities[Gmeans.m02.test.ix],'/','sepc:',m02.test.roc.obj$specificities[Gmeans.m02.test.ix]))
print(paste('model02 tcga sens:',m02.tcga.roc.obj$sensitivities[Gmeans.m02.tcga.ix],'/','sepc:',m02.tcga.roc.obj$specificities[Gmeans.m02.tcga.ix]))

print(paste('Sig25 test sens:',sig25.test.roc.obj$sensitivities[Gmeans.sig25.test.ix],'/','sepc:',sig25.test.roc.obj$specificities[Gmeans.sig25.test.ix]))
print(paste('Sig25 tcga sens:',sig25.tcga.roc.obj$sensitivities[Gmeans.sig25.tcga.ix],'/','sepc:',sig25.tcga.roc.obj$specificities[Gmeans.sig25.tcga.ix]))

#get the PR thresholding
m01.test.pr.f.ix<-((2 * m01.test.pr.0$precision* m01.test.pr.0$recall) / (m01.test.pr.0$precision + m01.test.pr.0$recall)) %>% which.max() 
m01.tcga.pr.f.ix<-((2 * m01.tcga.pr.0$precision* m01.tcga.pr.0$recall) / (m01.tcga.pr.0$precision + m01.tcga.pr.0$recall)) %>% which.max() 
m02.test.pr.f.ix<-((2 * m02.test.pr.0$precision* m02.test.pr.0$recall) / (m02.test.pr.0$precision + m02.test.pr.0$recall)) %>% which.max() 
m02.tcga.pr.f.ix<-((2 * m02.tcga.pr.0$precision* m02.tcga.pr.0$recall) / (m02.tcga.pr.0$precision + m02.tcga.pr.0$recall)) %>% which.max() 
sig25.test.pr.f.ix<-((2 * sig25.test.pr.0$precision* sig25.test.pr.0$recall) / (sig25.test.pr.0$precision + sig25.test.pr.0$recall)) %>% which.max() 
sig25.tcga.pr.f.ix<-((2 * sig25.tcga.pr.0$precision* sig25.tcga.pr.0$recall) / (sig25.tcga.pr.0$precision + sig25.tcga.pr.0$recall)) %>% which.max() 

print(paste('model01 test precision:',m01.test.pr.0$precision[m01.test.pr.f.ix],'/','Recall(sensitivity):',m01.test.pr.0$recall[m01.test.pr.f.ix]))
print(paste('model01 tcga precision:',m01.tcga.pr.0$precision[m01.tcga.pr.f.ix],'/','Recall(sensitivity):',m01.tcga.pr.0$recall[m01.tcga.pr.f.ix]))

print(paste('model02 test precision:',m02.test.pr.0$precision[m02.test.pr.f.ix],'/','Recall(sensitivity):',m02.test.pr.0$recall[m02.test.pr.f.ix]))
print(paste('model02 tcga precision:',m02.tcga.pr.0$precision[m02.tcga.pr.f.ix],'/','Recall(sensitivity):',m02.tcga.pr.0$recall[m02.tcga.pr.f.ix]))

print(paste('sig25 test precision:',sig25.test.pr.0$precision[sig25.test.pr.f.ix],'/','Recall(sensitivity):',sig25.test.pr.0$recall[sig25.test.pr.f.ix]))
print(paste('sig25 tcga precision:',sig25.tcga.pr.0$precision[sig25.tcga.pr.f.ix],'/','Recall(sensitivity):',sig25.tcga.pr.0$recall[sig25.tcga.pr.f.ix]))


sink('manuscript_model_perfm_threshold.csv')

print(paste('model01 test sens:',round(m01.test.roc.obj$sensitivities[Gmeans.m01.test.ix],3),'/','sepc:',round(m01.test.roc.obj$specificities[Gmeans.m01.test.ix],3)))
print(paste('model01 tcga sens:',round(m01.tcga.roc.obj$sensitivities[Gmeans.m01.tcga.ix],3),'/','sepc:',round(m01.tcga.roc.obj$specificities[Gmeans.m01.tcga.ix],3)))

print(paste('model02 test sens:',round(m02.test.roc.obj$sensitivities[Gmeans.m02.test.ix],3),'/','sepc:',round(m02.test.roc.obj$specificities[Gmeans.m02.test.ix],3)))
print(paste('model02 tcga sens:',round(m02.tcga.roc.obj$sensitivities[Gmeans.m02.tcga.ix],3),'/','sepc:',round(m02.tcga.roc.obj$specificities[Gmeans.m02.tcga.ix],3)))

print(paste('Sig25 test sens:',round(sig25.test.roc.obj$sensitivities[Gmeans.sig25.test.ix],3),'/','sepc:',round(sig25.test.roc.obj$specificities[Gmeans.sig25.test.ix],3)))
print(paste('Sig25 tcga sens:',round(sig25.tcga.roc.obj$sensitivities[Gmeans.sig25.tcga.ix],3),'/','sepc:',round(sig25.tcga.roc.obj$specificities[Gmeans.sig25.tcga.ix],3)))


print(paste('model01 test precision:',round(m01.test.pr.0$precision[m01.test.pr.f.ix],3),'/','Recall(sensitivity):',round(m01.test.pr.0$recall[m01.test.pr.f.ix],3)))
print(paste('model01 tcga precision:',round(m01.tcga.pr.0$precision[m01.tcga.pr.f.ix],3),'/','Recall(sensitivity):',round(m01.tcga.pr.0$recall[m01.tcga.pr.f.ix],3)))

print(paste('model02 test precision:',round(m02.test.pr.0$precision[m02.test.pr.f.ix],3),'/','Recall(sensitivity):',round(m02.test.pr.0$recall[m02.test.pr.f.ix],3)))
print(paste('model02 tcga precision:',round(m02.tcga.pr.0$precision[m02.tcga.pr.f.ix],3),'/','Recall(sensitivity):',round(m02.tcga.pr.0$recall[m02.tcga.pr.f.ix],3)))

print(paste('sig25 test precision:',round(sig25.test.pr.0$precision[sig25.test.pr.f.ix],3),'/','Recall(sensitivity):',round(sig25.test.pr.0$recall[sig25.test.pr.f.ix],3)))
print(paste('sig25 tcga precision:',round(sig25.tcga.pr.0$precision[sig25.tcga.pr.f.ix],3),'/','Recall(sensitivity):',round(sig25.tcga.pr.0$recall[sig25.tcga.pr.f.ix],3)))

sink()



library(pROC)
m01.test.auroc<-ggroc(m01.test.roc.obj, colour = 'steelblue', size = 2,legacy.axes = T) +
  ggtitle(paste0('ROC Curve ', '(AUC = ', round(m01.test.roc.obj$auc,3), ')')) +
  theme_minimal()+theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                     axis.text=element_text(size=15))+geom_abline(intercept =0 , slope = 1,linetype=2)+coord_equal()

m01.test.auprc<-ggplot(data.frame(m01.test.pr.0$curve),aes(x=X1,y=X2)) + 
  geom_path(colour = 'steelblue',size=2) +
  labs(x="Recall",y="Precision",
       title=format(m01.test.pr.0$auc.integral,digits=3),
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
