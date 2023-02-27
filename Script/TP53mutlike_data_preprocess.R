#2022.11.09
#Final version of the data: GEP (BEAT and TCGA) and Clinical data

packages<-c('edgeR','limma','dplyr','tidyr','sva','tidyverse',
            'ggplot2','survival','survminer','HGNChelper','stringr','ComplexHeatmap')


for(p in packages){
  library(p,character.only = T)
}
setwd("C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/")
BEAT_RNA_counts.461.clinic<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/BEAT_AML/BEAT_AML_preprocessed/BEAT_RNA_counts_461_clinic.rds")
TCGA_AML_counts<-readRDS("C:/Users/yklee/Desktop/Sachs Lab/R_project/RNAseq/TCGA/TCGA_preprocessed/TCGA_RNA_counts.rds")
BEAT_clinical_final.ELN<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/BEAT2022_clinical_final_mtwt_ELN.rds')
TCGA_clinical_final.ELN<-readRDS('C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')
AML_id<-BEAT_clinical_final.ELN$LabId[BEAT_clinical_final.ELN$ELN2017!='NonAML']
BEAT_RNA_counts.AML.mtmtlwt<-BEAT_RNA_counts.461.clinic[,AML_id]
BEAT_clinical_final.ELN.AML<-BEAT_clinical_final.ELN[BEAT_clinical_final.ELN$ELN2017!='NonAML',]


BEAT_genes<-rownames(BEAT_RNA_counts.461.clinic)
TCGA_genes<-rownames(TCGA_AML_counts)


t<-getCurrentHumanMap()
BEAT_HGNC<-checkGeneSymbols(BEAT_genes,map = t)
TCGA_HGNC<-checkGeneSymbols(TCGA_genes,map = t)
#get limma genes for the same genes
BEAT_HGNC$limma<-alias2SymbolTable(BEAT_genes,species='Hs')
TCGA_HGNC$limma<-alias2SymbolTable(TCGA_genes,species='Hs')
BEAT_HGNC$corrected<-BEAT_HGNC$x
TCGA_HGNC$corrected<-TCGA_HGNC$x

F_ix.beat<-which(BEAT_HGNC$Approved==FALSE)
F_ix.tcga<-which(TCGA_HGNC$Approved==FALSE)
for(i in F_ix.beat){
  if (BEAT_HGNC$Approved[i]==FALSE){
    if(is.na(BEAT_HGNC$Suggested.Symbol[i]) &is.na(BEAT_HGNC$limma[i])){#both NA
      BEAT_HGNC$corrected[i]<-BEAT_HGNC$x[i]
    }else if(is.na(BEAT_HGNC$Suggested.Symbol[i]) | is.na(BEAT_HGNC$limma[i])){# either one is NA
      
      if(is.na(BEAT_HGNC$Suggested.Symbol[i])){
        BEAT_HGNC$corrected[i]<-BEAT_HGNC$limma[i]
      }else{
        BEAT_HGNC$corrected[i]<-BEAT_HGNC$Suggested.Symbol[i]
      }
    }else{
      BEAT_HGNC$corrected[i]<-BEAT_HGNC$limma[i]
    }
  }else{#when Approved ==TRUE
    BEAT_HGNC$corrected[i]<-BEAT_HGNC$x[i]
  }
}
which(duplicated(BEAT_HGNC$corrected))
for(i in F_ix.tcga){
  if (TCGA_HGNC$Approved[i]==FALSE){
    if(is.na(TCGA_HGNC$Suggested.Symbol[i]) &is.na(TCGA_HGNC$limma[i])){#both NA
      TCGA_HGNC$corrected[i]<-TCGA_HGNC$x[i]
    }else if(is.na(TCGA_HGNC$Suggested.Symbol[i]) | is.na(TCGA_HGNC$limma[i])){# either one is NA
      
      if(is.na(TCGA_HGNC$Suggested.Symbol[i])){
        TCGA_HGNC$corrected[i]<-TCGA_HGNC$limma[i]
      }else{
        TCGA_HGNC$corrected[i]<-TCGA_HGNC$Suggested.Symbol[i]
      }
    }else{
      TCGA_HGNC$corrected[i]<-TCGA_HGNC$limma[i]
    }
  }else{#when Approved ==TRUE
    TCGA_HGNC$corrected[i]<-TCGA_HGNC$x[i]
  }
}
TCGA_HGNC$corrected[which(duplicated(TCGA_HGNC$corrected))]
#c('a','a','b','b')%>% make.names(.,unique=T)
BEAT_HGNC$corrected<-BEAT_HGNC$corrected%>% make.names(.,unique=T)
TCGA_HGNC$corrected<-TCGA_HGNC$corrected%>% make.names(.,unique=T)
length(intersect(BEAT_HGNC$corrected,TCGA_HGNC$corrected))
rownames(BEAT_RNA_counts.461.clinic)<-BEAT_HGNC$corrected
rownames(TCGA_AML_counts)<-TCGA_HGNC$corrected
saveRDS(BEAT_RNA_counts.AML.mtmtlwt,"C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_counts_mtmtlwt.rds")
saveRDS(TCGA_AML_counts,"C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/TCGA_AML_counts_mtmtlwt.rds")
saveRDS(BEAT_clinical_final.ELN.AML,"C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.rds")
saveRDS(TCGA_clinical_final.ELN,'C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.rds')
write.csv(BEAT_clinical_final.ELN.AML,"C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Final Version/BEAT_AML_403_clinical_final_ELN.csv")
write.csv(TCGA_clinical_final.ELN,'C:/Users/yklee/Desktop/Sachs Lab/Collaboration works/TP53 analysis/Data/Clinical data final 2022/TCGA_clinical_final_ELN.csv')
