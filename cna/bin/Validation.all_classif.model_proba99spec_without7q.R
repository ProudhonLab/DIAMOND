#!/usr/bin/env Rscript
# Script by Klaus von Grafenstein, 23_12_15

#@@@@@@@@@@@@@ LOADING LIBRAIRIES @@@@@@@@@@@@@
library(tidyverse)
library(caret)
library("xlsx")

#@@@@@@@@@@@@@ FUNCTIONS @@@@@@@@@@@@@

prepGzscore <- function(x) {
  names(x) <- c("amp", "Gzscore", "sample")
  return(x %>% filter(amp == "Mixed"))
}

Model = function(classif,CnaInfo,df){
  print(paste0("Classification: ",classif))
  print(CnaInfo)
  
  df = df %>%filter(versus==classif)
  cancer= (df%>%filter( biological_class != "healthy_plasma")%>%select(biological_class))[1,1]
  probath = ProbaCancerMaxSenAtMaxSep[which(ProbaCancerMaxSenAtMaxSep$classif==classif),2]
  
  # Predict with predCancer (and sometimes the Gzscore cutoff)
  df=df%>%mutate(predict_biological_class=ifelse(probaCancer <= probath,"healthy_plasma",cancer))
  if(CnaInfo=="WithCna"){
    
    df=df%>%mutate(predict_biological_class_with_CNA=ifelse(predict_biological_class=="healthy_plasma"& Gzscore <= cutoff,"healthy_plasma",cancer),probaCancer_Theshold=probath)
    
    
    #write.xlsx2(df%>%select(sample,biological_class,probaCancer,probaCancer_Theshold,predict_biological_class,Gzscore,predict_biological_class_with_CNA),file = paste0("figures/tablePaperProba99_without7q/validation/",classif,".xlsx"),row.names=FALSE,append=FALSE)
   # write.xlsx2(df%>%select(sample,biological_class,probaCancer,probaCancer_Theshold,predCancer,predict_biological_class,Gzscore,predict_biological_class_with_CNA),file = "figures/tablePaperProba99_without7q/validation/AllTables.xlsx",row.names=FALSE,append=TRUE,sheetName = classif)
    
    
    
    df=df%>%mutate(predict_biological_class=predict_biological_class_with_CNA)
  }
  
  confusionM=caret::confusionMatrix(as.factor(df$predict_biological_class),as.factor( df$biological_class),positive=cancer) # confusion Matrix of the model
  stat=as.matrix(confusionM, what = "classes")
  rownames(stat)=paste(rownames(stat),":")
  
  # If model dir do not exists, create it
  path=paste("figures","ModelProbaWithout7q",CnaInfo,classif,sep = "/")
  if(!dir.exists(path)){dir.create(path)}
  
  # Write results
  write.table(df%>%select(sample,biological_class,predict_biological_class),paste0(path,"/PredictionsTable.csv"),row.names = FALSE,quote = FALSE)
  write.table(confusionM$table,paste0(path,"/ConfusionMatrix.csv"),quote = FALSE,col.names = TRUE)
  write.table( stat,paste0(path,"/statistics.csv"),quote = FALSE,col.names = FALSE,sep=",")
  write.table( round(stat,digits = 3),paste0(path,"/Rounded_statistics.csv"),quote = FALSE,col.names = FALSE,sep=",")
  
  return(as.data.frame(t(c(classif,CnaInfo,stat[1:2],probath))))
  
}

#@@@@@@@@@@@@@ CONSTANTS @@@@@@@@@@@@@
cutoff=121

#@@@@@@@@@@@@@ LOADING FILES @@@@@@@@@@@@@

#@@@ Haplotype Classification @@@
PredictionsScores <- read.csv("PredictionsScoresHaploClassif.csv", sep = ",", header = FALSE)

#@@@ Gzscore @@@

GzscoreArm <-
  read.csv("Zscore/All_Gzscore.csv", header = FALSE,
           
           sep = ";") %>% dplyr::rename(c("Gzscore" = "V1", "sample" = "V2"))

ProbaCancerMaxSenAtMaxSep <- read.csv("ProbaCancer99SepHaploClassif.csv", sep="")

#@@@@@@@@@@@@@ PREPARING DATA @@@@@@@@@@@@@
names(PredictionsScores) <- c("sample", "biological_class", "n", "probaHD", "predHD", "probaCancer", "predCancer", "classification")
PredictionsScores <- PredictionsScores %>%
  cbind(str_split_fixed(.$classification, pattern = "\\.", n = 5)) %>%
  rename(type = "2", versus = "3") %>%
  select(-c(classification, "1", "4", "5")) %>%
  filter(type == "validation_cohort")

names(ProbaCancerMaxSenAtMaxSep)=c("classif","probaTh")


PredictionsScores=PredictionsScores%>%left_join(GzscoreArm) 

#@@@@@@@@@@@@@ MAIN FILE @@@@@@@@@@@@@
index = expand_grid(unique(unique(PredictionsScores$versus)[unique(PredictionsScores$versus) %in% ProbaCancerMaxSenAtMaxSep$classif]),c("WithoutCna","WithCna")) # create a grid to loop trough the different classifications and to predict with and without CNA cutoff
classif_perf_all=map2_dfr(unlist(index[,1]),unlist(index[,2]), ~ Model(.x,.y,PredictionsScores)) # walk to the previous generated grid

names(classif_perf_all)=c("classification_validation",'Type',"sensitivity","specificity","thresholdProbaCancer")
write.table(classif_perf_all,"figures/perfProba99Cancer.csv",row.names = FALSE)
