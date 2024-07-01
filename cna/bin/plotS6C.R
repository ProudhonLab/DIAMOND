#!/usr/bin/env Rscript
# Script by Klaus von Grafenstein, 24_01_26

# @@@@@@@@@@@@@ LOADING LIBRAIRIES @@@@@@@@@@@@@
library(tidyverse)
library(grid)
library(ggnewscale)

# @@@@@@@@@@@@@ FUNCTIONS @@@@@@@@@@@@@


plotClasséesGzscore <- function( versusvalue, interval, dt, x) {
  print(versusvalue)
  print(ProbaCancer99Sep[which(ProbaCancer99Sep==versusvalue),2])
  if (any(dt$type == "validation_cohort"& dt$versus == versusvalue)) {
    df <- dt %>% filter(type == "validation_cohort" & versus == versusvalue)
    head(df)
    titlename <- unique(df$Title)
    filename <- paste0("figures/S6/",unique(df$file))
    ggplot(df, aes(x = df[[x]], y = df[["Gzscore"]], color = df[["color"]])) +
      
      labs(
        title = "",
      ) +
      xlab(x) +  geom_hline(yintercept = 121)+geom_vline(xintercept = ProbaCancer99Sep[which(ProbaCancer99Sep==versusvalue),2]) +geom_point(size=4)+
      ylab(paste0("Genome Wide ", " z-score")) +
      scale_color_identity("Biological class", labels = str_replace_all(unique(df[["biological_class"]]), "_", " "), breaks = unique(df[["color"]]), guide = "legend") +theme_light()+
      theme(legend.position = "none",text = element_text(size=14))+xlim(0,1)
    ggsave(filename)
    print("yup")
  } else {
    print("nope")
  }
}






# @@@@@@@@@@@@@ MAIN @@@@@@@@@@@@@


S3paper <- read.csv("S3paper.csv", sep=";")
GzscoreArm <-
  read.csv(
    "Zscore/All_Gzscore.csv",header = FALSE,
    
    sep = ";") %>% dplyr::rename(c("Gzscore"="V1","sample"="V2"))%>%
  left_join(S3paper, by=c("sample"="Sample_ID"))%>%filter(sample %in% S3paper$Sample_ID)


color_table_c2 <- read.csv("../../projects/l1pa_meth/5_classifications/color_table_c2_mgorse.csv", sep = ",")
color_table_c1 <- read.csv("../../projects/l1pa_meth/5_classifications/color_table_c1_mgorse.csv", sep = ",")


PredictionsScores <- read.csv("PredictionsScoresHaploClassif.csv", sep = ",", header = FALSE)
names(PredictionsScores) <- c("sample", "biological_class", "n", "probaHD", "predHD", "probaCancer", "predCancer", "classification")

PredictionsScores <- PredictionsScores %>%
  cbind(str_split_fixed(.$classification, pattern = "\\.", n = 5)) %>%
  rename(type = "2", versus = "3") %>%
  select(-c(classification, "1", "4", "5")) %>%
  mutate(Title = paste(str_to_title(str_extract(type, ".*(?=_)")), ":", str_replace_all(versus, "_", " ")), file = paste0( "99speTh.",versus, ".png"),file2 = paste0(type, ".", versus, ".html")) %>%
  left_join(color_table_c1 %>% mutate(type = "discovery_cohort") %>% rbind(color_table_c2 %>% mutate(type = "validation_cohort")))





ProbaCancer99Sep <- read.csv("ProbaCancer99SepHaploClassif.csv", sep="")



frame <- PredictionsScores %>%
  left_join(GzscoreArm %>% select(sample, Gzscore)) %>%
  filter(!is.na(Gzscore))

notAll=frame%>%filter(!versus=="healthy_vs_cancer")
All=frame%>%filter(versus=="healthy_vs_cancer")%>%select(-color,-biological_class)%>%left_join(notAll%>%select(sample,color,biological_class)%>%distinct())

walk(unique(notAll$versus), ~ plotClasséesGzscore(.x,  "Arm", notAll, "probaCancer"))
walk("healthy_vs_cancer", ~ plotClasséesGzscore(.x,  "Arm", All, "probaCancer"))