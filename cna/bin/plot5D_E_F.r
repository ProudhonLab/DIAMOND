library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggrepel)
library(ggnewscale)
library(viridis)
library(glue)


S3paper <- read.csv("S3paper.csv", sep=";")

cg_methyl.all_samples <-
  read.csv(
    "../../projects/l1pa_meth/4_methylation_data/default_scores/10_largest/cg_methyl.all_samples.final.csv"
  )
cg_methyl.all_samples.newHD <-
  read.csv(
    "../../projects/l1pa_meth_v2/3_methylation_data_vs_all/cg_methyl.all_samples.csv"
  )%>%rbind( read.csv(
    "../../projects/extraction_methods_comparison/3_methylation_data/cg_methyl.all_samples.csv"))

GzscoreArm <-
  read.csv(
    "Zscore/All_Gzscore.csv",header = FALSE,
   
    sep = ";") %>% dplyr::rename(c("Gzscore"="V1","sample"="V2"))%>%
  left_join(S3paper, by=c("sample"="Sample_ID"))

`GzscoreArm` <- `GzscoreArm` %>%filter(sample %in% S3paper$Sample_ID)%>%
  mutate(cancer = if_else(Disease_status != "healthy" & !(is.na(Disease_status)), "Cancer", "Healthy")) %>%
  mutate(
    Disease_status = fct_recode(
      factor(
        Disease_status,
        levels = c(
          "colorectal_cancer",
          "breast_cancer",
          "uveal_melanoma_cancer",
          "lung_cancer",
          "healthy",
          "ovarian_cancer",
          "gastric_cancer"
        )
      ),
      "CRC" = "colorectal_cancer",
      "BRC" = "breast_cancer",
      "UVC" = "uveal_melanoma_cancer",
      "LC" = "lung_cancer",
      "HD" = "healthy",
      "OVC" = "ovarian_cancer",
      "GAC" = "gastric_cancer"
    )
  ) %>%
  group_by(Disease_status) %>%
  mutate(sample_size = n()) %>%
  ungroup()


GlobalMethData <- cg_methyl.all_samples %>% select(         -X,)%>%rbind(cg_methyl.all_samples.newHD)%>%
  select(-starts_with("L1HS_15"),

         -biological_class) %>%
  pivot_longer(!sample, names_to = "pos", values_to = "meth") %>%
  group_by(sample) %>%
  mutate(GlobalMeth = mean(meth)) %>%
  select(GlobalMeth, sample) %>%
  distinct()



`ArmMergedGlobalZscoreMeta`  <-  GzscoreArm  %>% select(-Stage) %>%arrange(Disease_status, Metastasis_status) %>%filter(!(is.na(Metastasis_status) & Disease_status == "OVC"))%>%
  
  mutate(Disease_status = factor(ifelse(
    cancer == "Cancer",
    paste(Disease_status, Metastasis_status),
    "HD"
  ))) %>%
  group_by(Disease_status) %>%
  mutate(sample_size = n()) %>%ungroup()





pwcArmMeta <-
  wilcox_test(
    ArmMergedGlobalZscoreMeta,
    Gzscore ~ Disease_status,
    ref.group = "HD",
    p.adjust.method = "none"
  )
pArmMeta <-
  ggboxplot(
    ArmMergedGlobalZscoreMeta,
    x = "Disease_status",
    y = "Gzscore",
    col = "Disease_status",
    add = "jitter",
    order = c("CRC M+",  "BRC M+", "UVC M+", "LC M+","OVC M+", "GAC M+","HD","BRC M0","OVC M0","GAC M0"),
    palette = c(
      "#FEB700",
      "#0060CB",
      
      "#ff9bdd",
      "#DDB4FC","#BF5858",  "#008847",
      "grey","#78C8FE",
      
      "#FEABAE",
      
      
      "#02c065"
    )
    
  ) + geom_segment(aes(
    x = 6.5,
    xend = 6.5,
    y = 0,
    yend = 3000
  ),
  linetype = "dashed",
  color = "gray") + geom_segment(aes(
    x = 7.5,
    xend = 7.5,
    y = 0,
    yend = 3000
  ),
  linetype = "dashed",
  color = "gray")+
  stat_pvalue_manual(
    pwcArmMeta %>% mutate(group2 = factor(
      group2, c("CRC M+",  "BRC M+", "UVC M+", "LC M+","OVC M+", "GAC M+","HD","BRC M0","OVC M0","GAC M0")
    )) %>% arrange(group2),
    x = "group2", y.position = 4000,
    label = "p.adj.signif",angle=90
  ) +
  theme(legend.position = "none",
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1
        )) + labs(
          title = "",
         
          y = "Genome-wide z-score",
          x = ""
        )

ggsave(
  "figures/Fig_5E_without7q.png",
  pArmMeta,
  height = 6,
  width = 4.2
)
write.table(pwcArmMeta, "figures/tableStat/Fig5E_withot7q.csv", row.names = FALSE)



palette2 = c(
  "CRC M+" = "#FEB700",
  "BRC M+" = "#0060CB",
  "BRC M0" = "#78C8FE",
  "UVC M+" = "#ff9bdd",
  "LC M+" = "#DDB4FC",
  "HD" = "grey",
  "OVC M+" = "#BF5858",
  "OVC M0" = "#FEABAE",
  "OVC NA" = "#FEDDDE",
  "GAC M+" = "#008847",
  "GAC M0" = "#02c065"
)


ggplot(
  `ArmMergedGlobalZscoreMeta`
  %>% left_join(GlobalMethData),
  aes(Gzscore, GlobalMeth, col = Disease_status)
) +
  stat_smooth(method = "lm", se = FALSE) +
  geom_point() +
  labs(title = "",
       y = "Global methylation",
       x = "Genome-wide z-score") +
  scale_colour_manual(values = palette2) + theme_light() +theme(legend.position = "none")

ggsave("figures/Fig_5F_without7q.png",height = 4,width = 4)

pwcCancerArm=  wilcox_test(
  GzscoreArm,
  Gzscore ~ cancer,
)
  
  ggboxplot(
    `GzscoreArm` %>% group_by(cancer) %>% mutate(sample_size = n()),
    x = "cancer",
    y = "Gzscore",
    add = "jitter"
  ) +
  stat_pvalue_manual(pwcCancerArm, y.position = c(4000)) +
  theme(legend.position = "none") + labs(
    title = "",
  
    y = "Genome-wide z-score",
    x = ""
  ) +  geom_text(aes(x = cancer,
                     y = -100,
                     label = sample_size),
                 check_overlap = TRUE)

ggsave("figures/Fig5D_without7q.png",units = "cm",width = 10.6, height = 14)
write.table(pwcCancerArm, "figures/tableStat/Fig5D_without7q.csv", row.names = FALSE)


