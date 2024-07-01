library(tidyverse)
library(ggpubr)




palette <-
  c(
    "CRC M+" = "#FEB700",
    "BRC M+" = "#0060CB",
    "LC M+" = "#DDB4FC",
    "OVC M+" = "#BF5858",
    "All"="grey",
    "BRC M0"= "#78C8FE",
    "OVC M0" = "#FEABAE"
  )


perfProba99Cancer <-
  read.csv("figures/perfProba99Cancer.csv",
           sep = "")

perfProba99Cancer=perfProba99Cancer %>% mutate(
  classification_validation = case_when(
    classification_validation == "healthy_vs_cancer" ~ "All",
    
    classification_validation == "healthy_vs_ovarian_cancer_plasma" ~ "OVC M+",
    classification_validation == "healthy_vs_early_ovarian_cancer_plasma" ~ "OVC M0",
    classification_validation == "healthy_vs_early_breast_cancer_plasma" ~ "BRC M0",
    classification_validation == "healthy_vs_breast_cancer_plasma" ~ "BRC M+",
    classification_validation == "healthy_vs_colorectal_cancer_plasma" ~ "CRC M+",
    classification_validation == "healthy_vs_lung_cancer_plasma" ~ "LC M+",TRUE~NA
    
  )
)%>%filter((!is.na(classification_validation)))



palette3= c(palette,rep("white",length(palette)))
names(palette3)=c(paste0(names(palette),"WithCna"),paste0(names(palette),"WithoutCna"))


palette4= rep(palette,2)
names(palette4)=names(palette3)

ggplot(
  perfProba99Cancer%>%mutate(Type=factor(Type,c("WithoutCna","WithCna")))%>%arrange(classification_validation,Type),
  aes(
    x = classification_validation,
    y = sensitivity,
    fill = paste0(classification_validation, Type),
    group = Type
  )
) + geom_bar(
  aes(color = paste0(classification_validation, Type)),
  stat = "identity",
  width = 0.7,
  position = position_dodge(0.7),
  linewidth = 1.5,
) +geom_text(
  aes(label=paste0( round((sensitivity)*100),"%"), y = sensitivity + 0.01),position = position_dodge(0.7), vjust=0,color="grey10",fontface="bold"
)+ scale_fill_manual(values = palette3) + scale_color_manual(values = palette4) +  theme_pubr()+
  theme(legend.position = "none",text = element_text(size=16),axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) + scale_y_continuous(
    labels = scales::percent,
    expand = c(0, 0),
    limits = c(0, 1.05)
  ) +
  labs(x = "", y = "Sensitivity (%)") + geom_text(
    aes(label=paste0( "Spec=",round((specificity)*100),"%"),y=0.01,),color="grey10", angle=90,hjust="left",fontface="bold", size=4, position = position_dodge(0.7)
    
  ) 

ggsave("figures/5G_Dual.png", width = 8.5, height = 7)


ggplot(
  perfProba99Cancer%>%filter(Type=="WithCna")%>%arrange(classification_validation),
  aes(
    x = classification_validation,
    y = sensitivity,
    fill = classification_validation,
    group = Type
  )
) + geom_bar(
 
  stat = "identity",
  width = 0.7,
  position = position_dodge(0.7),
  linewidth = 1.5,
) +geom_text(
  aes(label=paste0( round((sensitivity)*100),"%"), y = sensitivity + 0.01),position = position_dodge(0.7), vjust=0,color="grey10",fontface="bold"
)+ scale_fill_manual(values = palette) +  theme_pubr()+
  theme(legend.position = "none",text = element_text(size=16),axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  )) + scale_y_continuous(
    labels = scales::percent,
    expand = c(0, 0),
    limits = c(0, 1.05)
  ) +
  labs(x = "", y = "Sensitivity (%)") + geom_text(
    aes(label=paste0( "Spec=",round((specificity)*100),"%"),y=0.01,),color="grey10", angle=90,hjust="left",fontface="bold", size=4, position = position_dodge(0.7)
    
  ) 

ggsave("figures/5G_WithCNA.png", width = 4.5, height = 5)
