library(tidyverse)
library(caret)

set.seed(42)


GzscoreArm <-
  read.csv("Zscore/All_Gzscore.csv", header = FALSE,
           
           sep = ";") %>% dplyr::rename(c("Gzscore" = "V1", "sample" = "V2"))

discovery_cohort.healthy_vs_cancer <-
  read.csv(
    "../../projects/l1pa_meth_v2/4_classification/24_05_21_classification_paper_rebutal_cancer_discovery/inputs/discovery/discovery_cohort.healthy_vs_cancer.haplotypes.no_l1hs15_correct.csv")


discovery_cohort.healthy_vs_cancer = discovery_cohort.healthy_vs_cancer %>%
  select(sample, biological_class)



GzscoreArm =  GzscoreArm %>% left_join(discovery_cohort.healthy_vs_cancer) %>%
  filter(sample %in% discovery_cohort.healthy_vs_cancer$sample)



TestThreshold = function(th, df) {
  df = df %>% mutate(Predicted_biological_class = ifelse(Gzscore < th, "healthy_plasma", "cancer_plasma"))
  Sensitivity_Specificity = confusionMatrix(
    as.factor(df$Predicted_biological_class),
    as.factor(df$biological_class),
    positive = "cancer_plasma"
  )$byClass[1:2]
  tab <-
    caret::confusionMatrix(
      as.factor(df$Predicted_biological_class),
      as.factor(df$biological_class),
      positive = "cancer_plasma"
    )
  res <-
    c(tab$byClass, tab$overall[c("Accuracy", "Kappa")]) #code from caret thresholder function
  res <- c(res,
           res["Sensitivity"] + res["Specificity"] - 1,
           sqrt((res["Sensitivity"] - 1) ^ 2 + (res["Specificity"] - 1) ^ 2))
  
  names(res)[-seq_len(length(res) - 2)] <- c("J", "Dist")
  names(res) = str_replace_all(names(res), pattern = " ", replacement = "_")
  rs = as.data.frame(t(c(th, res)))
  colnames(rs) = c("Threshold", names(res))
  return(rs)
}
ForAllFold = function(foldnumber, folds, df) {
  foldindex = folds[foldnumber]
  dfTraining = df[unlist(foldindex), ]
  dfTesting =  df[-unlist(foldindex), ]
  Thresholds = seq(0, round(max(dfTraining$Gzscore)), by = 0.5)
  ROC_values = map_dfr(Thresholds,  ~ TestThreshold(th = .x, df = dfTraining))
  
  th_Max_se_at_max_spe = ROC_values[(which.max(ROC_values[which(ROC_values$Specificity ==
                                                                  max(ROC_values$Specificity)), "Sensitivity"]) - 1) + which.max(ROC_values$Specificity), 1]
  th_99_spec = ROC_values[which.min(abs(0.99 - ROC_values$Specificity)), 1]
  others_max_th = unlist(map(colnames(ROC_values)[4:16], function(stat) {
    th_selected = ROC_values[which.max(ROC_values[, stat]), 1]
    names(th_selected) = paste0("th_", stat)
    ggplot(ROC_values[, c(1, which(names(ROC_values) == stat))], aes_string(y =
                                                                              stat, x = Thresholds)) + geom_point(color = "gray52") + xlab("Gzscore Threshold") +
      geom_line(color = "gray") + labs(
        title = paste0("Training ROC Fold #", foldnumber),
        subtitle = paste0("Theshold at max ", stat, " : ", th_selected)
      )
    ggsave(
      paste0(
        "figures/DetermingCutoffDiscovery7q/",
        stat,
        "/Fold",
        foldnumber,
        ".png"
      )
    )
    return(th_selected)
  }))
  
  
  ggplot(ROC_values, aes(1 - Specificity, Sensitivity)) + geom_point() +
    geom_abline(slope = 1, color = "gray") + geom_line() + labs(
      title = paste0("Training ROC Fold #", foldnumber),
      subtitle = paste0(
        "Max Se at 100% Sp: ",
        ROC_values[which(ROC_values$Threshold == th_Max_se_at_max_spe), 2],
        "\nThreshold: ",
        th_Max_se_at_max_spe
      )
    )
  ggsave(
    paste0(
      "figures/DetermingCutoffDiscovery7q/Training_ROC_Fold#",
      foldnumber,
      ".png"
    )
  )
  
  
  rs_th = as.data.frame(t(
    c(foldnumber, th_Max_se_at_max_spe, th_99_spec, others_max_th)
  ))
  colnames(rs_th) = c("fold",
                      "th_max_Se_at_max_Spe",
                      "th_99_spec",
                      names(others_max_th))
  return(rs_th)
}

folds = groupKFold(GzscoreArm$sample, 5)
final_df = map_dfr(1:5,  ~ ForAllFold(.x, folds, GzscoreArm))
write.table(
  final_df,
  "figures/DetermingCutoffDiscovery7q/thresholdByFoldAndMaxStats.csv",
  row.names = FALSE
)