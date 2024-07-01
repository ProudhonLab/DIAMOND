##################
# October 9th, 2023 #
##################
# We want to do classifications after adding new HD samples in C1 (30 pools) + HD MHOVC25 in C2 
# This version is made for MAC version (run in GenOuest)

# !!! change this !!!
user = "mgorse"
# !!!

##############################################
#  MHOVC25 : batch 3 IN DISCOVERY, BATCH 1,2 IN VALIDATION #
##############################################

workDir_d = paste0("/scratch/",user,"/20231009_classification_batch3_MHOVC25_C1_and_batch1_2_MHOVC25_C2/discovery/inputs/")
workDir_v = paste0("/scratch/",user,"/20231009_classification_batch3_MHOVC25_C1_and_batch1_2_MHOVC25_C2/validation/inputs/")


# Discovery

c1_infos = read.csv( paste0("/groups/proudhon_lab/projects/l1pa_meth/0_data/c1_plasma_samples_info.csv"))
rownames(c1_infos) = c1_infos$sample
c1_infos$biological_class[which(c1_infos$biological_class=="gastric_cancer_plasma")] = "early_gastric_cancer_plasma"
c1_infos$biological_class[which(c1_infos$biological_class=="ovarian_cancer_plasma" & c1_infos$stage=="3")] = "early_ovarian_cancer_plasma"
c1_infos$biological_class[which(c1_infos$biological_class=="ovarian_cancer_plasma" & is.na(c1_infos$stage))] = "ovarian_cancer_plasma_ND"


c1_methyl = read.csv( paste0("/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/cg_methyl.c1.csv"))
c1_methyl = c1_methyl[,-1]
rownames(c1_methyl) = c1_methyl$sample
c1_methyl = c1_methyl[rownames(c1_infos),-grep("L1HS_15",colnames(c1_methyl))]
c1_methyl$biological_class = c1_infos$biological_class

# Ajout des eUVM de MHOVC25 en discovery et du batch 3
infos_all= read.csv( paste0("/groups/proudhon_lab/projects/l1pa_meth/0_data/samples_infos.csv"))
rownames(infos_all) = infos_all$sample
infos_eUVM = infos_all[which(infos_all$mhovc == "mhovc25" & infos_all$type == "plasma" & infos_all$class != "healthy" & infos_all$paper_cohort == "discovery"),]  
methyl_mhovc25 = read.csv( paste0("/groups/proudhon_lab/projects/l1pa_meth_v2/3_methylation_data_vs_all/cg_methyl.all_samples.csv"))
rownames(methyl_mhovc25) = methyl_mhovc25$sample
c1d_methyl=methyl_mhovc25[which(methyl_mhovc25$sample %in% infos_eUVM$sample ),]
c1d_methyl = c1d_methyl [,-grep("L1HS_15",colnames(c1d_methyl))]

list_batch3 = c('D1525S085', 'D1525S086', 'D1525S087','D1525S088','D1525S089','D1525S090','D1525S091','D1525S092','D1525S093','D1525S094','D1525S095','D1525S096','D1525S097','D1525S098','D1525S099','D1525S100','D1525S101','D1525S102', 'D1525S103')
info_batch3 = infos_all[which(infos_all$sample %in% list_batch3),]  
c1e_methyl=methyl_mhovc25[which(methyl_mhovc25$sample %in% info_batch3$sample ),]
c1e_methyl = c1e_methyl [,-grep("L1HS_15",colnames(c1e_methyl))]


# Merge de toutes les data de discovery
c1_methyl = rbind(c1_methyl, c1d_methyl, c1e_methyl)


# Séparation des datat en HD - cancer
healthy_methyl = c1_methyl[which(c1_methyl$biological_class=="healthy_plasma"),]

methyl_cancer = c1_methyl[which(c1_methyl$biological_class!="healthy_plasma"),]
methyl_cancer$biological_class = "cancer_plasma"
m = rbind(healthy_methyl,methyl_cancer)
write.csv(m,  paste0(workDir_d,"discovery_cohort.healthy_vs_cancer.methylation.no_l1hs15.csv"), row.names=FALSE)

for (cancer in unique(c1_methyl$biological_class)) {
  if (cancer != "healthy_plasma" & cancer != "ovarian_cancer_plasma_ND") {
    df = rbind(healthy_methyl,c1_methyl[which(c1_methyl$biological_class==cancer),])
    cancer_type = gsub(" ","",cancer)
    cancer_type = gsub("\\+","p",cancer_type)
    cancer_type = gsub("0","z",cancer_type)
    write.csv(df,  paste0(workDir_d,"discovery_cohort.healthy_vs_",cancer_type,".methylation.no_l1hs15.csv"), row.names=FALSE)
  }
}

methyl_cancer = c1_methyl[grep("ovarian_cancer_plasma",c1_methyl$biological_class),]
methyl_cancer$biological_class = "ovarian_cancer_plasma"
m = rbind(healthy_methyl,methyl_cancer)
write.csv(m,  paste0(workDir_d,"discovery_cohort.healthy_vs_all_ovarian_cancer_plasma.methylation.no_l1hs15.csv"), row.names=FALSE)

# validation

c2_infos = read.csv( paste0("/groups/proudhon_lab/projects/l1pa_meth/0_data/c2_plasma_samples_info.csv"))
rownames(c2_infos) = c2_infos$sample
c2_infos$biological_class[which(c2_infos$biological_class=="ovarian_cancer_plasma" & c2_infos$stage=="3")] = "early_ovarian_cancer_plasma"

c2_methyl = read.csv( paste0("/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/cg_methyl.c2.csv"))
c2_methyl = c2_methyl[,-1]
rownames(c2_methyl) = c2_methyl$sample
c2_methyl = c2_methyl[rownames(c2_infos),-grep("L1HS_15",colnames(c2_methyl))]
c2_methyl$biological_class = c2_infos$biological_class



c2_infos_mhovc25 = infos_all[which(infos_all$mhovc == "mhovc25" & infos_all$type == "plasma" & infos_all$paper_cohort == "validation"),]

c2d_methyl=methyl_mhovc25[which(methyl_mhovc25$sample %in% c2_infos_mhovc25$sample ),]
c2d_methyl=c2d_methyl[which(!(c2d_methyl$sample %in% list_batch3)),]
c2d_methyl = c2d_methyl[,-grep("L1HS_15",colnames(c2d_methyl))]
c2d_methyl = c2d_methyl[-grep("D1525S028", c2d_methyl$sample),]
c2d_methyl = c2d_methyl[-grep("D1525S020", c2d_methyl$sample),]
c2d_methyl = c2d_methyl[-grep("D1525S021", c2d_methyl$sample),]
c2d_methyl = c2d_methyl[-grep("D1525S038", c2d_methyl$sample),]
c2d_methyl = c2d_methyl[-grep("D1525S039", c2d_methyl$sample),]



c2_methyl = rbind(c2_methyl,c2c_methyl,c2d_methyl)
c2_methyl$biological_class = factor(c2_methyl$biological_class)

healthy_methyl = c2_methyl[which(c2_methyl$biological_class=="healthy_plasma"),]

methyl_cancer = c2_methyl[which(c2_methyl$biological_class!="healthy_plasma"),]
methyl_cancer$biological_class = "cancer_plasma"
m = rbind(healthy_methyl,methyl_cancer)
write.csv(m,  paste0(workDir_v,"validation_cohort.healthy_vs_cancer.methylation.no_l1hs15.csv"), row.names=FALSE)

for (cancer in levels(c2_methyl$biological_class)) {
  if (cancer != "healthy_plasma") {
    df = rbind(healthy_methyl,c2_methyl[which(c2_methyl$biological_class==cancer),])
    cancer_type = gsub(" ","",cancer)
    cancer_type = gsub("\\+","p",cancer_type)
    cancer_type = gsub("0","z",cancer_type)
    write.csv(df,  paste0(workDir_v,"validation_cohort.healthy_vs_",cancer_type,".methylation.no_l1hs15.csv"), row.names=FALSE)
  }
}

methyl_cancer = c2_methyl[grep("ovarian_cancer_plasma",c2_methyl$biological_class),]
methyl_cancer$biological_class = "ovarian_cancer_plasma"
m = rbind(healthy_methyl,methyl_cancer)
write.csv(m,  paste0(workDir_v,"validation_cohort.healthy_vs_all_ovarian_cancer_plasma.methylation.no_l1hs15.csv"), row.names=FALSE)



