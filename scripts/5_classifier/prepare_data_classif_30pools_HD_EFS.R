##################
# October 3th, 2023 #
##################
# We want to do classifications after addin new HD samples in C2 (30 pools) , plus MHOVC25 cancers and EFS HD
# This version is made for genouest running

# !!! change this !!!
user = "kvongrafenst"
# !!!


workDir_d = paste0("/scratch/",user,"/20231004_classification_with_30pools_in_C2_EFS_samples/discovery/inputs/")
workDir_v = paste0("/scratch/",user,"/20231004_classification_with_30pools_in_C2_EFS_samples/validation/inputs/")


# Discovery

c1_infos = read.csv(paste0("/groups/proudhon_lab/projects/l1pa_meth/0_data/c1_plasma_samples_info.csv"))
rownames(c1_infos) = c1_infos$sample
c1_infos$biological_class[which(c1_infos$biological_class=="gastric_cancer_plasma")] = "early_gastric_cancer_plasma"
c1_infos$biological_class[which(c1_infos$biological_class=="ovarian_cancer_plasma" & c1_infos$stage=="3")] = "early_ovarian_cancer_plasma"
c1_infos$biological_class[which(c1_infos$biological_class=="ovarian_cancer_plasma" & is.na(c1_infos$stage))] = "ovarian_cancer_plasma_ND"

c1_methyl = read.csv(paste0("/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/cg_methyl.c1.csv"))
c1_methyl = c1_methyl[,-1]
rownames(c1_methyl) = c1_methyl$sample
c1_methyl = c1_methyl[rownames(c1_infos),-grep("L1HS_15",colnames(c1_methyl))]
c1_methyl$biological_class = c1_infos$biological_class

infos_all= read.csv(paste0("/groups/proudhon_lab/projects/l1pa_meth/0_data/samples_infos.csv"))

healthy_methyl = c1_methyl[which(c1_methyl$biological_class=="healthy_plasma"),]

methyl_cancer = c1_methyl[which(c1_methyl$biological_class!="healthy_plasma"),]
methyl_cancer$biological_class = "cancer_plasma"


rownames(infos_all) = infos_all$sample
c1_infos_mhovc25 = infos_all[which(infos_all$mhovc == "mhovc25" & infos_all$type == "plasma" & infos_all$class != "healthy" & infos_all$paper_cohort == "discovery"),]

methyl_mhovc25 = read.csv(paste0("/groups/proudhon_lab/projects/l1pa_meth_v2/3_methylation_data_vs_all/cg_methyl.all_samples.csv"))
rownames(methyl_mhovc25) = methyl_mhovc25$sample

c1d_methyl=methyl_mhovc25[which(methyl_mhovc25$sample %in% c1_infos_mhovc25$sample ),]

c1d_methyl = c1d_methyl[,-grep("L1HS_15",colnames(c1d_methyl))]

mhovc25methy = c1d_methyl
mhovc25methy$biological_class = "cancer_plasma"


m = rbind(healthy_methyl,methyl_cancer,mhovc25methy)
write.csv(m, paste0(workDir_d,"discovery_cohort.healthy_vs_cancer.methylation.no_l1hs15.csv"), row.names=FALSE)

c1_methyl = rbind(c1_methyl,c1d_methyl)
for (cancer in unique(c1_methyl$biological_class)) {
  if (cancer != "healthy_plasma" & cancer != "ovarian_cancer_plasma_ND") {
    df = rbind(healthy_methyl,c1_methyl[which(c1_methyl$biological_class==cancer),])
    cancer_type = gsub(" ","",cancer)
    cancer_type = gsub("\\+","p",cancer_type)
    cancer_type = gsub("0","z",cancer_type)
    write.csv(df, paste0(workDir_d,"discovery_cohort.healthy_vs_",cancer_type,".methylation.no_l1hs15.csv"), row.names=FALSE)
  }
}

methyl_cancer = c1_methyl[grep("ovarian_cancer_plasma",c1_methyl$biological_class),]
methyl_cancer$biological_class = "ovarian_cancer_plasma"
m = rbind(healthy_methyl,methyl_cancer)
write.csv(m, paste0(workDir_d,"discovery_cohort.healthy_vs_all_ovarian_cancer_plasma.methylation.no_l1hs15.csv"), row.names=FALSE)

# validation

c2_infos = read.csv(paste0("/groups/proudhon_lab/projects/l1pa_meth/0_data/c2_plasma_samples_info.csv"))
rownames(c2_infos) = c2_infos$sample
c2_infos$biological_class[which(c2_infos$biological_class=="ovarian_cancer_plasma" & c2_infos$stage=="3")] = "early_ovarian_cancer_plasma"

c2_methyl = read.csv(paste0("/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/cg_methyl.c2.csv"))
c2_methyl = c2_methyl[,-1]
rownames(c2_methyl) = c2_methyl$sample
c2_methyl = c2_methyl[rownames(c2_infos),-grep("L1HS_15",colnames(c2_methyl))]
c2_methyl$biological_class = c2_infos$biological_class

methyl = read.csv(paste0("/groups/proudhon_lab/projects/extraction_methods_comparison/3_methylation_data/cg_methyl.all_samples.csv"))
infos = read.table(paste0("/groups/proudhon_lab/projects/extraction_methods_comparison/0_data/sequencing_data/D1357.sampleDescription.txt"),sep="|")
methyl = methyl[,-grep("L1HS_15",colnames(methyl))]
rownames(methyl) = methyl$sample
methyl$name = infos$V2[match(methyl$sample,infos$V1)]
methyl$individual = gsub("(.*)_.*","\\1",methyl$name)
methyl$method = gsub(".*_(.*)","\\1",methyl$name)
methyl = methyl[-which(methyl$method=="P" | methyl$method=="100" | methyl$method=="50" | methyl$method=="0"),]
methyl$biological_class = "healthy_plasma"
c2c_methyl = methyl[which(methyl$method=="Mn"),1:32]

rownames(infos_all) = infos_all$sample
c2_infos_mhovc25 = infos_all[which((infos_all$mhovc == "mhovc25" & infos_all$type == "plasma" & infos_all$patient_cohort == "EFS") | (infos_all$mhovc == "mhovc25" & infos_all$type == "plasma" & infos_all$class != "healthy" & infos_all$paper_cohort == "validation")),]
c2_infos_mhovc25 = c2_infos_mhovc25[which(c2_infos_mhovc25 != "D1525S028"),]

methyl_mhovc25_c2 = read.csv(paste0("/groups/proudhon_lab/projects/l1pa_meth_v2/3_methylation_data_vs_all/cg_methyl.all_samples.csv"))
rownames(methyl_mhovc25_c2) = methyl_mhovc25_c2$sample

c2d_methyl=methyl_mhovc25_c2[which(methyl_mhovc25_c2$sample %in% c2_infos_mhovc25$sample ),]

c2d_methyl = c2d_methyl[,-grep("L1HS_15",colnames(c2d_methyl))]



c2_methyl = rbind(c2_methyl,c2c_methyl,c2d_methyl)
c2_methyl$biological_class = factor(c2_methyl$biological_class)

healthy_methyl = c2_methyl[which(c2_methyl$biological_class=="healthy_plasma"),]

methyl_cancer = c2_methyl[which(c2_methyl$biological_class!="healthy_plasma"),]
methyl_cancer$biological_class = "cancer_plasma"
m = rbind(healthy_methyl,methyl_cancer)
write.csv(m, paste0(workDir_v,"validation_cohort.healthy_vs_cancer.methylation.no_l1hs15.csv"), row.names=FALSE)

for (cancer in levels(c2_methyl$biological_class)) {
  if (cancer != "healthy_plasma") {
    df = rbind(healthy_methyl,c2_methyl[which(c2_methyl$biological_class==cancer),])
    cancer_type = gsub(" ","",cancer)
    cancer_type = gsub("\\+","p",cancer_type)
    cancer_type = gsub("0","z",cancer_type)
    write.csv(df, paste0(workDir_v,"validation_cohort.healthy_vs_",cancer_type,".methylation.no_l1hs15.csv"), row.names=FALSE)
  }
}

methyl_cancer = c2_methyl[grep("ovarian_cancer_plasma",c2_methyl$biological_class),]
methyl_cancer$biological_class = "ovarian_cancer_plasma"
m = rbind(healthy_methyl,methyl_cancer)
write.csv(m, paste0(workDir_v,"validation_cohort.healthy_vs_all_ovarian_cancer_plasma.methylation.no_l1hs15.csv"), row.names=FALSE)
