# !!! change this !!!
user = "kdasilva"
inputDir = paste0("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=",user,"/groups/proudhon_lab/projects/l1pa_meth/")
workDir_d = paste0("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=",user,"/scratch/kdasilva/20230623_trainingMz_testingMp/inputs/discovery/")
workDir_v = paste0("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=",user,"/scratch/kdasilva/20230623_trainingMz_testingMp/inputs/validation/")

#############
# Discovery #
#############

# Load data

# curated list C1
c1 = read.csv(paste0(inputDir,"0_data/c1_plasma_samples_info.csv"))
rownames(c1) = c1$sample
c1$biological_class[which(c1$biological_class=="gastric_cancer_plasma")] = "early_gastric_cancer_plasma"
c1$biological_class[which(c1$biological_class=="ovarian_cancer_plasma" & c1$stage=="3")] = "early_ovarian_cancer_plasma"
c1$biological_class[which(c1$biological_class=="ovarian_cancer_plasma" & is.na(c1$stage))] = "ovarian_cancer_plasma_ND"

# methylation
methylation = read.csv(paste0(inputDir,"4_methylation_data/default_scores/10_largest/cg_methyl.all_samples.csv"))
rownames(methylation) = methylation$sample
methylation = methylation[rownames(c1),]
methylation$biological_class = c1$biological_class

# healthy only
methyl_healthy = methylation[which(methylation$biological_class=="healthy_plasma"),]

# OVC
methyl_cancer = methylation[which(methylation$biological_class=="early_ovarian_cancer_plasma"),]
m = rbind(methyl_healthy,methyl_cancer)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_d,"discovery_cohort.healthy_and_early_ovarian_cancer_plasma.methylation.no_l1hs15.csv"), row.names=FALSE)

# BRC
methyl_cancer = methylation[which(methylation$biological_class=="early_breast_cancer_plasma"),]
m = rbind(methyl_healthy,methyl_cancer)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_d,"discovery_cohort.healthy_and_early_breast_cancer_plasma.methylation.no_l1hs15.csv"), row.names=FALSE)

# GAC
methyl_cancer = methylation[which(methylation$biological_class=="early_gastric_cancer_plasma"),]
m = rbind(methyl_healthy,methyl_cancer)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_d,"discovery_cohort.healthy_and_early_gastric_cancer_plasma.methylation.no_l1hs15.csv"), row.names=FALSE)


##############
# Validation #
##############

# curated list C2
c2 = read.csv(paste0(inputDir,"0_data/c2_plasma_samples_info.csv"))
rownames(c2) = c2$sample
c2$biological_class[which(c2$biological_class=="ovarian_cancer_plasma" & c2$stage=="3")] = "early_ovarian_cancer_plasma"

# methylation
methylation = read.csv(paste0(inputDir,"4_methylation_data/default_scores/10_largest/cg_methyl.all_samples.csv"))
rownames(methylation) = methylation$sample
methylation = methylation[rownames(c2),]
methylation$biological_class = c2$biological_class

# healthy only
methyl_healthy = methylation[which(methylation$biological_class=="healthy_plasma"),]

# OVC
methyl_cancer = methylation[which(methylation$biological_class=="ovarian_cancer_plasma"),]
m = rbind(methyl_healthy,methyl_cancer)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_v,"validation_cohort.healthy_and_ovarian_cancer_plasma.methylation.no_l1hs15.csv"), row.names=FALSE)

# BRC
methyl_cancer = methylation[which(methylation$biological_class=="breast_cancer_plasma"),]
m = rbind(methyl_healthy,methyl_cancer)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_v,"validation_cohort.healthy_and_breast_cancer_plasma.methylation.no_l1hs15.csv"), row.names=FALSE)

# GAC
methyl_cancer = methylation[which(methylation$biological_class=="gastric_cancer_plasma"),]
m = rbind(methyl_healthy,methyl_cancer)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_v,"validation_cohort.healthy_and_gastric_cancer_plasma.methylation.no_l1hs15.csv"), row.names=FALSE)
