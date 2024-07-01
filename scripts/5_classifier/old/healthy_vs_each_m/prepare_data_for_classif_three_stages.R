##################
# May 25th, 2023 #
##################
# We want to do classifications based on the cancer stages
# with stage 1-2 = early
# stage 3 = locally advanced
# stage 4 = meta
# biological_class becomes the simplified stages

# !!! change this !!!
user = "kdasilva"
# !!!

workDir_d = paste0("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=",user,"/scratch/",user,"/three_stages_classification/discovery/inputs/")
workDir_v = paste0("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=",user,"/scratch/",user,"/three_stages_classification/validation/inputs/")

#############
# DISCOVERY #
#############

c1_infos = read.csv(paste0("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=",user,"/groups/proudhon_lab/projects/l1pa_meth/0_data/c1_plasma_samples_info.csv"))
rownames(c1_infos) = c1_infos$sample
c1_methyl = read.csv(paste0("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=",user,"/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/cg_methyl.c1.csv"))
c1_methyl = c1_methyl[,-1]
rownames(c1_methyl) = c1_methyl$sample
c1_methyl = c1_methyl[rownames(c1_infos),-grep("L1HS_15",colnames(c1_methyl))]
c1_methyl$biological_class[which(c1_methyl$biological_class!="healthy_plasma")] = c1_infos$stage[which(c1_infos$biological_class!="healthy_plasma")]
c1_methyl$biological_class[which(c1_methyl$biological_class=="healthy_plasma")] = "healthy"
c1_methyl$biological_class[which(c1_methyl$biological_class=="1" | c1_methyl$biological_class=="2")] = "early"
c1_methyl$biological_class[which(c1_methyl$biological_class=="3")] = "advanced"
c1_methyl$biological_class[which(c1_methyl$biological_class=="4")] = "meta"
c1_methyl = rbind(c1_methyl[which(c1_methyl$biological_class=="healthy"),],c1_methyl[which(c1_methyl$biological_class=="early"),],c1_methyl[which(c1_methyl$biological_class=="advanced"),],c1_methyl[which(c1_methyl$biological_class=="meta"),]) # healthy class first
c1_methyl$biological_class = factor(c1_methyl$biological_class,levels=c("healthy","early","advanced","meta"))
healthy_vs_early = rbind(c1_methyl[which(c1_methyl$biological_class=="healthy"),],c1_methyl[which(c1_methyl$biological_class=="early"),])
healthy_vs_advanced = rbind(c1_methyl[which(c1_methyl$biological_class=="healthy"),],c1_methyl[which(c1_methyl$biological_class=="advanced"),])
healthy_vs_meta = rbind(c1_methyl[which(c1_methyl$biological_class=="healthy"),],c1_methyl[which(c1_methyl$biological_class=="meta"),])

write.csv(healthy_vs_early, paste0(workDir_d,"discovery_cohort.healthy_vs_early.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(healthy_vs_advanced, paste0(workDir_d,"discovery_cohort.healthy_vs_advanced.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(healthy_vs_meta, paste0(workDir_d,"discovery_cohort.healthy_vs_meta.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(c1_methyl, paste0(workDir_d,"discovery_cohort.multiclass_stages.methylation.no_l1hs15.csv"), row.names=FALSE)

##############
# VALIDATION #
##############

c2_infos = read.csv(paste0("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=",user,"/groups/proudhon_lab/projects/l1pa_meth/0_data/c2_plasma_samples_info.csv"))
rownames(c2_infos) = c2_infos$sample
c2_methyl = read.csv(paste0("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=",user,"/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/cg_methyl.c2.csv"))
c2_methyl = c2_methyl[,-1]
rownames(c2_methyl) = c2_methyl$sample
c2_methyl = c2_methyl[rownames(c2_infos),-grep("L1HS_15",colnames(c2_methyl))]
c2_methyl$biological_class[which(c2_methyl$biological_class!="healthy_plasma")] = c2_infos$stage[which(c2_infos$biological_class!="healthy_plasma")]
c2_methyl$biological_class[which(c2_methyl$biological_class=="healthy_plasma")] = "healthy"
c2_methyl$biological_class[which(c2_methyl$biological_class=="1" | c2_methyl$biological_class=="2")] = "early"
c2_methyl$biological_class[which(c2_methyl$biological_class=="3")] = "advanced"
c2_methyl$biological_class[which(c2_methyl$biological_class=="4")] = "meta"
c2_methyl = rbind(c2_methyl[which(c2_methyl$biological_class=="healthy"),],c2_methyl[which(c2_methyl$biological_class=="early"),],c2_methyl[which(c2_methyl$biological_class=="advanced"),],c2_methyl[which(c2_methyl$biological_class=="meta"),]) # healthy class first
c2_methyl$biological_class = factor(c2_methyl$biological_class,levels=c("healthy","early","advanced","meta"))
healthy_vs_early = rbind(c2_methyl[which(c2_methyl$biological_class=="healthy"),],c2_methyl[which(c2_methyl$biological_class=="early"),])
healthy_vs_advanced = rbind(c2_methyl[which(c2_methyl$biological_class=="healthy"),],c2_methyl[which(c2_methyl$biological_class=="advanced"),])
healthy_vs_meta = rbind(c2_methyl[which(c2_methyl$biological_class=="healthy"),],c2_methyl[which(c2_methyl$biological_class=="meta"),])

write.csv(healthy_vs_early, paste0(workDir_v,"validation_cohort.healthy_vs_early.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(healthy_vs_advanced, paste0(workDir_v,"validation_cohort.healthy_vs_advanced.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(healthy_vs_meta, paste0(workDir_v,"validation_cohort.healthy_vs_meta.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(c2_methyl, paste0(workDir_v,"validation_cohort.multiclass_stages.methylation.no_l1hs15.csv"), row.names=FALSE)
