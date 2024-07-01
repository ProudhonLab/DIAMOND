#############
# Discovery #
#############

# Load data

# curated list C1
c1 = read.csv("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/groups/proudhon_lab/projects/l1pa_meth/0_data/c1_plasma_samples_info.csv")
rownames(c1) = c1$sample
to_del = c("D541S034","D541S036","D541S044","D541S058","D541S062","D541S065","D541S066","D541S067")
c1 = c1[-which(c1$sample %in% to_del),]

# methylation
methylation = read.csv("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/cg_methyl.all_samples.csv")
rownames(methylation) = methylation$sample
methylation = methylation[rownames(c1),]
methylation$biological_class = c1$biological_class

#ebrc = c1[which(c1$biological_class=="early_breast_cancer_plasma"),]
#meth_ebrc = methylation[rownames(ebrc),]
#write.csv(meth_ebrc,"Téléchargements/meth_ebrc.csv")

# haplotypes
haplotypes = read.csv("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/haplotypes.all_samples.csv")
rownames(haplotypes) = haplotypes$sample
haplotypes = haplotypes[rownames(c1),]
haplotypes$biological_class = c1$biological_class

# meth+haplo
meth_haplo = cbind(methylation,haplotypes[3:ncol(haplotypes)])
meth_haplo$biological_class = c1$biological_class

# healthy only
methyl_healthy = methylation[which(methylation$biological_class=="healthy_plasma"),]
haplo_healthy = haplotypes[which(haplotypes$biological_class=="healthy_plasma"),]
meth_haplo_healthy = meth_haplo[which(meth_haplo$biological_class=="healthy_plasma"),]

# Healthy vs multiCancer

methyl_cancer = methylation[which(methylation$biological_class!="healthy_plasma"),]
haplo_cancer = haplotypes[which(haplotypes$biological_class!="healthy_plasma"),]
meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class!="healthy_plasma"),]
m = rbind(methyl_healthy,methyl_cancer)
h = rbind(haplo_healthy,haplo_cancer)
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)
write.csv(m, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_multicancers.methylation.all_primers.csv", row.names=FALSE)
write.csv(h, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_multicancers.haplotypes.all_primers.csv", row.names=FALSE)
write.csv(mh, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_multicancers.methylation_and_haplotypes.all_primers.csv", row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_multicancers.methylation.no_l1hs15.csv", row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_multicancers.haplotypes.no_l1hs15.csv", row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_multicancers.methylation_and_haplotypes.no_l1hs15.csv", row.names=FALSE)

# Healthy vs Cancer

methyl_cancer = methylation[which(methylation$biological_class!="healthy_plasma"),]
haplo_cancer = haplotypes[which(haplotypes$biological_class!="healthy_plasma"),]
meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class!="healthy_plasma"),]
methyl_cancer$biological_class = haplo_cancer$biological_class = meth_haplo_cancer$biological_class = "cancer_plasma"
m = rbind(methyl_healthy,methyl_cancer)
h = rbind(haplo_healthy,haplo_cancer)
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)
write.csv(m, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_cancer.methylation.all_primers.csv", row.names=FALSE)
write.csv(h, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_cancer.haplotypes.all_primers.csv", row.names=FALSE)
write.csv(mh, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_cancer.methylation_and_haplotypes.all_primers.csv", row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_cancer.methylation.no_l1hs15.csv", row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_cancer.haplotypes.no_l1hs15.csv", row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_cancer.methylation_and_haplotypes.no_l1hs15.csv", row.names=FALSE)

# Healthy vs ebrc

cancer_type = "early_breast_cancer_plasma"
methyl_cancer = methylation[which(methylation$biological_class==cancer_type),]
haplo_cancer = haplotypes[which(haplotypes$biological_class==cancer_type),]
meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class==cancer_type),]

m = rbind(methyl_healthy,methyl_cancer)
h = rbind(haplo_healthy,haplo_cancer)
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)

write.csv(m, paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_",cancer_type,".methylation.all_primers.csv",sep=""), row.names=FALSE)
write.csv(h, paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_",cancer_type,".haplotypes.all_primers.csv",sep=""), row.names=FALSE)
write.csv(mh, paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_",cancer_type,".methylation_and_haplotypes.all_primers.csv",sep=""), row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_",cancer_type,".methylation.no_l1hs15.csv",sep=""), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_",cancer_type,".haplotypes.no_l1hs15.csv",sep=""), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs/discovery_cohort.healthy_and_",cancer_type,".methylation_and_haplotypes.no_l1hs15.csv",sep=""), row.names=FALSE)

# Healthy vs M0 vs M+

methyl_cancer = methylation[which(methylation$biological_class!="healthy_plasma"),]
haplo_cancer = haplotypes[which(haplotypes$biological_class!="healthy_plasma"),]
meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class!="healthy_plasma"),]
m = rbind(methyl_healthy,methyl_cancer)
m = m[rownames(c1),]
m$biological_class[which(m$biological_class=="breast_cancer_plasma")] = "mplus_plasma"
m$biological_class[which(m$biological_class=="colorectal_cancer_plasma")] = "mplus_plasma"
m$biological_class[which(m$biological_class=="early_breast_cancer_plasma")] = "mzero_plasma"
m$biological_class[which(m$biological_class=="gastric_cancer_plasma")] = "mzero_plasma"
m$biological_class[which(m$biological_class=="lung_cancer_plasma")] = "mplus_plasma"
m$biological_class[which(m$biological_class=="ovarian_cancer_plasma" & c1$stage=="III")] = "mzero_plasma"
m$biological_class[which(m$biological_class=="ovarian_cancer_plasma" & c1$stage=="IV")] = "mplus_plasma"
m = m[-which(m$biological_class=="ovarian_cancer_plasma"),]
m$biological_class[which(m$biological_class=="uveal_melanoma_cancer_plasma")] = "mplus_plasma"
h = rbind(haplo_healthy,haplo_cancer)
h = h[rownames(c1),]
h$biological_class[which(h$biological_class=="breast_cancer_plasma")] = "mplus_plasma"
h$biological_class[which(h$biological_class=="colorectal_cancer_plasma")] = "mplus_plasma"
h$biological_class[which(h$biological_class=="early_breast_cancer_plasma")] = "mzero_plasma"
h$biological_class[which(h$biological_class=="gastric_cancer_plasma")] = "mzero_plasma"
h$biological_class[which(h$biological_class=="lung_cancer_plasma")] = "mplus_plasma"
h$biological_class[which(h$biological_class=="ovarian_cancer_plasma" & c1$stage=="III")] = "mzero_plasma"
h$biological_class[which(h$biological_class=="ovarian_cancer_plasma" & c1$stage=="IV")] = "mplus_plasma"
h = h[-which(h$biological_class=="ovarian_cancer_plasma"),]
h$biological_class[which(h$biological_class=="uveal_melanoma_cancer_plasma")] = "mplus_plasma"
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)
mh = mh[rownames(c1),]
mh$biological_class[which(mh$biological_class=="breast_cancer_plasma")] = "mplus_plasma"
mh$biological_class[which(mh$biological_class=="colorectal_cancer_plasma")] = "mplus_plasma"
mh$biological_class[which(mh$biological_class=="early_breast_cancer_plasma")] = "mzero_plasma"
mh$biological_class[which(mh$biological_class=="gastric_cancer_plasma")] = "mzero_plasma"
mh$biological_class[which(mh$biological_class=="lung_cancer_plasma")] = "mplus_plasma"
mh$biological_class[which(mh$biological_class=="ovarian_cancer_plasma" & c1$stage=="III")] = "mzero_plasma"
mh$biological_class[which(mh$biological_class=="ovarian_cancer_plasma" & c1$stage=="IV")] = "mplus_plasma"
mh = mh[-which(mh$biological_class=="ovarian_cancer_plasma"),]
mh$biological_class[which(mh$biological_class=="uveal_melanoma_cancer_plasma")] = "mplus_plasma"

write.csv(m, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs_three/discovery_cohort.healthy_mzero_mplus.methylation.all_primers.csv", row.names=FALSE)
write.csv(h, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs_three/discovery_cohort.healthy_mzero_mplus.haplotypes.all_primers.csv", row.names=FALSE)
write.csv(mh, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs_three/discovery_cohort.healthy_mzero_mplus.methylation_and_haplotypes.all_primers.csv", row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs_three/discovery_cohort.healthy_mzero_mplus.methylation.no_l1hs15.csv", row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs_three/discovery_cohort.healthy_mzero_mplus.haplotypes.no_l1hs15.csv", row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/discovery/inputs_three/discovery_cohort.healthy_mzero_mplus.methylation_and_haplotypes.no_l1hs15.csv", row.names=FALSE)

##############
# Validation #
##############

# curated list C2
c2 = read.csv("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/groups/proudhon_lab/projects/l1pa_meth/0_data/c2_plasma_samples_info.csv")
rownames(c2) = c2$sample

# methylation
methylation = read.csv("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/cg_methyl.all_samples.csv")
rownames(methylation) = methylation$sample
methylation = methylation[rownames(c2),]

# haplotypes
haplotypes = read.csv("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/haplotypes.all_samples.csv")
rownames(haplotypes) = haplotypes$sample
haplotypes = haplotypes[rownames(c2),]

# meth+haplo
meth_haplo = cbind(methylation,haplotypes[3:ncol(haplotypes)])

# healthy only
methyl_healthy = methylation[which(methylation$biological_class=="healthy_plasma"),]
haplo_healthy = haplotypes[which(haplotypes$biological_class=="healthy_plasma"),]
meth_haplo_healthy = meth_haplo[which(meth_haplo$biological_class=="healthy_plasma"),]

# Multiclass

methyl_cancer = methylation[which(methylation$biological_class!="healthy_plasma"),]
haplo_cancer = haplotypes[which(haplotypes$biological_class!="healthy_plasma"),]
meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class!="healthy_plasma"),]
m = rbind(methyl_healthy,methyl_cancer)
h = rbind(haplo_healthy,haplo_cancer)
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)
write.csv(m, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs/validation_cohort.healthy_and_multicancers.methylation.all_primers.csv", row.names=FALSE)
write.csv(h, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs/validation_cohort.healthy_and_multicancers.haplotypes.all_primers.csv", row.names=FALSE)
write.csv(mh, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs/validation_cohort.healthy_and_multicancers.methylation_and_haplotypes.all_primers.csv", row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs/validation_cohort.healthy_and_multicancers.methylation.no_l1hs15.csv", row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs/validation_cohort.healthy_and_multicancers.haplotypes.no_l1hs15.csv", row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs/validation_cohort.healthy_and_multicancers.methylation_and_haplotypes.no_l1hs15.csv", row.names=FALSE)

# Healthy vs Cancer

methyl_cancer = methylation[which(methylation$biological_class!="healthy_plasma"),]
haplo_cancer = haplotypes[which(haplotypes$biological_class!="healthy_plasma"),]
meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class!="healthy_plasma"),]
methyl_cancer$biological_class = haplo_cancer$biological_class = meth_haplo_cancer$biological_class = "cancer_plasma"
m = rbind(methyl_healthy,methyl_cancer)
h = rbind(haplo_healthy,haplo_cancer)
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)
write.csv(m, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs_bin/validation_cohort.healthy_and_cancer.methylation.all_primers.csv", row.names=FALSE)
write.csv(h, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs_bin/validation_cohort.healthy_and_cancer.haplotypes.all_primers.csv", row.names=FALSE)
write.csv(mh, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs_bin/validation_cohort.healthy_and_cancer.methylation_and_haplotypes.all_primers.csv", row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs_bin/validation_cohort.healthy_and_cancer.methylation.no_l1hs15.csv", row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs_bin/validation_cohort.healthy_and_cancer.haplotypes.no_l1hs15.csv", row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs_bin/validation_cohort.healthy_and_cancer.methylation_and_haplotypes.no_l1hs15.csv", row.names=FALSE)

# Healthy vs ebrc

cancer_type = "early_breast_cancer_plasma"
methyl_cancer = methylation[which(methylation$biological_class==cancer_type),]
haplo_cancer = haplotypes[which(haplotypes$biological_class==cancer_type),]
meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class==cancer_type),]

m = rbind(methyl_healthy,methyl_cancer)
h = rbind(haplo_healthy,haplo_cancer)
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)

write.csv(m, paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs/validation_cohort.healthy_and_",cancer_type,".methylation.all_primers.csv",sep=""), row.names=FALSE)
write.csv(h, paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs/validation_cohort.healthy_and_",cancer_type,".haplotypes.all_primers.csv",sep=""), row.names=FALSE)
write.csv(mh, paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs/validation_cohort.healthy_and_",cancer_type,".methylation_and_haplotypes.all_primers.csv",sep=""), row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs/validation_cohort.healthy_and_",cancer_type,".methylation.no_l1hs15.csv",sep=""), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs/validation_cohort.healthy_and_",cancer_type,".haplotypes.no_l1hs15.csv",sep=""), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs/validation_cohort.healthy_and_",cancer_type,".methylation_and_haplotypes.no_l1hs15.csv",sep=""), row.names=FALSE)

# Healthy vs M0 vs M+

methyl_cancer = methylation[which(methylation$biological_class!="healthy_plasma"),]
haplo_cancer = haplotypes[which(haplotypes$biological_class!="healthy_plasma"),]
meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class!="healthy_plasma"),]
m = rbind(methyl_healthy,methyl_cancer)
m = m[rownames(c2),]
m$biological_class[which(m$biological_class=="breast_cancer_plasma")] = "mplus_plasma"
m$biological_class[which(m$biological_class=="colorectal_cancer_plasma")] = "mplus_plasma"
m$biological_class[which(m$biological_class=="early_breast_cancer_plasma")] = "mzero_plasma"
m$biological_class[which(m$biological_class=="gastric_cancer_plasma")] = "mplus_plasma"
m$biological_class[which(m$biological_class=="lung_cancer_plasma")] = "mplus_plasma"
m$biological_class[which(m$biological_class=="ovarian_cancer_plasma" & c2$stage=="3")] = "mzero_plasma"
m$biological_class[which(m$biological_class=="ovarian_cancer_plasma" & c2$stage=="4")] = "mplus_plasma"
m$biological_class[which(m$biological_class=="uveal_melanoma_cancer_plasma")] = "mplus_plasma"
h = rbind(haplo_healthy,haplo_cancer)
h = h[rownames(c2),]
h$biological_class[which(h$biological_class=="breast_cancer_plasma")] = "mplus_plasma"
h$biological_class[which(h$biological_class=="colorectal_cancer_plasma")] = "mplus_plasma"
h$biological_class[which(h$biological_class=="early_breast_cancer_plasma")] = "mzero_plasma"
h$biological_class[which(h$biological_class=="gastric_cancer_plasma")] = "mplus_plasma"
h$biological_class[which(h$biological_class=="lung_cancer_plasma")] = "mplus_plasma"
h$biological_class[which(h$biological_class=="ovarian_cancer_plasma" & c2$stage=="3")] = "mzero_plasma"
h$biological_class[which(h$biological_class=="ovarian_cancer_plasma" & c2$stage=="4")] = "mplus_plasma"
h$biological_class[which(h$biological_class=="uveal_melanoma_cancer_plasma")] = "mplus_plasma"
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)
mh = mh[rownames(c2),]
mh$biological_class[which(mh$biological_class=="breast_cancer_plasma")] = "mplus_plasma"
mh$biological_class[which(mh$biological_class=="colorectal_cancer_plasma")] = "mplus_plasma"
mh$biological_class[which(mh$biological_class=="early_breast_cancer_plasma")] = "mzero_plasma"
mh$biological_class[which(mh$biological_class=="gastric_cancer_plasma")] = "mplus_plasma"
mh$biological_class[which(mh$biological_class=="lung_cancer_plasma")] = "mplus_plasma"
mh$biological_class[which(mh$biological_class=="ovarian_cancer_plasma" & c2$stage=="3")] = "mzero_plasma"
mh$biological_class[which(mh$biological_class=="ovarian_cancer_plasma" & c2$stage=="4")] = "mplus_plasma"
mh$biological_class[which(mh$biological_class=="uveal_melanoma_cancer_plasma")] = "mplus_plasma"

write.csv(m, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs_three/validation_cohort.healthy_mzero_mplus.methylation.all_primers.csv", row.names=FALSE)
write.csv(h, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs_three/validation_cohort.healthy_mzero_mplus.haplotypes.all_primers.csv", row.names=FALSE)
write.csv(mh, "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs_three/validation_cohort.healthy_mzero_mplus.methylation_and_haplotypes.all_primers.csv", row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs_three/validation_cohort.healthy_mzero_mplus.methylation.no_l1hs15.csv", row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs_three/validation_cohort.healthy_mzero_mplus.haplotypes.no_l1hs15.csv", row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], "/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/rerun_classif_noBadEbrc2/validation/inputs_three/validation_cohort.healthy_mzero_mplus.methylation_and_haplotypes.no_l1hs15.csv", row.names=FALSE)
