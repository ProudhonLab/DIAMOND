# !!! change this !!!
user = "kdasilva"
inputDir = paste0("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=",user,"/groups/proudhon_lab/projects/l1pa_meth/")
workDir_d = paste0("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=",user,"/scratch/kdasilva/20230612_classic_classifications/discovery/inputs/")
workDir_v = paste0("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=",user,"/scratch/kdasilva/20230612_classic_classifications/validation/inputs/")

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

# haplotypes
haplotypes = read.csv(paste0(inputDir,"4_methylation_data/default_scores/10_largest/haplotypes.all_samples.csv"))
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
#write.csv(m, paste0(workDir_d,"discovery_cohort.healthy_and_multicancers.methylation.all_primers.csv"), row.names=FALSE)
#write.csv(h, paste0(workDir_d,"discovery_cohort.healthy_and_multicancers.haplotypes.all_primers.csv"), row.names=FALSE)
#write.csv(mh, paste0(workDir_d,"discovery_cohort.healthy_and_multicancers.methylation_and_haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_d,"discovery_cohort.healthy_and_multicancers.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste0(workDir_d,"discovery_cohort.healthy_and_multicancers.haplotypes.no_l1hs15.csv"), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste0(workDir_d,"discovery_cohort.healthy_and_multicancers.methylation_and_haplotypes.no_l1hs15.csv"), row.names=FALSE)

# Healthy vs multiCancer (without uveal)

methyl_cancer = methylation[which(methylation$biological_class!="healthy_plasma" & methylation$biological_class!="uveal_melanoma_cancer_plasma"),]
haplo_cancer = haplotypes[which(haplotypes$biological_class!="healthy_plasma" & methylation$biological_class!="uveal_melanoma_cancer_plasma"),]
meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class!="healthy_plasma" & methylation$biological_class!="uveal_melanoma_cancer_plasma"),]
m = rbind(methyl_healthy,methyl_cancer)
h = rbind(haplo_healthy,haplo_cancer)
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)
write.csv(m, paste0(workDir_d,"discovery_cohort.healthy_and_multicancers_without_uveal.methylation.all_primers.csv"), row.names=FALSE)
write.csv(h, paste0(workDir_d,"discovery_cohort.healthy_and_multicancers_without_uveal.haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(mh, paste0(workDir_d,"discovery_cohort.healthy_and_multicancers_without_uveal.methylation_and_haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_d,"discovery_cohort.healthy_and_multicancers_without_uveal.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste0(workDir_d,"discovery_cohort.healthy_and_multicancers_without_uveal.haplotypes.no_l1hs15.csv"), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste0(workDir_d,"discovery_cohort.healthy_and_multicancers_without_uveal.methylation_and_haplotypes.no_l1hs15.csv"), row.names=FALSE)

# Healthy vs Cancer

methyl_cancer = methylation[which(methylation$biological_class!="healthy_plasma"),]
haplo_cancer = haplotypes[which(haplotypes$biological_class!="healthy_plasma"),]
meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class!="healthy_plasma"),]
methyl_cancer$biological_class = haplo_cancer$biological_class = meth_haplo_cancer$biological_class = "cancer_plasma"
m = rbind(methyl_healthy,methyl_cancer)
h = rbind(haplo_healthy,haplo_cancer)
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)
#write.csv(m, paste0(workDir_d,"discovery_cohort.healthy_and_cancer.methylation.all_primers.csv"), row.names=FALSE)
#write.csv(h, paste0(workDir_d,"discovery_cohort.healthy_and_cancer.haplotypes.all_primers.csv"), row.names=FALSE)
#write.csv(mh, paste0(workDir_d,"discovery_cohort.healthy_and_cancer.methylation_and_haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_d,"discovery_cohort.healthy_and_cancer.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste0(workDir_d,"discovery_cohort.healthy_and_cancer.haplotypes.no_l1hs15.csv"), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste0(workDir_d,"discovery_cohort.healthy_and_cancer.methylation_and_haplotypes.no_l1hs15.csv"), row.names=FALSE)

# Healthy vs each Cancer 

bioclass = unique(c1$biological_class)[-c(6,10)] # remove healthy and ovarian cancer ND

for (i in 1:length(bioclass)) {
  cancer_type = bioclass[i]
  methyl_cancer = methylation[which(methylation$biological_class==cancer_type),]
  haplo_cancer = haplotypes[which(haplotypes$biological_class==cancer_type),]
  meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class==cancer_type),]
  
  m = rbind(methyl_healthy,methyl_cancer)
  h = rbind(haplo_healthy,haplo_cancer)
  mh = rbind(meth_haplo_healthy,meth_haplo_cancer)
  
  #write.csv(m, paste0(workDir_d,"discovery_cohort.healthy_and_",cancer_type,".methylation.all_primers.csv"), row.names=FALSE)
  #write.csv(h, paste0(workDir_d,"discovery_cohort.healthy_and_",cancer_type,".haplotypes.all_primers.csv"), row.names=FALSE)
  #write.csv(mh, paste0(workDir_d,"discovery_cohort.healthy_and_",cancer_type,".methylation_and_haplotypes.all_primers.csv"), row.names=FALSE)
  write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_d,"discovery_cohort.healthy_and_",cancer_type,".methylation.no_l1hs15.csv"), row.names=FALSE)
  write.csv(h[,-grep("L1HS_15",colnames(h))], paste0(workDir_d,"discovery_cohort.healthy_and_",cancer_type,".haplotypes.no_l1hs15.csv"), row.names=FALSE)
  write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste0(workDir_d,"discovery_cohort.healthy_and_",cancer_type,".methylation_and_haplotypes.no_l1hs15.csv"), row.names=FALSE)
  
}

# All OVC (M0 and M+ and ND)

methyl_cancer = methylation[grep("ovarian_cancer_plasma",methylation$biological_class),]
haplo_cancer = haplotypes[grep("ovarian_cancer_plasma",haplotypes$biological_class),]
meth_haplo_cancer = meth_haplo[grep("ovarian_cancer_plasma",meth_haplo$biological_class),]
methyl_cancer$biological_class = haplo_cancer$biological_class = meth_haplo_cancer$biological_class = "ovarian_cancer_plasma"

m = rbind(methyl_healthy,methyl_cancer)
h = rbind(haplo_healthy,haplo_cancer)
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)

#write.csv(m, paste0(workDir_d,"discovery_cohort.healthy_and_ovarian_cancer_plasma_stage3_only.methylation.all_primers.csv"), row.names=FALSE)
#write.csv(h, paste0(workDir_d,"discovery_cohort.healthy_and_ovarian_cancer_plasma_stage3_only.haplotypes.all_primers.csv"), row.names=FALSE)
#write.csv(mh, paste0(workDir_d,"discovery_cohort.healthy_and_ovarian_cancer_plasma_stage3_only.methylation_and_haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_d,"discovery_cohort.healthy_and_all_ovarian_cancer_plasma.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste0(workDir_d,"discovery_cohort.healthy_and_all_ovarian_cancer_plasma.haplotypes.no_l1hs15.csv"), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste0(workDir_d,"discovery_cohort.healthy_and_all_ovarian_cancer_plasma.methylation_and_haplotypes.no_l1hs15.csv"), row.names=FALSE)

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
mh$biological_class[which(mh$biological_class=="uveal_melanoma_cancer_plasma")] = "mplus_plasma"

write.csv(m, paste0(workDir_d,"discovery_cohort.healthy_mzero_mplus.methylation.all_primers.csv"), row.names=FALSE)
write.csv(h, paste0(workDir_d,"discovery_cohort.healthy_mzero_mplus.haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(mh, paste0(workDir_d,"discovery_cohort.healthy_mzero_mplus.methylation_and_haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_d,"discovery_cohort.healthy_mzero_mplus.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste0(workDir_d,"discovery_cohort.healthy_mzero_mplus.haplotypes.no_l1hs15.csv"), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste0(workDir_d,"discovery_cohort.healthy_mzero_mplus.methylation_and_haplotypes.no_l1hs15.csv"), row.names=FALSE)

# Healthy, eBRC, BRC, CRC

m = rbind(methyl_healthy,methylation[which(methylation$biological_class=="early_breast_cancer_plasma" | methylation$biological_class=="breast_cancer_plasma" | methylation$biological_class=="colorectal_cancer_plasma"),])
h = rbind(haplo_healthy,haplotypes[which(haplotypes$biological_class=="early_breast_cancer_plasma" | haplotypes$biological_class=="breast_cancer_plasma" | haplotypes$biological_class=="colorectal_cancer_plasma"),])
mh = rbind(meth_haplo_healthy,meth_haplo[which(meth_haplo$biological_class=="early_breast_cancer_plasma" | meth_haplo$biological_class=="breast_cancer_plasma" | meth_haplo$biological_class=="colorectal_cancer_plasma"),])

write.csv(m, paste0(workDir_d,"discovery_cohort.healthy_ebrc_brc_crc.methylation.all_primers.csv"), row.names=FALSE)
write.csv(h, paste0(workDir_d,"discovery_cohort.healthy_ebrc_brc_crc.haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(mh, paste0(workDir_d,"discovery_cohort.healthy_ebrc_brc_crc.methylation_and_haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_d,"discovery_cohort.healthy_ebrc_brc_crc.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste0(workDir_d,"discovery_cohort.healthy_ebrc_brc_crc.haplotypes.no_l1hs15.csv"), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste0(workDir_d,"discovery_cohort.healthy_ebrc_brc_crc.methylation_and_haplotypes.no_l1hs15.csv"), row.names=FALSE)

# eBRC only

m = rbind(methyl_healthy,methylation[which(methylation$biological_class=="early_breast_cancer_plasma"),])
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_d,"discovery_cohort.healthy_ebrc.methylation.no_l1hs15.csv"), row.names=FALSE)

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

# haplotypes
haplotypes = read.csv(paste0(inputDir,"4_methylation_data/default_scores/10_largest/haplotypes.all_samples.csv"))
rownames(haplotypes) = haplotypes$sample
haplotypes = haplotypes[rownames(c2),]
haplotypes$biological_class = c2$biological_class

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
#write.csv(m, paste0(workDir_v,"validation_cohort.healthy_and_multicancers.methylation.all_primers.csv"), row.names=FALSE)
#write.csv(h, paste0(workDir_v,"validation_cohort.healthy_and_multicancers.haplotypes.all_primers.csv"), row.names=FALSE)
#write.csv(mh, paste0(workDir_v,"validation_cohort.healthy_and_multicancers.methylation_and_haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_v,"validation_cohort.healthy_and_multicancers.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste0(workDir_v,"validation_cohort.healthy_and_multicancers.haplotypes.no_l1hs15.csv"), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste0(workDir_v,"validation_cohort.healthy_and_multicancers.methylation_and_haplotypes.no_l1hs15.csv"), row.names=FALSE)

# Multiclass (remove gastric, lung, ovarian, uveal)

methyl_cancer = methylation[which(methylation$biological_class!="healthy_plasma" & methylation$biological_class!="gastric_cancer_plasma" & methylation$biological_class!="lung_cancer_plasma" & methylation$biological_class!="ovarian_cancer_plasma" & methylation$biological_class!="uveal_melanoma_plasma"),]
haplo_cancer = haplotypes[which(methylation$biological_class!="healthy_plasma" & methylation$biological_class!="gastric_cancer_plasma" & methylation$biological_class!="lung_cancer_plasma" & methylation$biological_class!="ovarian_cancer_plasma" & methylation$biological_class!="uveal_melanoma_plasma"),]
meth_haplo_cancer = meth_haplo[which(methylation$biological_class!="healthy_plasma" & methylation$biological_class!="gastric_cancer_plasma" & methylation$biological_class!="lung_cancer_plasma" & methylation$biological_class!="ovarian_cancer_plasma" & methylation$biological_class!="uveal_melanoma_plasma"),]
m = rbind(methyl_healthy,methyl_cancer)
h = rbind(haplo_healthy,haplo_cancer)
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)
write.csv(m, paste0(workDir_v,"validation_cohort.healthy_and_multicancers_some_removed.methylation.all_primers.csv"), row.names=FALSE)
write.csv(h, paste0(workDir_v,"validation_cohort.healthy_and_multicancers_some_removed.haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(mh, paste0(workDir_v,"validation_cohort.healthy_and_multicancers_some_removed.methylation_and_haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_v,"validation_cohort.healthy_and_multicancers_some_removed.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste0(workDir_v,"validation_cohort.healthy_and_multicancers_some_removed.haplotypes.no_l1hs15.csv"), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste0(workDir_v,"validation_cohort.healthy_and_multicancers_some_removed.methylation_and_haplotypes.no_l1hs15.csv"), row.names=FALSE)

# Healthy vs Cancer

methyl_cancer = methylation[which(methylation$biological_class!="healthy_plasma"),]
haplo_cancer = haplotypes[which(haplotypes$biological_class!="healthy_plasma"),]
meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class!="healthy_plasma"),]
methyl_cancer$biological_class = haplo_cancer$biological_class = meth_haplo_cancer$biological_class = "cancer_plasma"
m = rbind(methyl_healthy,methyl_cancer)
h = rbind(haplo_healthy,haplo_cancer)
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)
#write.csv(m, paste0(workDir_v,"validation_cohort.healthy_and_cancer.methylation.all_primers.csv"), row.names=FALSE)
#write.csv(h, paste0(workDir_v,"validation_cohort.healthy_and_cancer.haplotypes.all_primers.csv"), row.names=FALSE)
#write.csv(mh, paste0(workDir_v,"validation_cohort.healthy_and_cancer.methylation_and_haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_v,"validation_cohort.healthy_and_cancer.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste0(workDir_v,"validation_cohort.healthy_and_cancer.haplotypes.no_l1hs15.csv"), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste0(workDir_v,"validation_cohort.healthy_and_cancer.methylation_and_haplotypes.no_l1hs15.csv"), row.names=FALSE)

# Healthy vs each Cancer 

bioclass = unique(c2$biological_class)[-6] # remove healthy

for (i in 1:length(bioclass)) {
  cancer_type = bioclass[i]
  methyl_cancer = methylation[which(methylation$biological_class==cancer_type),]
  haplo_cancer = haplotypes[which(haplotypes$biological_class==cancer_type),]
  meth_haplo_cancer = meth_haplo[which(meth_haplo$biological_class==cancer_type),]
  
  m = rbind(methyl_healthy,methyl_cancer)
  h = rbind(haplo_healthy,haplo_cancer)
  mh = rbind(meth_haplo_healthy,meth_haplo_cancer)
  
  #write.csv(m, paste(paste0(workDir_v,"validation_cohort.healthy_and_",cancer_type,".methylation.all_primers.csv"),sep=""), row.names=FALSE)
  #write.csv(h, paste(paste0(workDir_v,"validation_cohort.healthy_and_",cancer_type,".haplotypes.all_primers.csv"),sep=""), row.names=FALSE)
  #write.csv(mh, paste(paste0(workDir_v,"validation_cohort.healthy_and_",cancer_type,".methylation_and_haplotypes.all_primers.csv"),sep=""), row.names=FALSE)
  write.csv(m[,-grep("L1HS_15",colnames(m))], paste(paste0(workDir_v,"validation_cohort.healthy_and_",cancer_type,".methylation.no_l1hs15.csv"),sep=""), row.names=FALSE)
  write.csv(h[,-grep("L1HS_15",colnames(h))], paste(paste0(workDir_v,"validation_cohort.healthy_and_",cancer_type,".haplotypes.no_l1hs15.csv"),sep=""), row.names=FALSE)
  write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste(paste0(workDir_v,"validation_cohort.healthy_and_",cancer_type,".methylation_and_haplotypes.no_l1hs15.csv"),sep=""), row.names=FALSE)
  
}

# All OVC (M0 and M+ and ND)

methyl_cancer = methylation[grep("ovarian_cancer_plasma",methylation$biological_class),]
haplo_cancer = haplotypes[grep("ovarian_cancer_plasma",haplotypes$biological_class),]
meth_haplo_cancer = meth_haplo[grep("ovarian_cancer_plasma",meth_haplo$biological_class),]
methyl_cancer$biological_class = haplo_cancer$biological_class = meth_haplo_cancer$biological_class = "ovarian_cancer_plasma"

m = rbind(methyl_healthy,methyl_cancer)
h = rbind(haplo_healthy,haplo_cancer)
mh = rbind(meth_haplo_healthy,meth_haplo_cancer)

write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_v,"validation_cohort.healthy_and_all_ovarian_cancer_plasma.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste0(workDir_v,"validation_cohort.healthy_and_all_ovarian_cancer_plasma.haplotypes.no_l1hs15.csv"), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste0(workDir_v,"validation_cohort.healthy_and_all_ovarian_cancer_plasma.methylation_and_haplotypes.no_l1hs15.csv"), row.names=FALSE)


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

write.csv(m, paste0(workDir_v,"validation_cohort.healthy_mzero_mplus.methylation.all_primers.csv"), row.names=FALSE)
write.csv(h, paste0(workDir_v,"validation_cohort.healthy_mzero_mplus.haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(mh, paste0(workDir_v,"validation_cohort.healthy_mzero_mplus.methylation_and_haplotypes.all_primers.csv"), row.names=FALSE)
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_v,"validation_cohort.healthy_mzero_mplus.methylation.no_l1hs15.csv"), row.names=FALSE)
write.csv(h[,-grep("L1HS_15",colnames(h))], paste0(workDir_v,"validation_cohort.healthy_mzero_mplus.haplotypes.no_l1hs15.csv"), row.names=FALSE)
write.csv(mh[,-grep("L1HS_15",colnames(mh))], paste0(workDir_v,"validation_cohort.healthy_mzero_mplus.methylation_and_haplotypes.no_l1hs15.csv"), row.names=FALSE)

# only eBRC stade3

m = rbind(methyl_healthy,methylation[which(methylation$biological_class=="early_breast_cancer_plasma" & (c2$stage=="IIIA" | c2$stage=="IIIB" | c2$stage=="IIIC")),])
write.csv(m[,-grep("L1HS_15",colnames(m))], paste0(workDir_v,"validation_cohort.healthy_and_ebrc_stageIII.methylation.no_l1hs15.csv"), row.names=FALSE)
