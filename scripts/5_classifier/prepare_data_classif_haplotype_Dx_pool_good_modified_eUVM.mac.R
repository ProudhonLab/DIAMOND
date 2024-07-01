##################
# October 19th, 2023 #
##################
# We want to do classifications with haplotype only + Dx + pools 
# This version is made for MAC version (run in GenOuest)

# !!! change this !!!
user = "mgorse"
# !!!

###############
# haplotypes #
##############

workDir_d = paste0("/scratch/",user,"/20231025_classification_haplotype_Dx_pool_good/discovery/inputs/")
workDir_v = paste0("/scratch/",user,"/20231025_classification_haplotype_Dx_pool_good/validation/inputs/")

# Discovery

c1_infos = read.csv("/groups/proudhon_lab/projects/l1pa_meth/0_data/c1_plasma_samples_info.csv")

rownames(c1_infos) = c1_infos$sample
c1_infos$biological_class[which(c1_infos$biological_class=="gastric_cancer_plasma")] = "early_gastric_cancer_plasma"
c1_infos$biological_class[which(c1_infos$biological_class=="ovarian_cancer_plasma" & c1_infos$stage=="3")] = "early_ovarian_cancer_plasma"
c1_infos$biological_class[which(c1_infos$biological_class=="ovarian_cancer_plasma" & is.na(c1_infos$stage))] = "ovarian_cancer_plasma_ND"


# Ajout des haplotypes C1 
haplo_c1 = read.csv("/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/haplotypes.c1.csv")
#haplo_c1 = read.csv("haplotypes.c1.csv")

rownames(haplo_c1) = haplo_c1$sample
haplo_c1 = haplo_c1[,-grep("L1HS_15",colnames(haplo_c1))]
haplo_c1 = haplo_c1[which(haplo_c1$sample %in% c1_infos$sample),]
haplo_c1 = haplo_c1[,-1]
haplo_c1$biological_class <- c1_infos$biological_class

# Ajout des eUVM de MHOVC25 en discovery
infos_all= read.csv("/groups/proudhon_lab/projects/l1pa_meth/0_data/samples_infos.csv")
#infos_all= read.csv("samples_infos.csv")
rownames(infos_all) = infos_all$sample
infos_eUVM = infos_all[which(infos_all$mhovc == "mhovc25" & infos_all$type == "plasma" & infos_all$class != "healthy" & infos_all$paper_cohort == "discovery"),]  


haplo_mhovc25 = read.csv("/groups/proudhon_lab/projects/l1pa_meth_v2/3_methylation_data_vs_all/haplotypes.all_samples.csv")
#haplo_mhovc25 = read.csv("haplotypes.all_samples.csv")

#c1d_haplo = haplo_mhovc25[which(haplo_mhovc25$sample %in% infos_eUVM$sample),]
#c1d_haplo = c1d_haplo[,-grep("L1HS_15",colnames(c1d_haplo))]
#rownames(c1d_haplo) = c1d_haplo$sample


# Merge de toutes les data de discovery
#c1_methyl =  rbind(haplo_c1, c1d_haplo)
c1_methyl =  haplo_c1

# Séparation des data en HD - cancer
healthy_methyl = c1_methyl[which(c1_methyl$biological_class=="healthy_plasma"),]

methyl_cancer = c1_methyl[which(c1_methyl$biological_class!="healthy_plasma" & c1_methyl$biological_class != "early_uveal_melanoma_cancer_plasma" & c1_methyl$biological_class != "early_gastric_cancer_plasma"),]
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

c2_infos = read.csv("/groups/proudhon_lab/projects/l1pa_meth/0_data/c2_plasma_samples_info.csv")
#c2_infos = read.csv("c2_plasma_samples_info.csv")

rownames(c2_infos) = c2_infos$sample
c2_infos$biological_class[which(c2_infos$biological_class=="ovarian_cancer_plasma" & c2_infos$stage=="3")] = "early_ovarian_cancer_plasma"


haplo_c2 = read.csv("/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/haplotypes.c2.csv")
#haplo_c2 = read.csv("haplotypes.c2.csv")
haplo_c2 = haplo_c2[,-grep("L1HS_15",colnames(haplo_c2))]
haplo_c2$biological_class <- c2_infos$biological_class


# Ajout des pools non douteux
methyl = read.csv(paste0("/groups/proudhon_lab/projects/extraction_methods_comparison/3_methylation_data/haplotypes.all_samples.csv"))
#methyl = read.csv(paste0("haplotypes.all_samples_pool.csv"))
infos = read.table(paste0("/groups/proudhon_lab/projects/extraction_methods_comparison/0_data/sequencing_data/D1357.sampleDescription.txt"),sep="|")
#infos = read.table(paste0("D1357.sampleDescription.txt"),sep="|")
methyl = methyl[,-grep("L1HS_15",colnames(methyl))]
rownames(methyl) = methyl$sample
methyl$name = infos$V2[match(methyl$sample,infos$V1)]
methyl$individual = gsub("(.*)_.*","\\1",methyl$name)
methyl$method = gsub(".*_(.*)","\\1",methyl$name)
methyl = methyl[-which(methyl$method=="P" | methyl$method=="100" | methyl$method=="50" | methyl$method=="0"),]
methyl$biological_class = "healthy_plasma"
c2c_methyl = methyl[which(methyl$method=="Mn"),1:374]

c2c_methyl = c2c_methyl[which(!(rownames(c2c_methyl) %in% c("D1357S148", "D1357S126", "D1357S128", "D1357S135", "D1357S127", "D1357S142", "D1357S136", "D1357S150", "D1357S141", "D1357S143", "D1357S149"))),]

# Ajout des cancers eUVM et eGAC

c2_infos_mhovc25 = infos_all[which(infos_all$mhovc == "mhovc25" & infos_all$type == "plasma" & infos_all$paper_cohort == "validation" & infos_all$patient_cohort != "EFS"),]
print(unique(c2_infos_mhovc25$biological_class))

c2_infos_mhovc25=c2_infos_mhovc25[which(!(rownames(c2_infos_mhovc25) %in% c("D1525S028","D1525S020","D1525S021","D1525S038","D1525S039"))),]

c2d_haplo = haplo_mhovc25[which(haplo_mhovc25$sample %in% c2_infos_mhovc25$sample),]
c2d_haplo = c2d_haplo[,-grep("L1HS_15",colnames(c2d_haplo))]


c2_methyl = rbind(haplo_c2[,-1], c2d_haplo, c2c_methyl)
c2_methyl$biological_class = factor(c2_methyl$biological_class)



healthy_methyl = c2_methyl[which(c2_methyl$biological_class=="healthy_plasma"),]

methyl_cancer = c2_methyl[which(c2_methyl$biological_class!="healthy_plasma" & c2_methyl$biological_class != "early_uveal_melanoma_cancer_plasma" & c2_methyl$biological_class != "early_gastric_cancer_plasma"),]
print(unique(methyl_cancer$biological_class))
methyl_cancer$biological_class = "cancer_plasma"
m = rbind(healthy_methyl,methyl_cancer)
write.csv(m,  paste0(workDir_v,"validation_cohort.healthy_vs_cancer.methylation.no_l1hs15.csv"), row.names=FALSE)

print(levels(c2_methyl$biological_class))

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




