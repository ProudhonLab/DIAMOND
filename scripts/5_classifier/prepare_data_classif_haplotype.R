##################
# October 9th, 2023 #
##################
# We want to do classifications with methylation and haplotypes
# This version is made for MAC version (run in GenOuest)

# !!! change this !!!
user = "mgorse"
# !!!

##############################################
#  Methylation + haplotypes #
##############################################

workDir_d = paste0("/scratch/",user,"/20231009_classification_haplotype/discovery/inputs/")
workDir_v = paste0("/scratch/",user,"/20231009_classification_haplotype/validation/inputs/")

# Discovery

c1_infos = read.csv("/groups/proudhon_lab/projects/l1pa_meth/0_data/c1_plasma_samples_info.csv")
#c1_infos = read.csv("c1_plasma_samples_info.csv")

rownames(c1_infos) = c1_infos$sample
c1_infos$biological_class[which(c1_infos$biological_class=="gastric_cancer_plasma")] = "early_gastric_cancer_plasma"
c1_infos$biological_class[which(c1_infos$biological_class=="ovarian_cancer_plasma" & c1_infos$stage=="3")] = "early_ovarian_cancer_plasma"
c1_infos$biological_class[which(c1_infos$biological_class=="ovarian_cancer_plasma" & is.na(c1_infos$stage))] = "ovarian_cancer_plasma_ND"


c1_methyl = read.csv("/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/cg_methyl.c1.csv")
#c1_methyl = read.csv("cg_methyl.c1.csv")

c1_methyl = c1_methyl[,-1]
rownames(c1_methyl) = c1_methyl$sample
c1_methyl = c1_methyl[rownames(c1_infos),-grep("L1HS_15",colnames(c1_methyl))]
c1_methyl$biological_class = c1_infos$biological_class

# Ajout des haplotypes C1 
haplo_c1 = read.csv("/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/haplotypes.c1.csv")
rownames(haplo_c1) = haplo_c1$sample
haplo_c1 = haplo_c1[,-grep("L1HS_15",colnames(haplo_c1))]
haplo_c1 = haplo_c1[which(haplo_c1$sample %in% c1_infos$sample),]
haplo_c1 = haplo_c1[,-1]

# Ajout des eUVM de MHOVC25 en discovery
infos_all= read.csv("/groups/proudhon_lab/projects/l1pa_meth/0_data/samples_infos.csv")
#infos_all= read.csv("samples_infos.csv")
rownames(infos_all) = infos_all$sample
infos_eUVM = infos_all[which(infos_all$mhovc == "mhovc25" & infos_all$type == "plasma" & infos_all$class != "healthy" & infos_all$paper_cohort == "discovery"),]  

methyl_mhovc25 = read.csv("/groups/proudhon_lab/projects/l1pa_meth_v2/3_methylation_data_vs_all/cg_methyl.all_samples.csv")
#methyl_mhovc25 = read.csv("cg_methyl.all_samples-1.csv")
rownames(methyl_mhovc25) = methyl_mhovc25$sample
c1d_methyl=methyl_mhovc25[which(methyl_mhovc25$sample %in% infos_eUVM$sample ),]
c1d_methyl = c1d_methyl [,-grep("L1HS_15",colnames(c1d_methyl))]

haplo_mhovc25 = read.csv("/groups/proudhon_lab/projects/l1pa_meth_v2/3_methylation_data_vs_all/haplotypes.all_samples.csv")
#haplo_mhovc25 = read.csv("haplotypes.all_samples.csv")

c1d_haplo = haplo_mhovc25[which(haplo_mhovc25$sample %in% infos_eUVM$sample),]
c1d_haplo = c1d_haplo[,-grep("L1HS_15",colnames(c1d_haplo))]
rownames(c1d_haplo) = c1d_haplo$sample


# Merge de toutes les data de discovery
c1_methyl_tmp = rbind(c1_methyl, c1d_methyl)
c1_haplo =  rbind(haplo_c1, c1d_haplo)

c1_methyl = cbind(c1_methyl_tmp, c1_haplo[,3:374])


# Séparation des data en HD - cancer
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

c2_infos = read.csv("/groups/proudhon_lab/projects/l1pa_meth/0_data/c2_plasma_samples_info.csv")
#c2_infos = read.csv("c2_plasma_samples_info.csv")

rownames(c2_infos) = c2_infos$sample
c2_infos$biological_class[which(c2_infos$biological_class=="ovarian_cancer_plasma" & c2_infos$stage=="3")] = "early_ovarian_cancer_plasma"

c2_methyl = read.csv("/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/cg_methyl.c2.csv")
#c2_methyl = read.csv("cg_methyl.c2.csv")

c2_methyl = c2_methyl[,-1]
rownames(c2_methyl) = c2_methyl$sample
c2_methyl = c2_methyl[rownames(c2_infos),-grep("L1HS_15",colnames(c2_methyl))]
c2_methyl$biological_class = c2_infos$biological_class

haplo_c2 = read.csv("/groups/proudhon_lab/projects/l1pa_meth/4_methylation_data/default_scores/10_largest/haplotypes.c2.csv")
#haplo_c2 = read.csv("haplotypes.c2.csv")
haplo_c2 = haplo_c2[,-grep("L1HS_15",colnames(haplo_c2))]

c2_methyl = cbind(c2_methyl, haplo_c2[,4:375])

# Ajout des MHOVC25 
c2_infos_mhovc25 = infos_all[which(infos_all$mhovc == "mhovc25" & infos_all$type == "plasma" & infos_all$paper_cohort == "validation"),]
c2_infos_mhovc25 = c2_infos_mhovc25[-grep("D1525S028", c2_infos_mhovc25$sample),]
c2_infos_mhovc25 = c2_infos_mhovc25[-grep("D1525S020", c2_infos_mhovc25$sample),]
c2_infos_mhovc25 = c2_infos_mhovc25[-grep("D1525S021", c2_infos_mhovc25$sample),]
c2_infos_mhovc25 = c2_infos_mhovc25[-grep("D1525S038", c2_infos_mhovc25$sample),]
c2_infos_mhovc25 = c2_infos_mhovc25[-grep("D1525S039", c2_infos_mhovc25$sample),]

c2d_methyl=methyl_mhovc25[which(methyl_mhovc25$sample %in% c2_infos_mhovc25$sample ),]
c2d_methyl = c2d_methyl[,-grep("L1HS_15",colnames(c2d_methyl))]

c2d_haplo = haplo_mhovc25[which(methyl_mhovc25$sample %in% c2_infos_mhovc25$sample ),]
c2d_haplo = c2d_haplo[,-grep("L1HS_15",colnames(c2d_haplo))]

c2e_methyl = cbind(c2d_methyl, c2d_haplo[,3:374])

c2_methyl = rbind(c2_methyl, c2e_methyl)
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

