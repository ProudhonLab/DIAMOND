# ! CHANGE THIS !
setwd("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/groups/proudhon_lab/projects/l1pa_meth/")
# !

library("readxl")

#####################
# Preprocessed data #
#####################

file_list = c("0_data/fastq_to_preproc_mhovc1_mergeMode.txt",
              "0_data/fastq_to_preproc_mhovc2_mergeMode.txt",
              "0_data/fastq_to_preproc_mhovc8_r2Mode.txt",
              "0_data/fastq_to_preproc_mhovc9_r2Mode.txt",
              "0_data/fastq_to_preproc_mhovc10_r2Mode.txt",
              "0_data/fastq_to_preproc_mhovc12_r2Mode.txt",
              "0_data/fastq_to_preproc_mhovc13_r2Mode.txt",
              "0_data/fastq_to_preproc_mhovc14_r2Mode.txt",
              "0_data/fastq_to_preproc_mhovc15_r2Mode.txt",
              "0_data/fastq_to_preproc_mhovc16_r2Mode.txt",
              "0_data/fastq_to_preproc_mhovc17_r2Mode.txt",
              "0_data/fastq_to_preproc_mhovc20_r2Mode.txt",
              "0_data/fastq_to_preproc_mhovc22_r2Mode.txt")

list_samples = c()
list_runcode = c()
list_mhovc = c()

for (f in file_list) {
  samples = read.table(f,header = F)
  list_mhovc = c(list_mhovc,rep(gsub(".*_(mhovc.+)_.*","\\1",f),nrow(samples)))
  list_runcode = c(list_runcode,gsub("/groups/proudhon_lab/sequencing_data/(.+)/.+","\\1",samples$V1))
  list_samples = c(list_samples,gsub("/groups/proudhon_lab/sequencing_data/.+/(.+)","\\1",samples$V1))
}

df = data.frame(sample = list_samples, mhovc = list_mhovc, run_code = list_runcode)
rownames(df) = df$sample

##########################
# Biological annotations #
##########################

# Infos from Marc's table

bioclass_annot = read.csv("0_data/samples/sample_biological-class_annotation.csv")
df$biological_class = bioclass_annot$biological_class[match(df$sample,bioclass_annot$sample)]

# Add manually remaining samples (from Book1, temporarily created by Charlotte)
# Update: Anissa discovered that L53S12 is not a methylated control but an unmethylated control
# Update 22/02/2023: Anissa discovered that V317R30 is a breast cancer cell line and L53S15 is a normal breast cell line

df[c("V317R1","V317R2","V317R3","V317R4","V317R5","V317R6","V317R7","V317R8","V317R9","V317R10","V317R31","V317R33","V317R34"),"biological_class"] = "uveal_melanoma_cellLine"
df[c("V317R11","V317R12","V317R13","V317R14","V317R15","V317R16","V317R17","V317R18","V317R19","V317R20","V317R21","V317R22","V317R23","V317R24","V317R30"),"biological_class"] = "breast_cancer_cellLine"
df[c("V317R25","V317R26","V317R27","V317R28","V317R29"),"biological_class"] = "healthy_breast_tissue"
df[c("L53S15"),"biological_class"] = "normal_breast_cellLine"
df[c("V317R32"),"biological_class"] = "normal_retin_cellLine"
df[c("V317R35","V317R36","V317R37","V317R38","V317R39","V317R40","V317R41","V317R42","V317R43","V317R44"),"biological_class"] = "healthy_plasma"
df[c("V317R45","V317R46","V317R47","L53S1"),"biological_class"] = "healthy_whiteBloodCells"
df[c("L53S2","L53S3","L53S4"),"biological_class"] = "ovarian_cancer_cellLine"
df[c("L53S5","L53S6","L53S7","L53S12","L53S14"),"biological_class"] = "unmethylated_control"
df[c("L53S8"),"biological_class"] = "colorectal_cancer_cellLine"
df[c("L53S9"),"biological_class"] = "immortalized_non_cancer_cellLine"
df[c("L53S10","L53S11","L53S13"),"biological_class"] = "fully_methylated_control"

# divide biological annotations into class and type

df$class = gsub("(.+)_[a-zA-Z]+","\\1",df$biological_class)
df$type = gsub(".+_([a-zA-Z]+)","\\1",df$biological_class)

# age
df$age = bioclass_annot$age[match(df$sample,bioclass_annot$sample)]

# sex
df$sex = bioclass_annot$sex[match(df$sample,bioclass_annot$sample)]

##########################
# Infos from SEQ_SUMMARY #
##########################
# Sheet 4 (SampleIDs) has been decomposed into separate files because its organization is not suited for analysis

df$nip = NA
df$patient_cohort = NA
df$comments = NA

# Part 1
# contains infos on NIPs about MHOVC16
  
data = read_excel("0_data/samples/seq_summary_part1.xlsx", col_names = F)
df$nip[match(data$...2,df$sample)] = data$...5

# Part 3
# contains comments about MHOVC17
# !!! mhovc17 has redo of mhovc13 == duplicates !!!

data = read_excel("0_data/samples/seq_summary_part3.xlsx", col_names = F)
df$comments[match(data$...2,df$sample)] = data$...5

# Part 2
# contains infos on NIPs and cohorts
# contains infos on MHOVC3 that I did not preprocessed / is not used now => remove from "data"
# UM A1195S380-399 are not available in the sequencing_data => remove from "data"

data = read_excel("0_data/samples/seq_summary_part2.xlsx", col_names = T)
data = data[-which(data$`seq sample` %in% setdiff(data$`seq sample`,df$sample)),]

## nips

df$nip[match(data$`seq sample`,df$sample)] = data$NIP

## cohorts

df$patient_cohort[match(data$`seq sample`[which(data$RUN=="MHOVC12")],df$sample)] = "ALCINA6"
df$patient_cohort[match(data$`seq sample`[which(data$RUN=="MHOVC13")],df$sample)] = gsub("(.+)-.+","\\1",data$`sample ID`[which(data$RUN=="MHOVC13")])
df$patient_cohort[match(data$`seq sample`[which(data$RUN=="MHOVC15")],df$sample)] = "ALCINA6"
df$patient_cohort[match(data$`seq sample`[which(data$RUN=="MHOVC16")],df$sample)] = gsub("ALC9_.+","ALCINA9",data$`sample ID`[which(data$RUN=="MHOVC16")])
df$patient_cohort[match(data$`seq sample`[which(data$RUN=="MHOVC17")],df$sample)] = gsub("ALC9.+","ALCINA9",data$`sample ID`[which(data$RUN=="MHOVC17")])

##############################
# info from discovery_cohort #
##############################
# ages, timing, localisation, stage

df$paper_cohort = NA
df$timing = NA
df$localisation = NA
df$stage = NA

# Sheet 4 (OVC)

data = read_excel("0_data/samples/DISCOVERY_COHORT.xlsx", sheet=4, col_names = T)
data = data[-which(is.na(data$sample)),]
df$paper_cohort[match(data$sample,df$sample)] = "discovery"
df$age[match(data$sample,df$sample)] = data$age
df$sex[match(data$sample,df$sample)] = data$sex
df$timing[match(data$sample,df$sample)] = data$timing
df$nip[match(data$sample,df$sample)] = data$nip
df$localisation[match(data$sample,df$sample)] = data$localisation
df$stage[match(data$sample,df$sample)] = data$STAGE...12

# Sheet 6 (BRC)

data = read_excel("0_data/samples/DISCOVERY_COHORT.xlsx", sheet=6, col_names = T)
data = data[-which(is.na(data$sample)),]
df$paper_cohort[match(data$sample,df$sample)] = "discovery"
df$age[match(data$sample,df$sample)] = data$age
df$sex[match(data$sample,df$sample)] = data$sex
df$timing[match(data$sample,df$sample)] = data$timing
df$nip[match(data$sample,df$sample)] = data$nip
df$localisation[match(data$sample,df$sample)] = data$localisation

# Sheet 7 (CRC)

data = read_excel("0_data/samples/DISCOVERY_COHORT.xlsx", sheet=7, col_names = T)
df$paper_cohort[match(data$sample,df$sample)] = "discovery"
df$age[match(data$sample,df$sample)] = data$age
df$sex[match(data$sample,df$sample)] = data$sex
df$nip[match(data$sample,df$sample)] = data$nip

# Sheet 8 (GAC)

data = read_excel("0_data/samples/DISCOVERY_COHORT.xlsx", sheet=8, col_names = T)
df$paper_cohort[match(data$sample,df$sample)] = "discovery"
df$age[match(data$sample,df$sample)] = data$age
df$sex[match(data$sample,df$sample)] = data$sex
df$nip[match(data$sample,df$sample)] = data$nip
df$stage[match(data$sample,df$sample)] = data$STAGE

# Sheet 9 (LC)

data = read_excel("0_data/samples/DISCOVERY_COHORT.xlsx", sheet=9, col_names = T)
df$paper_cohort[match(data$sample,df$sample)] = "discovery"
df$nip[match(data$sample,df$sample)] = data$nip

# Sheet 10 (UM)
# UM A1195S380-399 are not available in the sequencing_data => remove from "data"

data = read_excel("0_data/samples/DISCOVERY_COHORT.xlsx", sheet=10, col_names = T)
data = data[-which(is.na(data$sample)),]
data = data[which(data$sample %in% intersect(data$sample,df$sample)),] # où sont les A1195S380-399 ?
df$paper_cohort[match(data$sample,df$sample)] = "discovery"
df$age[match(data$sample,df$sample)] = data$age
df$sex[match(data$sample,df$sample)] = data$sex
df$timing[match(data$sample,df$sample)] = data$timing
df$nip[match(data$sample,df$sample)] = data$nip
df$localisation[match(data$sample,df$sample)] = data$localisation
df$stage[match(data$sample,df$sample)] = data$STAGE

# Sheet 11 (Healthy)
# contains MHOVC 11 that is not used now
# !!! mhovc17 has redo of mhovc13 and mhovc16 == duplicates !!!

data = read_excel("0_data/samples/DISCOVERY_COHORT.xlsx", sheet=11, col_names = T)
data = data[-which(is.na(data$sample)),]
data = data[which(data$sample %in% intersect(data$sample,df$sample)),] # on utilise pas MHOVC 11 (A1186)
df$paper_cohort[match(data$sample,df$sample)] = "discovery"
df$age[match(data$sample,df$sample)] = data$age
df$sex[match(data$sample,df$sample)] = data$sex
df$nip[match(data$sample,df$sample)] = data$nip
df$comments[match(data$sample[which(data$STAGE=="MHOVC16")],df$sample)] = "redo of MHOVC16"

########################
# info from list_c1_c2 #
########################
# contains corrected classes!!!
# contains new info age

# Sheet 1 (C1)

data = read_excel("0_data/samples/Samples_list_C1_C2_Jan10.xlsx", sheet=1, col_names = T)
data = data[-which(is.na(data$sample)),]
df$paper_cohort[match(data$sample,df$sample)] = "discovery"
# df$nip[match(data$sample,df$sample)] = data$NIP # we already have all nips here
df$stage[match(data$sample,df$sample)] = data$stage
## correcting biological class
df$biological_class[match(data$sample,df$sample)] = data$biological_class
df$class = gsub("(.+)_[a-zA-Z]+","\\1",df$biological_class)
df$type = gsub(".+_([a-zA-Z]+)","\\1",df$biological_class)
## get age
df$age[match(data$sample,df$sample)] = data$age
df$age[which(df$age=="charlotte recherche les infos exactes")] = NA
## get sex
df$sex[match(data$sample,df$sample)] = data$sex

# Sheet 2 (C2)

data = read_excel("0_data/samples/Samples_list_C1_C2_Jan10.xlsx", sheet=2, col_names = T)
data = data[-which(is.na(data$sample)),]
df$paper_cohort[match(data$sample,df$sample)] = "validation"
df$nip[match(data$sample,df$sample)] = data$NIP
df$patient_cohort[match(data$sample,df$sample)] = data$Cohorte
df$stage[match(data$sample,df$sample)] = data$stage
df$age[match(data$sample,df$sample)] = data$`age à l'inclusion`
df$sex[match(data$sample,df$sample)] = data$sex

#########################
# info from Table S1-S3 #
#########################
# Last infos, from the Tables of the paper Michel et al
# adding missing Gastric age

data = read_excel("0_data/samples/Tables_S1-S3.xlsx", sheet=3, col_names = T)
colnames(data) = data[2,]
data = data[-c(1,2,nrow(data)-1,nrow(data)),]
df$age[match(data$Sample_ID,df$sample)] = data$Age

# update 11/05/2023: adding updated sex and stages infos
df$sex[match(data$Sample_ID,df$sample)] = data$Sex
df$stage[match(data$Sample_ID,df$sample)] = data$Stage

#####################
# Number of inserts #
#####################
# useful to remove samples with not enough inserts for analyses

df$nb_inserts = NA

# deduplicated inserts (r2 mode)

data = read.csv("1_preprocessing/preprocessing_r2/allSamplesR2.summary.csv")
df$nb_inserts[match(data$sample,df$sample)] = data$deduplicated.length.filtered.inserts

# inserts (merge mode, no UMIs)
data = read.csv("1_preprocessing/preprocessing_merge/allSamplesMerge.summary.csv")
df$nb_inserts[match(data$sample,df$sample)] = data$length.filtered.inserts

#####################
# Final adjustments #
#####################

# Anissa confirmed D541S056 is T1 and D541S057 is T2
df$timing[which(df$sample=="D541S056")] = "T1"
df$timing[which(df$sample=="D541S057")] = "T2"

df$mhovc = factor(df$mhovc)
df$run_code = factor(df$run_code)
df$biological_class = factor(df$biological_class)
df$age[which(df$age=="NA")] = NA
df$age = floor(as.numeric(df$age))
df$sex[which(df$sex=="NA")] = NA
df$sex[which(df$sex=="H")] = "M"
df$sex = factor(df$sex)
df$nip[which(df$nip=="?")] = NA
df$nip[which(df$nip=="NA")] = NA
df$patient_cohort = factor(df$patient_cohort)
df$paper_cohort = factor(df$paper_cohort)
df$localisation = factor(df$localisation)
df$stage[which(df$stage=="NA")] = NA
df$stage[which(df$stage=="ND")] = NA

df$comments[which(df$nip=="1206209")] = "BRCA1 mutation, maybe associated with ovarian cancer or even OVC+eBRC"

write.csv(df,"0_data/samples_infos.csv",row.names = FALSE)


###############################
# Only C1 C2 plasma for paper #
###############################
#  (note: D541S061 and D541S062 are the same patient but different cancer)

c1 = droplevels(df[which(df$type=="plasma" & df$paper_cohort=="discovery"),])
c2 = droplevels(df[which(df$type=="plasma" & df$paper_cohort=="validation"),])

# C1

## remove MHOVC15, not enough reads
c1 = c1[-which(c1$mhovc=="mhovc15"),]

## breast correction
## remove A1190S193, reason unknown but confirmed with Anissa
c1 = c1[-which(c1$sample=="A1190S193"),]

## colorectal correction
## remove A1179S50 because not enough reads
c1 = c1[-which(c1$sample=="A1179S50"),]

## early breast correction
## D541S056 and D541S057 are duplicates, keep T1
c1 = c1[-which(c1$sample=="D541S057"),]
## remove D541S038 because not enough reads
c1 = c1[-which(c1$sample=="D541S038"),]
## remove the following samples, unknown reason
c1 = c1[-which(c1$sample %in% c("D541S033","D541S045","D541S049","D541S051","D541S055","D541S063")),]
## update 01/03/2023: several samples are from patients treated with chimio => need to be excluded
c1 = c1[-which(c1$sample %in% c("D541S034","D541S036","D541S044","D541S058","D541S062","D541S065","D541S066","D541S067")),]

## lung duplicates correction
## Duplicates rules defined by Charlotte: keep D541S018 D541S021, potentially the first samples (to confirm) == remove D541S023 D541S026
c1 = c1[-which(c1$sample=="D541S023" | c1$sample=="D541S026"),]

## ovarian duplicates correction
## keep A1146S1-23 == remove ovarian_cancer_plasma from mhovc14
c1 = c1[-which(c1$mhovc=="mhovc14" & c1$biological_class=="ovarian_cancer_plasma"),]
## keep only T1
c1 = c1[-which(c1$timing=="T2" & c1$biological_class=="ovarian_cancer_plasma"),]
## remove OVC_23 because bad sample (reason unknown), no NIP
c1 = c1[-which(c1$nip=="OVC_23"),]
## remove (TEMPORARELY!!!) D541S061 because Anissa needs more info about it
c1 = c1[-which(c1$sample=="D541S061"),]

## uveal melanoma duplicates correction
## depends on the cohort
## importance defined by Charlotte: prioritize MU 2010-01 Cohorte 12, then MU 2010-02, then T1 between two MU 2012-01 (Novartis)
duplicates_to_remove = c("D473S459","D473S450","D473S447","D473S460","D473S449","D473S446","D473S444","D473S453","D473S454","D473S488","D473S462","D473S451","D473S445","D473S452","D473S438")
c1 = c1[-which(c1$sample %in% duplicates_to_remove),]
## we added new samples from MHOVC15, but they might already exist in MHOVC16 = keep MHOVC16 only / update: no more MHOVC 15
## nip 1084538 has T2 and T3 = keep T2, remove A1195S378
#duplicates_to_remove = c("A1195S378")
## nip 1104936, 1201023, 883411, 1106086, 1213273, 781449, 883411, 988679, 990720, 980730 already in MHOVC16
#duplicates_to_remove = c(duplicates_to_remove,"A1195S376","A1195S413","A1195S379","A1195S406","A1195S407","A1195S411","A1195S379","A1195S401","A1195S400","A1195S410")
## nip 992020 has two ALCINA6 with same timing ?
#duplicates_to_remove = c(duplicates_to_remove,"A1195S402","A1195S403")
#c1 = c1[-which(c1$sample %in% duplicates_to_remove),]

## Anissa asked to add (december 14th, 2022):
## 781653 -> OK
## 683857 -> OK
## 1086127 -> but missing A1195S380
## 785778 -> A1195S383
## 1119595 -> A1195S384
## 1212356 -> A1195S386
## 9481125 -> A1195S389
## 506322 -> A1195S390
## 488360 -> A1195S395
## 1213228 -> A1195S396
## 1213249 -> A1195S397
## 1092277 -> OK
## 880688 -> OK
## 1106189 -> OK
## 1215435 -> OK
## 1084538 -> duplicates OK
## 580507 -> ?? already in MHOVC16 and missing A1195S381 and A1195S399
## 992020 -> ?? already in MHOVC16 + in MHOVC15, same cohort same timing

## healthy plasma duplicates correction
## remove the redos
c1 = c1[-which(c1$comments=="redo from MH (MHOVC13)" | c1$comments == "redo of MHOVC16"),]
## remove mhovc 8, 10, 12
c1 = c1[-which(c1$biological_class=="healthy_plasma" & (c1$mhovc == "mhovc8"|c1$mhovc=="mhovc10"|c1$mhovc=="mhovc12")),]
## remove bad sample D473S427
c1 = c1[-which(c1$sample=="D473S427"),]

write.csv(c1,"0_data/c1_plasma_samples_info.csv",row.names = FALSE)

# C2

# old correction (from list C1C2 Dec22):
## remove D1000S080, not enough reads
#c2 = c2[-which(c2$sample=="D1000S080"),]
## remove the Qias vs Manual test samples
#to_remove = c("D1158S57","D1158S58","D1158S59","D1158S60","D1158S61","D1158S62","D1158S63","D1158S64","D1158S65","D1158S66")
#c2 = c2[-which(c2$sample %in% to_remove),]

write.csv(c2,"0_data/c2_plasma_samples_info.csv",row.names = FALSE)

##########################
# Only tissues for paper #
##########################

tissues = droplevels(df[which(df$type=="tissue"|df$type=="control"|df$type=="whiteBloodCells"),])

# remove bad sample L53S5 (because how the KO was done)
tissues = tissues[-which(tissues$sample=="L53S5"),]

write.csv(tissues,"0_data/tissues_samples_info.csv")

#############################
# Only cell lines for paper #
#############################

cell_lines = droplevels(df[which(df$type=="cellLine"|df$type=="control"|df$type=="whiteBloodCells"),])

write.csv(cell_lines,"0_data/cellLines_samples_info.csv")
