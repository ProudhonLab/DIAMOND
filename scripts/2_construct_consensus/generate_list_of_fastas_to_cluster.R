# We want to cluster only C1/C2 samples
# only plasma and no duplicates

c1 = read.csv("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/groups/proudhon_lab/projects/l1pa_meth/0_data/c1_plasma_samples_info.csv")
c2 = read.csv("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/groups/proudhon_lab/projects/l1pa_meth/0_data/c2_plasma_samples_info.csv")
c1c2 = rbind(c1[,"sample",F],c2[,"sample",F])

for (i in c(2,4,6,7,8,9,14,15)) {
  fastas_to_cluster = c1c2
  fastas_to_cluster$sample = paste("/groups/proudhon_lab/projects/l1pa_meth/1_preprocessing/preprocessing_r2/",fastas_to_cluster$sample,"/inserts_dereplicated.length_filtered.",fastas_to_cluster$sample,".L1HS_",i,".fasta",sep="")
  write.table(fastas_to_cluster,paste("/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/groups/proudhon_lab/projects/l1pa_meth/0_data/test/fastas_to_cluster_c1c2_L1HS_",i,".txt",sep=""),col.names = FALSE,row.names = FALSE, quote = FALSE)
}