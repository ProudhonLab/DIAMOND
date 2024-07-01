library(ggplot2)

workDir="/run/user/1001/gvfs/sftp:host=genossh.genouest.org,user=kdasilva/scratch/kdasilva/l1pameth_2023/consensus_sequences_abund/"

for (primer in c("L1HS_2","L1HS_4","L1HS_6","L1HS_7","L1HS_8","L1HS_9","L1HS_14","L1HS_15")) {
  
  data = read.csv(paste(workDir,"c1c2_",primer,".consensus.nb_seqs.csv",sep=""), header = F)
  data$cluster = factor(seq(1,nrow(data)),levels=seq(1,nrow(data)))
  colnames(data)[1] = "Raw number of sequences"
  data["Proportion of total sequences"] = data$`Raw number of sequences`/sum(data$`Raw number of sequences`)*100
  data$cumsum=cumsum(data$`Proportion of total sequences`)
  
  g1 = ggplot(data.frame(data[1:10,]),aes(x=cluster,y=Raw.number.of.sequences)) +
  geom_bar(stat="identity")
  g2 = ggplot(data.frame(data[1:10,]),aes(x=cluster,y=Proportion.of.total.sequences)) +
  geom_bar(stat="identity") +labs(title=primer)
  ggsave(paste("/home/kevin/sharedVM/clustering_exploration/barplot_proportion_sequences_cluster_",primer,".png",sep=""),g2)

}
