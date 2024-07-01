#!/usr/bin/env Rscript


#@@@@@@@@@@@@@ LOADING LIBRAIRIES @@@@@@@@@@@@@
if (!require("pacman"))
  install.packages("pacman") #install pacman if not installed already on the conda env and then load it

pacman::p_load("tidyverse",
               "glue" ) #way faster load via pacman of the packages

#@@@@@@@@@@@@@ FUNCTION @@@@@@@@@@@@@

ZscoreCalulating = function(DataNorm,  Amp, ref) {
  dfSample = suppressMessages(DataNorm %>% filter( amp == Amp) %>% select( chr:reads, normReads, amp,sample) %>%
                                left_join(ref, by = c("chr", "start", "end", "amp")
                                ) %>% mutate(zscore = (normReads - MeanC) / sdC) %>% mutate(
                                  colour = case_when(zscore <= -5 ~ "blue", zscore >= 5 ~ "red", TRUE ~ "darkolivegreen2"),
                                  zscore = case_when(is.na(zscore) ~ 0, TRUE ~ zscore)
                                ))
  assign(x = Amp
         ,
         value = dfSample,
         envir = .GlobalEnv)
}


#@@@@@@@@@@@@@ LOADING ARGUMENTS @@@@@@@@@@@@@
args = commandArgs(trailingOnly = TRUE)
print(args)
input = args[1]
output = args[2]
aggregated_Amp= args[3]
chromosomeArmExclude=unlist(strsplit(args[4],","))


#@@@@@@@@@@@@@ LOADING FILES @@@@@@@@@@@@@
data <- read.csv(input, header = FALSE,sep = ";") 

# "chr1" "0" "155527" "R2VS45" "L1HS_2" "1558" "Exp"
# "chr1" 0 123400000 561 "L53S10" "L1HS_2" "Exp"

#@@@@@@@@@@@@@ MAIN CODE @@@@@@@@@@@@@

names(data) = c("chr", "start", "end", "reads", "sample", "amp",  "group") 


if(aggregated_Amp==TRUE){
  aggregatedData= data%>%group_by(chr,start,end,sample)%>%mutate(reads2=sum(reads))%>%select(-reads)%>%dplyr::rename("reads"="reads2")%>%ungroup()%>%mutate(amp="Mixed")%>%distinct()
  data= aggregatedData
}

DataNorm = data%>%mutate(chrdetails = case_when(start == 0 ~ paste0(chr, "p"), TRUE ~ paste0(chr, "q"))) %>%filter(chr!="chrY",chr!="chrX",!(chrdetails %in% chromosomeArmExclude) )%>% group_by(sample, amp) %>% mutate(SumampSample = sum(reads)) %>% #change maybe
  ungroup() %>% mutate(normReads = ((reads / SumampSample)))

RefData = DataNorm %>% filter(group == TRUE) %>% group_by(amp, start, end, chr) %>%
  mutate(MeanC = mean(normReads), sdC = sd(normReads)) %>% ungroup() %>% select(chr, start, end, sdC, MeanC, amp,) %>%
  distinct()

ExpData = DataNorm %>% filter(group == FALSE)
zscoreByChrArm =   suppressMessages(ExpData %>% group_by( amp ) %>% select( chr:reads, normReads, amp,sample) %>%
                                                left_join(RefData, by = c("chr", "start", "end", "amp")
                                                ) %>% mutate(zscore = ifelse(test = is.na((normReads - MeanC) / sdC ),0,(normReads - MeanC) / sdC ))%>%mutate(Gzscore=sum((zscore)^2))%>%ungroup())
  

print(names(zscoreByChrArm))
write.table(zscoreByChrArm%>%select(Gzscore,sample)%>%distinct()%>%arrange(Gzscore),glue("{output}_Genome_wideZscoreByChrArm.csv"),col.names=FALSE,row.names=FALSE,sep=";")
write.table(zscoreByChrArm%>%select(!c(Gzscore,group,reads,normReads,sdC,MeanC)),glue("{output}_ZscoreByChrArm.csv"),col.names=FALSE,row.names=FALSE,sep=";")

