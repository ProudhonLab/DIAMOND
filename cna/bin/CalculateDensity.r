#!/usr/bin/env Rscript


#@@@@@@@@@@@@@ LOADING LIBRAIRIES @@@@@@@@@@@@@
if (!require("pacman"))
  install.packages("pacman")
#install pacman if not installed already on the conda env and then load it
pacman::p_load("tidyverse",
               "karyoploteR",
               "glue") #way faster load via pacman of the packages

#@@@@@@@@@@@@@ LOADING ARGUMENTS @@@@@@@@@@@@@
args = commandArgs(trailingOnly = TRUE)

sample = args[1]
Amp= args[2]
data = args[3]
ctrlStatus = args[4]
ArmOrInterval = args[5]
print(ArmOrInterval)
print(ArmOrInterval==TRUE)
if(ArmOrInterval){
  chrPos = args[6]
}else{
  window.size = args[6]
}

#@@@@@@@@@@@@@ MAIN CODE @@@@@@@@@@@@@

x <- rtracklayer::import.bed(data)

kpWin = plotKaryotype(
  plot.type = 4,
  ideogram.plotter = NULL,
  labels.plotter = NULL,
  genome = "hg38"
)
if (ArmOrInterval == FALSE) {
  windows = tileGenome(
    setNames(kpWin$chromosome.lengths, names(kpWin$chromosome.lengths)),
    tilewidth = window.size,
    cut.last.tile.in.chrom = TRUE
  )
  posdataframe = data.frame(
    chr = as.character(seqnames(windows)),
    start = as.numeric(start(windows)),
    end = as.numeric(end(windows))
  )
} else{
  posdataframe = read.delim(chrPos,
    sep = ";"
  ) 
  windows = makeGRangesFromDataFrame(posdataframe)
  
}

countCNA = cbind(posdataframe,
                 countOverlaps(windows, x),sample,Amp,ctrlStatus)

write.table(countCNA,
            paste0(sample,"_L1HS_",Amp,
                    ".csv"),sep = ";",
            col.names = FALSE,row.names = FALSE)
