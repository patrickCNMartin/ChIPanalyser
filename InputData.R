## Creating new input data for ChIPanalyser
#!/usr/bin/Rscript

### Loading Libraries and Scripts

direc<-getwd()

library(BSgenome.Dmelanogaster.UCSC.dm6)
library(BSgenome)
library(RcppRoll)
library(GenomicRanges)
library(ROCR)


setwd("/home/patrickmartin/ChIPanalyser/ChIPanalyserFinal/ChIPdev")
files <- dir()
for (i in files) source(i)

setwd(direc)
##
load("BG3_modEncode_921_BEAF-32_5cat_reduce100MSEKingOfTheHilloutputValidation.Rda")


#### top 10 Granges in that data set
top <- lapply(KingOfTheHill[[1]]@profiles[[1]][1:4], function(k){
      gr <- GRanges(seqnames = as.character(seqnames(k))[1L],
                    ranges = IRanges(start = start(k)[1L], end = end(k)[length(k)]))
})
top <- unlist(as(top, "GRangesList"))



#### Getting Chip profile data for each of those regions
chip <- import("/home/patrickmartin/ChIP/modEncode_BG3/modEncode/modEncode_921/signal_data_files/BEAF-32:Cell-Line=ML-DmBG3-c2#Developmental-Stage=Larvae-3rd-instar#Tissue=CNS-derived-cell-line:ChIP-chip:Rep-1::Dmel_r5.32:modENCODE_921:repset.4620429.smoothedM.bed")
chip <- subsetByOverlaps(chip, top)

#### Getting CS for those ranges
cs <- import("/home/patrickmartin/ChIPanalyser/ChIPanalyserFinal/CS/data/BG3-11states.bed.gz")
cs <- subsetByOverlaps(cs, top)

#### Getting Access
Access <- get(load("/home/patrickmartin/DNAaccess/cellAccess/BG3_DHS_005.Rda"))
Access <- subsetByOverlaps(Access,top)

save(top,chip,cs,Access, file = "/home/patrickmartin/ChIPanalyser/ChIPanalyserData.Rda")
