setwd("/Users/pamalamilla/Desktop/Fragmentomics_GOOD/")

library(stringr)
library(dplyr)

#Limma----
lymphoma_vs_control_limma = read.table("2024_10_09_lymphoma_vs_control_limma.txt",
                                 header = T)
nrow(lymphoma_vs_control_limma) #30,505 windows

### p-val <= 0.05 ----
limma_p0.05 = filter(lymphoma_vs_control_limma, adj.P.Val <= 0.05)
nrow(limma_p0.05) #19,019 windows with p-val <= 0.05

####logFC >= 0.585 -------------------------------------------------------------
#Represents a ~50% (or greater) increase in methylation in lymphoma
limma_p0.05_logFC0.585 = filter(limma_p0.05, logFC >= 0.585)
nrow(limma_p0.05_logFC0.585) #17,631 windows with logFC >= 0.585

#Output BED file
limma_p0.05_FC1.5 = limma_p0.05_logFC0.585$windows %>% unlist() %>% 
  sort() %>% str_replace_all("\\."," ")
write.table(limma_p0.05_FC1.5, file = "limma_p0.05_FC1.5.bed", 
            quote = F, row.names = F, col.names = F)

####logFC <= -1 ----------------------------------------------------------------
#Represents a 50% (or greater) decrease in methylation in lymphoma)
#limma_p0.05_logFCneg1 = filter(limma_p0.05, logFC <= -1)
#nrow(limma_p0.05_logFCneg1) #707 windows with logFC <= -1
