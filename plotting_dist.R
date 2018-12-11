working_dir <- "/Users/alanaalexander"
filename <- "dist_between_samples.txt"

library(readr)
library(dplyr)
library(ggplot2)

setwd(working_dir)  
temp <- read_table2(filename)

sample1col <- temp$sample1
sample2col <- temp$sample2
temp2 <- temp
temp2
temp <- rbind(temp,temp2)

temp <- temp %>% mutate(sampleplot=gsub("\\.[1-2]","",sample1)) %>% arrange(sampleplot)

temp <- temp %>% mutate(comparison=ifelse(gsub("\\.[1-2]","",sample1)==gsub("\\.[1-2]","",sample2),"within","between"))

temp <- temp %>% mutate(toplot=paste(sampleplot,"_",comparison,sep="")) %>% arrange(toplot)

between_only <- temp %>% filter(comparison=="between")



ggplot(data = between_only,mapping=aes(x = sampleplot,y=standtolocusmax)) + 
  geom_violin() + geom_jitter()
