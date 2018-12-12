working_dir <- "/Users/alanaalexander"
filename <- "dist_between_samples.txt"

# The dist_between_samples.txt results from https://github.com/laninsky/distance_calcs

# Plotting the distance functions
library(readr)
library(dplyr)
library(ggplot2)

setwd(working_dir)  
temp <- read_table2(filename)
## A tibble: 2,695 x 5
#   locus                            sample1    sample2    raw_dist standtolocusmax
#   <chr>                            <chr>      <chr>         <dbl>           <dbl>
# 1 locus_specific_fasta/1.fasta.out amphizoa.1 amphizoa.1      0             0    
# 2 locus_specific_fasta/1.fasta.out amphizoa.1 amphizoa.2      0             0    
# 3 locus_specific_fasta/1.fasta.out amphizoa.1 chlSer1.1      17.8           0.148
# 4 locus_specific_fasta/1.fasta.out amphizoa.1 chlSer1.2      17.8           0.148
# 5 locus_specific_fasta/1.fasta.out amphizoa.1 pterMel1.1    102.            0.846
# 6 locus_specific_fasta/1.fasta.out amphizoa.1 pterMel1.2    102.            0.844
# 7 locus_specific_fasta/1.fasta.out amphizoa.1 traGib1.1      17.4           0.145
# 8 locus_specific_fasta/1.fasta.out amphizoa.1 traGib1.2      17.4           0.145
# 9 locus_specific_fasta/1.fasta.out amphizoa.2 amphizoa.2      0             0    
#10 locus_specific_fasta/1.fasta.out amphizoa.2 chlSer1.1      17.8           0.148
## ... with 2,685 more rows

sample1col <- temp$sample1
sample2col <- temp$sample2

temp2 <- temp
temp2$sample1 <- sample2col
temp2$sample2 <- sample1col  
  
temp <- rbind(temp,temp2)

temp <- temp %>% mutate(sampleplot1=gsub("\\.[1-2]","",sample1))
temp <- temp %>% mutate(sampleplot2=gsub("\\.[1-2]","",sample2))
temp <- temp %>% mutate(comparison=ifelse(gsub("\\.[1-2]","",sample1)==gsub("\\.[1-2]","",sample2),"within","between"))
temp <- temp %>% mutate(toplot=paste(sampleplot1,"_",comparison,sep="")) %>% arrange(toplot)

# Potentially having only two different taxa present (e.g. temp) penalizes samples who have more loci present
# b/c when only two taxa are present, they are either 0 or 1
morethantwo <- temp %>% group_by(locus) %>% count(sampleplot1) %>% count(locus) %>% filter(nn > 2) %>% select(locus)
temp2 <- temp %>% filter(locus %in% morethantwo$locus)

# Same rational but filtering to loci just present across all individuals
morethansix <- temp %>% group_by(locus) %>% count(sampleplot1) %>% count(locus) %>% filter(nn > 6) %>% select(locus)
temp6 <- temp %>% filter(locus %in% morethansix$locus)

# Setting the dataframe we'll plot (could be temp, temp2, temp6)
datatoplot <- temp6

# Calculating error bars for standardized distances
filtertemp <- datatoplot %>% 
  filter(comparison=="between") %>% 
  group_by(sampleplot1) %>% 
  summarise(meanstand=mean(standtolocusmax,na.rm=TRUE),standstand=sd(standtolocusmax,na.rm=TRUE),standcount=n(),ci=((standstand/sqrt(standcount))*1.96)) %>% 
  arrange(meanstand)

# Ordering genomes on the x-axis by mean standardized distance
datatoplot$sampleplot1 <- factor(datatoplot$sampleplot1, levels = filtertemp$sampleplot1)

# Shoving together everything to do with standardized distance into a single plot
# Note: need to tweak the title if different datasets (e.g. temp, temp2, temp6) used
databoxplotstand <- ggplot() + 
  geom_jitter(data = (datatoplot %>% filter(comparison=="between")),alpha=0.5,mapping=aes(x=sampleplot1,y=standtolocusmax)) +
  geom_boxplot(data = (datatoplot %>% filter(comparison=="between")),alpha=0.5,color="blue",mapping=aes(x=sampleplot1,y=standtolocusmax)) + 
  geom_errorbar(data=filtertemp,aes(x=reorder(sampleplot1,meanstand),ymin=meanstand-ci, ymax=meanstand+ci), width=.2,size=1) +
  geom_point(data=filtertemp,aes(x=reorder(sampleplot1,meanstand),y=meanstand), size=2) +
  theme_bw() + 
  theme(legend.position="none") +
  xlab("Genome") +
  ylab("Genetic distance") +
  ggtitle("Loci found across all seven species",subtitle="All between-taxa pairwise comparisons for all loci, standardized to maximum distance per locus") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5))

# Printing graph to screen to save
databoxplotstand

# Doing the same steps above except with the raw rather than standardized distances
filtertemp_raw <- datatoplot %>% 
  filter(comparison=="between") %>% 
  filter(is.finite(raw_dist)) %>% 
  group_by(sampleplot1) %>% 
  summarise(meanstand=mean(raw_dist,na.rm=TRUE),standstand=sd(raw_dist,na.rm=TRUE),standcount=n(),ci=((standstand/sqrt(standcount))*1.96)) %>% 
  arrange(meanstand)
  
  datatoplot$sampleplot1 <- factor(datatoplot$sampleplot1, levels = filtertemp_raw$sampleplot1)  

  databoxplotraw <- ggplot() + 
  geom_jitter(data = (datatoplot %>% filter(comparison=="between")),alpha=0.5,mapping=aes(x=sampleplot1,y=raw_dist)) +
  geom_boxplot(data = (datatoplot %>% filter(comparison=="between")),alpha=0.5,color="blue",mapping=aes(x=sampleplot1,y=raw_dist)) + 
  geom_errorbar(data=filtertemp_raw,aes(x=reorder(sampleplot1,meanstand),ymin=meanstand-ci, ymax=meanstand+ci), width=.2,size=1) +
  geom_point(data=filtertemp_raw,aes(x=reorder(sampleplot1,meanstand),y=meanstand), size=2) +
  theme_bw() + 
  theme(legend.position="none") +
  xlab("Genome") +
  ylab("Genetic distance") +
  ggtitle("Loci found across two or more species",subtitle="All between-taxa pairwise comparisons for all loci, raw distances per locus") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5))

  databoxplotraw
