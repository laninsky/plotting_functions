working_dir <- "/Users/alanaalexander/Downloads"
distfilename <- "dist_between_samples.txt"
covfilename <- "coverage_summary.txt"

# The dist_between_samples.txt results from https://github.com/laninsky/distance_calcs
# The coverage_summary.txt results from https://github.com/laninsky/reference_aligning_to_established_loci/tree/master/phase_target

# Loading in libraries we need  
library(readr)

# Setting the variables to our bits and pieces  
setwd(working_dir)  
disttemp <- read_table2(distfilename)

# Doing some stuff on disttemp
disttemp
#locus                            sample1    sample2    raw_dist standtolocusmax
#<chr>                            <chr>      <chr>         <dbl>           <dbl>
#  1 locus_specific_fasta/1.fasta.out amphizoa.1 amphizoa.1      0             0    
#2 locus_specific_fasta/1.fasta.out amphizoa.1 amphizoa.2      0             0    
#3 locus_specific_fasta/1.fasta.out amphizoa.1 chlSer1.1      17.8           0.148
#4 locus_specific_fasta/1.fasta.out amphizoa.1 chlSer1.2      17.8           0.148
#5 locus_specific_fasta/1.fasta.out amphizoa.1 pterMel1.1    102.            0.846
#6 locus_specific_fasta/1.fasta.out amphizoa.1 pterMel1.2    102.            0.844
#7 locus_specific_fasta/1.fasta.out amphizoa.1 traGib1.1      17.4           0.145
#8 locus_specific_fasta/1.fasta.out amphizoa.1 traGib1.2      17.4           0.145
#9 locus_specific_fasta/1.fasta.out amphizoa.2 amphizoa.2      0             0    
#10 locus_specific_fasta/1.fasta.out amphizoa.2 chlSer1.1      17.8           0.148
# ... with 2,685 more rows

# This bit is so every comparison involving a given sample for a given locus is present with the same name in the sample1 column
sample1col <- disttemp$sample1
sample2col <- disttemp$sample2

disttemp2 <- disttemp
disttemp2$sample1 <- sample2col
disttemp2$sample2 <- sample1col  

disttemp <- rbind(disttemp,disttemp2)

# Pulling out the sample locus so it is equivalent to covtemp
disttemp <- disttemp %>% mutate(real_locus=gsub(".fasta.out","",gsub("locus_specific_fasta/","",locus)))

# Adding a column specifying if the comparison is between different samples or within (e.g. comparing alleles of the same genome)
disttemp <- disttemp %>% mutate(comparison=ifelse(gsub("\\.[1-2]","",sample1)==gsub("\\.[1-2]","",sample2),"within","between"))

# Adding a column that gets rid of the allele *.1 and *.2 for each sample
disttemp <- disttemp %>% mutate(sampleplot1=gsub("\\.[1-2]","",sample1))

# Summarizing distances on a per-locus, per-sample basis so we'll be able to join it with covtemp
sumdisttemp <- disttemp %>%  filter(comparison=="between") %>% group_by(sampleplot1,real_locus) %>% summarise(meanrawdist=mean(raw_dist,na.rm=TRUE),meanstddist=mean(standtolocusmax,na.rm=TRUE))

# Reading in the coverage information
covtemp<- read_table2(covfilename)

covtemp
## A tibble: 239 x 8
#samplename locus `ref_length(bp)` bp_covered_by_seq min_cov max_cov mean_inc_0 mean_exc_0
#<chr>      <int>            <int>             <int>   <int>   <int>      <dbl>      <dbl>
#  1 amphizoa      10               88                88       2       9       7.28       7.28
#2 amphizoa      11               60                60       3       6       5.45       5.45
#3 amphizoa      12              236               236       2       9       6.03       6.03
#4 amphizoa      13              142               142       2       5       3.79       3.79
#5 amphizoa      15              177               177       2       9       4.04       4.04
#6 amphizoa      16              129               129       4      11       8.50       8.50
#7 amphizoa       1               65                65       2       4       3.22       3.22
#8 amphizoa      20              298               298       3       8       5.54       5.54
#9 amphizoa      21              204               204       1       7       4.28       4.28
#10 amphizoa      22               40                40       3       5       4.85       4.85
# ... with 229 more rows

# Renaming some of column names so we can join the tibbles
names(covtemp)[1] <- "sampleplot1"
names(covtemp)[2] <- "real_locus"

# Converting the real_locus field of covtemp to a character so we can join tables
covtemp <- covtemp %>% mutate(real_locus=as.character(real_locus))

# anti_join will return any rows that can't be matched up. If they can be matched (e.g. the tibble is 0) up we are good to go.
anti_join(sumdisttemp,covtemp)
# Ones returned here are "singletons" that were not found across multiple species, or have 0 coverage
anti_join(covtemp,sumdisttemp)

covdist <- inner_join(sumdisttemp,covtemp)

covdist
#sampleplot1 real_locus meanrawdist meanstddist `ref_length(bp)` bp_covered_by_seq min_cov max_cov mean_inc_0 mean_exc_0
#<chr>       <chr>            <dbl>       <dbl>            <int>             <int>   <int>   <int>      <dbl>      <dbl>
#  1 amphizoa    1                 45.7       0.379               65                65       2       4       3.22       3.22
#2 amphizoa    10                51.4       0.444               88                88       2       9       7.28       7.28
#3 amphizoa    11                29.8       0.235               60                60       3       6       5.45       5.45
#4 amphizoa    12                16.4       0.606              236               236       2       9       6.03       6.03
#5 amphizoa    13               102.        0.716              142               142       2       5       3.79       3.79
#6 amphizoa    15                14.9       0.422              177               177       2       9       4.04       4.04
#7 amphizoa    16                20.2       0.652              129               129       4      11       8.50       8.50
#8 amphizoa    21                12.5       0.635              204               204       1       7       4.28       4.28
#9 amphizoa    22                79.2       0.827               40                40       3       5       4.85       4.85
#10 amphizoa    27                97.7       0.727              104               104       1       8       6.41       6.41


# Use this first inspection of data to pare down dataset to fewer variables
pairs(covdist[3:10])

# Paring down
paireddowncovdist <- covdist %>% select(sampleplot1,meanrawdist,meanstddist,bp_covered_by_seq,min_cov,max_cov)

# Double checking paring worked ok
pairs(paireddowncovdist[2:6])

# Going to get our data in order
filtertemp <- paireddowncovdist %>% 
  group_by(sampleplot1) %>% 
  summarise(meanstand=mean(meanstddist,na.rm=TRUE)) %>% 
  arrange(meanstand)

# Ordering genomes on the x-axis by mean standardized distance
paireddowncovdist$sampleplot1 <- factor(paireddowncovdist$sampleplot1, levels = filtertemp$sampleplot1)

# Renaming our factors so they look OK when plotted
names(paireddowncovdist) <- c("Genome","Mean raw dist.","Mean std. dist.","Locus length (bp)","Min. cov.","Max. cov.")

# Loading in GGally which is going to allow us to use the kick-ass ggpairs function
library(GGally)

#Combo options
#box: boxplots
#box_no_facet: listing samplenames beneath
#dot_no_facet: listing samplenames beneath
#facethist: histograms, listing y-scale beneath
# facetdensity: single density curve per sample line: listing y-scale beneath
# denstrip: summmarizes density - somewhat similar to boxplot

#continuous options
# density: density curves
# smooth: fits line with error shades
# smooth_loess: fits crazy curvy line with error shading

ggpairs(paireddowncovdist, mapping=(aes(colour = Genome, alpha = 0.4)),columns=c("Genome","Mean raw dist.","Mean std. dist.","Locus length (bp)","Min. cov.","Max. cov."),upper=list(continuous="density"),lower=list(combo="dot"))

