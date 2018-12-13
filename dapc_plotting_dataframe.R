#vDefine where your working directory is
# Double slashes if you are in windows
#working_dir <- "C:\\Users\\Alana\\Dropbox\\2017_data\\Transects"
# Single if you are on a mac
working_dir <- "/Users/alanaalexander/Dropbox/2017_data/Transects"
  
#Ctrl+Enter to send command down to console

# Loading libraries that have some packages we'll want
# If not already installed you'll need to install them first
# e.g. install.libraries("readr")
library(readr)
library(readxl)
library(dplyr)
library(ggplot2)
library(tibble)
library(factoextra)
library(ade4)
library(plot3D)
library(plot3Drgl)
library(RColorBrewer)
library(adegenet)

# Setting working directory
setwd(working_dir)

# Getting a list of files
files <- list.files(pattern=".xlsx")

# creating the fulldata object that we'll store the excel sheets in
fulldata <- NULL

# For each of the excel workbooks
for (i in files) {
  # Temp holds each sheet for workbook i
  temp <- excel_sheets(i)
  for (j in temp) {
    # For each sheet j from workbook i
    temptib <- read_excel(i, sheet = j)
    # Only including the actual data and not the groundtruth data
    temptib <- temptib[,1:14]
    # Using gsub to extract the site name from the workbook name
    sitename <- gsub("\\(transects\\).xlsx","",(gsub("2017_","",i)))
    # Adding a sample and site column to our data
    temptib <- add_column(temptib, sample = j, .before = 1)
    temptib <- add_column(temptib, site = sitename, .before = 1)
    # Binding the data for this worksheet into the fulldata object
    fulldata <- rbind(fulldata,temptib)
  }
}

# Making an individual row identification
indcol <- fulldata %>% transmute(indcol=paste(site,"_",sample,"_",ElapsedTime_s,sep=""))

# Adding the individual row identification
fulldata <- as_tibble(cbind(fulldata[1:4],indcol, fulldata[5:16]))

fulldata
# A tibble: 23,980 x 17
#   site  sample ElapsedTime_s `Distance (µm)` indcol Ca43_CPS `Li_mol/molCa`     Na
#   <chr> <chr>          <dbl>           <dbl> <chr>     <dbl>          <dbl>  <dbl>
# 1 Site1 E3             0               0     Site1…   11563.     0          0     
# 2 Site1 E3             0.183           0.128 Site1…   19597.     0          0     
# 3 Site1 E3             0.366           0.256 Site1…   14174.     0          0     
# 4 Site1 E3             0.549           0.384 Site1…   13334.     0.0000135  0.0170
# 5 Site1 E3             0.732           0.512 Site1…   14492.     0.00000532 0.0160
# 6 Site1 E3             0.915           0.641 Site1…   15090.     0.00000398 0.0158
# 7 Site1 E3             1.10            0.769 Site1…   15876.     0.00000702 0.0196
# 8 Site1 E3             1.28            0.897 Site1…   13546.     0.0000120  0.0159
# 9 Site1 E3             1.46            1.02  Site1…   14249.     0.00000662 0.0158
#10 Site1 E3             1.65            1.15  Site1…   15181.     0.00000508 0.0157
## ... with 23,970 more rows, and 9 more variables: `Mg_mol/molCa` <dbl>, Mg_25 <dbl>,
##   `P_mol/molCa` <dbl>, `K_mol/molCa` <dbl>, `Mn_mol/molCa` <dbl>, Rb85 <dbl>,
##   `Sr_mol/molCa` <dbl>, `Ba_mol/molCa137` <dbl>, Ba138 <dbl>

# Just taking the continuous variables for our data
contdata <- fulldata[,-1:-5]

# Running the PCA (selected 7 axes as this seemed to be where there was a break)
pcacont <- dudi.pca(contdata,center=TRUE,scale=TRUE)

# Show the eigenvalues
fviz_eig(pcacont)

# Show the PCA plot (note, may take a while to display, colored by sample site)
fviz_pca_biplot(pcacont, habillage=fulldata$site,palette="LocusZoom")

# Plotting against other PCs
fviz_pca_biplot(pcacont, axes=c(3,4),habillage=fulldata$site,palette="LocusZoom")

# Plotting against other PCs
fviz_pca_biplot(pcacont, axes=c(5,6),habillage=fulldata$site,palette="LocusZoom")

#Extracting first three PCs to plot
plotdata <- tibble(site=fulldata$site,sample=fulldata$sample,x=pcacont$l1$RS1*pcacont$eig[1],y=pcacont$l1$RS2*pcacont$eig[2],z=pcacont$l1$RS3*pcacont$eig[3])

# Colors to chose from (tried to get as divergent as possible)
sitecolors <- c("Greys","Greens","Blues","Purples","PuRd","YlOrBr","BuGn")

# count for number of sites
no_of_sites <- length(unique(plotdata$site))
colind <- 4

# for loop to go through adding samples based on colour
for (i in 1:no_of_sites) {
  tempcolorpalette <- brewer.pal(n = 9,sitecolors[i])
  plotsitedata <- filter(plotdata,site==unique(plotdata$site)[i])
  no_of_individuals <- length(unique(plotsitedata$sample))
  for (j in 1:no_of_individuals) {
    plotsiteinddata <- filter(plotsitedata,sample==unique(plotsitedata$sample)[j])
    if (i ==1 & j ==1 ) {
      lines3D(x=plotsiteinddata$x,y=plotsiteinddata$y,z=plotsiteinddata$z,col=tempcolorpalette[colind],bty="bl2",add=FALSE)
    } else {
      lines3D(x=plotsiteinddata$x,y=plotsiteinddata$y,z=plotsiteinddata$z,col=tempcolorpalette[colind],bty="bl2",add=TRUE)
    }
    if (colind==6) {
      colind <- 4
    } else {
      colind <- 6
    }  
  }
}

print("The following sites:")
print(unique(plotdata$site))
print("are shown as the following colours")
print(sitecolors)

# Displays 3D plot
plotrgl()

# Doing some flat 2D time series stuff
twoplotdata <- tibble(distance=fulldata$`Distance (µm)`)
twoplotdata <- as_tibble(cbind(twoplotdata,plotdata))

# Plotting PC1 against distance
ggplot(twoplotdata,aes(x=distance,y=x,group=sample,color=site)) + geom_line() +
  xlab("Distance") +
  ylab("PC1")

#Faceting
ggplot(twoplotdata,aes(x=distance,y=x,group=sample,color=site)) + geom_line() +
  xlab("Distance") +
  ylab("PC1")+
  facet_grid(rows=vars(site))

 
# Plotting PC2 against distance
ggplot(twoplotdata,aes(x=distance,y=y,group=sample,color=site)) + geom_line() +
  xlab("Distance") +
  ylab("PC2")

#Faceting
ggplot(twoplotdata,aes(x=distance,y=y,group=sample,color=site)) + geom_line() +
  xlab("Distance") +
  ylab("PC2")+
   facet_grid(rows=vars(site))

# Plotting PC3 against distance
ggplot(twoplotdata,aes(x=distance,y=z,group=sample,color=site)) + geom_line() +
  xlab("Distance") +
  ylab("PC3")

#Faceting
ggplot(twoplotdata,aes(x=distance,y=z,group=sample,color=site)) + geom_line() +
  xlab("Distance") +
  ylab("PC3")+
  facet_grid(rows=vars(site))

# Plotting PC1 against PC2, coloring by time
ggplot(twoplotdata,aes(x=y,y=x,group=sample,color=distance)) + geom_line() +
  facet_grid(rows=vars(site))+
  xlab("PC2") +
  ylab("PC1")

# Plotting PC2 against PC3, coloring by time
ggplot(twoplotdata,aes(x=z,y=y,group=sample,color=distance)) + geom_line() +
  facet_grid(rows=vars(site))+
  xlab("PC3") +
  ylab("PC2")

# Looking at nontransformed data
ggplot(fulldata,aes(x=`Distance (µm)`,y=Ca43_CPS,group=sample,color=site)) + geom_line() +
  xlab("Distance (µm)") +
  ylab("Ca43_CPS")+
  facet_grid(rows=vars(site))

# Doing a ylim on this one because some individuals are blippy
ggplot(fulldata,aes(x=`Distance (µm)`,y=`Li_mol/molCa`,group=sample,color=site)) + geom_line() +
  xlab("Distance (µm)") +
  ylab("Li_mol/molCa")+
  facet_grid(rows=vars(site))+
  ylim(0,0.00003)

# Doing a ylim on this one because some individuals are blippy
ggplot(fulldata,aes(x=`Distance (µm)`,y=`Na`,group=sample,color=site)) + geom_line() +
  xlab("Distance (µm)") +
  ylab("Na")+
  facet_grid(rows=vars(site))+
  ylim(0,0.04)

# Doing a ylim on this one because some individuals are blippy
ggplot(fulldata,aes(x=`Distance (µm)`,y=`Mg_mol/molCa`,group=sample,color=site)) + geom_line() +
  xlab("Distance (µm)") +
  ylab("Mg_mol/molCa")+
  facet_grid(rows=vars(site))+
  ylim(0,0.001)

ggplot(fulldata,aes(x=`Distance (µm)`,y=`Mg_25`,group=sample,color=site)) + geom_line() +
  xlab("Distance (µm)") +
  ylab("Mg_25")+
  facet_grid(rows=vars(site))+
  ylim(0,0.00075)

#Creating new dataframe with just the mean for the first five, middle five, and last five records per individual
without_indcol <- fulldata %>% select(-indcol)

firstfive <- without_indcol %>% group_by(sample,site) %>% arrange(`Distance (µm)`) %>% slice(1:5) %>% summarise_all(mean)
firstfive <- add_column(firstfive, timepoint = "first", .before = 1)

middlefive <- without_indcol %>% group_by(sample,site) %>% arrange(`Distance (µm)`) %>% slice((n()/2-4):(n()/2)) %>% summarise_all(mean)
middlefive <- add_column(middlefive, timepoint = "middle", .before = 1)

lastfive <- without_indcol %>% group_by(sample,site) %>% arrange(`Distance (µm)`) %>% slice((n()-4):n()) %>% summarise_all(mean)
lastfive <- add_column(lastfive, timepoint = "last", .before = 1)

meanfive <- as_tibble(rbind(firstfive,middlefive,lastfive))

# Extracting just the continuous data from this
contfive <- meanfive[6:17]

# Reunning the PCA (selected 7 axes to be consistent with last time
pcacontfive <- dudi.pca(contfive,center=TRUE,scale=TRUE)

#useful name for dots
meanfive <- meanfive %>% mutate(comboname=paste(timepoint,"_",site,sep=""))
write.table(meanfive,"summary_first_middle_last_five.txt",row.names=FALSE,quote=FALSE)

# Grouping the data
grp <- find.clusters(contfive, max.n.clust=100)

# Using our actual groups
indgrp <- as.factor(meanfive$comboname)
names(indgrp) <- names(grp$grp)
indgrp <- as.factor(meanfive$comboname)

# doing a dapc on full contfive
dapc1 <- dapc(pcacontfive$l1,grp=indgrp)

scatter(dapc1)

fviz_pca_biplot(pcacontfive, habillage=meanfive$comboname,palette="LocusZoom")

# Doing it just on firstfive to look for spatial structure
contfirstfive <- firstfive[6:17]
pcacontfirstfive <- dudi.pca(contfirstfive,center=TRUE,scale=TRUE)
fviz_pca_biplot(pcacontfirstfive, habillage=firstfive$site,palette="LocusZoom")
grp <- find.clusters(contfirstfive, max.n.clust=34)
indgrp <- as.factor(firstfive$site)
names(indgrp) <- names(grp$grp)
indgrp <- as.factor(firstfive$site)
dapc1 <- dapc(pcacontfirstfive$l1,grp=indgrp)
scatter(dapc1)

# Doing it just on middlefive to look for spatial structure
contmiddlefive <- middlefive[6:17]
pcacontmiddlefive <- dudi.pca(contmiddlefive,center=TRUE,scale=TRUE)
fviz_pca_biplot(pcacontmiddlefive, habillage=middlefive$site,palette="LocusZoom")
grp <- find.clusters(contmiddlefive, max.n.clust=34)
indgrp <- as.factor(middlefive$site)
names(indgrp) <- names(grp$grp)
indgrp <- as.factor(middlefive$site)
dapc1 <- dapc(pcacontmiddlefive$l1,grp=indgrp)
scatter(dapc1)

# Doing it just on lastfive to look for spatial structure
contlastfive <- lastfive[6:17]
pcacontlastfive <- dudi.pca(contlastfive,center=TRUE,scale=TRUE)
fviz_pca_biplot(pcacontlastfive, habillage=lastfive$site,palette="LocusZoom")
grp <- find.clusters(contlastfive, max.n.clust=34)
indgrp <- as.factor(lastfive$site)
names(indgrp) <- names(grp$grp)
indgrp <- as.factor(lastfive$site)
dapc1 <- dapc(pcacontlastfive$l1,grp=indgrp)
scatter(dapc1)
