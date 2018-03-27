#adegenet to generate a PCA
library(adegenet)
library(factoextra)
structure <- read.structure("lamprey.stru",n.ind=38,n.loc=21469,col.lab=1,col.pop=0,row.marknames=0)
replace_missing <- tab(structure, freq=TRUE,NA.method="mean")
pca1 <- dudi.pca(replace_missing,center=TRUE,scale=FALSE)
#select 2 axes
s.label(pca1$li)
eig.perc <- 100*pca1$eig/sum(pca1$eig)
add.scatter.eig(eig.perc,3,1,2, ratio=.3,posi="bottomright")
