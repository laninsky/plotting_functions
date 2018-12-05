#adegenet to generate a PCA
library(adegenet)
library(factoextra)
structure <- read.structure("lamprey.stru",n.ind=38,n.loc=21469,col.lab=1,col.pop=0,row.marknames=0)
replace_missing <- tab(structure, freq=TRUE,NA.method="mean")
pca1 <- dudi.pca(replace_missing,center=TRUE,scale=FALSE)
#select 2 axes

# Option 1 for plotting
fviz_eig(pca1)

fviz_pca_ind(pca1,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             )

fviz_pca_biplot(pca1, 
                col.ind = "#696969"  # Individuals color
                )

# Option 2 for plotting 
s.label(pca1$li,addaxes=TRUE)
eig.perc <- 100*pca1$eig/sum(pca1$eig)
add.scatter.eig(eig.perc,3,1,2, ratio=.3,posi="bottomright")
#add.scatter.eig(eig.perc,3,1,2, ratio=.3,posi="topright")

#DAPC
grp <- find.clusters(structure, max.n.clust = 10)
# This will display a graph of the variance explained by PCA and you can chose the number of PCs to retain
# select 1
# select 2 clusters
dapc1 <- dapc(structure, grp$grp)
# This will display a graph of the variance explained by PCA and you can chose the number of PCs to retain
#select 1
