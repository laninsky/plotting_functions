# plotting_functions
Some functions for using R to plot data prettily


#PCA plot of qPCR data, displaying stage/overlapping experiment as colours/symbols
```
# This function expects a pathway and name of a csv file with qPCR count data in it. In the same folder, appended by ".cols", it expects
#you on the first line of this plain text file to give the name of the column you would like to display by colour in the PCA plot,
#followed by the names of the columns with the count data

PCA_plot_qPCR <- function(csv_file) {

if(missing(csv_file)) {
stop('Please give the full path and name of your csv file with qPCR counts e.g. PCA_plot_qPCR("/Users/alanaalexander/Downloads/Practice data.csv")')
}

library(ggplot2)
library(ggfortify)
library(pcaMethods)

input_data <- read.csv(csv_file)
input_cols <- as.matrix(read.table(paste(csv_file,".cols",sep=""),sep="\t"))
input_cols <- gsub("\\s", ".", input_cols)

pca_data <- input_data[ , input_cols[2:(dim(input_cols)[1]),1] ]

q2PPCA <- Q2((pca(pca_data, method="ppca", center=FALSE, nPcs=2)), pca_data, fold=2)



resPCA <- completeObs(pca(pca_data, method="svd", center=FALSE, nPcs=2))
resPPCA <- completeObs(pca(pca_data, method="ppca", center=FALSE, nPcs=2))
resBPCA <- completeObs(pca(pca_data, method="bpca", center=FALSE, nPcs=2))
resSVDI <- completeObs(pca(pca_data, method="svdImpute", center=FALSE, nPcs=2))
resNipals <- completeObs(pca(pca_data, method="nipals", center=FALSE, nPcs=2))
resNLPCA <- completeObs(pca(pca_data, method="nlpca", center=FALSE, nPcs=2, maxSteps=300))


#na.omit (regular, not imputed PCA, dropping rows that do not have values in each column
autoplot(prcomp(na.omit(input_data[ , input_cols[2:(dim(input_cols)[1]),1] ])))

autoplot(prcomp(na.omit(input_data[ , input_cols[2:(dim(input_cols)[1]),1] ])), data = na.omit(input_data), colour = input_cols[1,1], loadings = TRUE, loadings.label = TRUE)

}

```
