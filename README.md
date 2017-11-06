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

#loading required libraries
library(ggplot2)
library(ggfortify)
library(pcaMethods)

#reading in data and column information
input_data <- read.csv(csv_file)
input_cols <- as.matrix(read.table(paste(csv_file,".cols",sep=""),sep="\t"))
input_cols <- gsub("\\s", ".", input_cols)

#extracting columns with numerical data and centering the values
pca_data <- input_data[ , input_cols[2:(dim(input_cols)[1]),1] ]
pca_data_C <- prep(pca_data, scale="none", center=TRUE)

#imputing missing data if necessary
if(!(dim(pca_data_C)[1]==dim(na.omit(pca_data_C))[1])) {

errPCA  <- kEstimate(na.omit(pca_data_C),method="svd",evalPcs=1:((dim(pca_data_C)[2])-1),allVariables=TRUE)
errPPCA <- kEstimate(pca_data_C,method="ppca",evalPcs=1:((dim(pca_data_C)[2])-1))
errBPCA <- kEstimate(pca_data_C,method="bpca",evalPcs=1:((dim(pca_data_C)[2])-1))
errSVDI <- kEstimate(pca_data_C,method="svdImpute",evalPcs=1:((dim(pca_data_C)[2])-1))
errNipals <- kEstimate(pca_data_C,method="nipals",evalPcs=1:((dim(pca_data_C)[2])-1))
errNLPCA <- kEstimate(pca_data_C,method="nlpca",evalPcs=1:((dim(pca_data_C)[2])-1))

header_row <- c("methods","PCA","PPCA","BPCA","SVDI","Nipals","NLPCA")
row_names <- c("bestNPcs","eError_1","eError_2")
error_matrix <- matrix(NA,nrow=4,ncol=7)
error_matrix[1,] <- header_row
error_matrix[2:(dim(error_matrix)[1]),1] <- row_names
error_matrix[2:(dim(error_matrix)[1]),2] <- c(errPCA$bestNPcs,errPCA$eError[1],errPCA$eError[2])
error_matrix[2:(dim(error_matrix)[1]),3] <- c(errPPCA$bestNPcs,errPPCA$eError[1],errPPCA$eError[2])
error_matrix[2:(dim(error_matrix)[1]),4] <- c(errBPCA$bestNPcs,errBPCA$eError[1],errBPCA$eError[2])

error_matrix[2:(dim(error_matrix)[1]),5] <- c(errPCA$bestNPcs,errPCA$eError[1],errPCA$eError[2])
error_matrix[2:(dim(error_matrix)[1]),6] <- c(errPCA$bestNPcs,errPCA$eError[1],errPCA$eError[2])
error_matrix[2:(dim(error_matrix)[1]),7] <- c(errPCA$bestNPcs,errPCA$eError[1],errPCA$eError[2])




# imputing missing data values
resPCA <- completeObs(pca(na.omit(pca_data_C), method="svd", center=FALSE, nPcs=(dim(pca_data_C)[2])))
resPPCA <- completeObs(pca(pca_data_C, method="ppca", center=FALSE, nPcs=(dim(pca_data_C)[2])))
resBPCA <- completeObs(pca(pca_data_C, method="bpca", center=FALSE, nPcs=(dim(pca_data_C)[2])))
resSVDI <- completeObs(pca(pca_data_C, method="svdImpute", center=FALSE, nPcs=(dim(pca_data_C)[2])))
resNipals <- completeObs(pca(pca_data_C, method="nipals", center=FALSE, nPcs=(dim(pca_data_C)[2])))
resNLPCA <- completeObs(pca(pca_data_C, method="nlpca", center=FALSE, nPcs=(dim(pca_data_C)[2])), maxSteps=300)

} else {
print("your dataset is complete, so no imputation has been performed")
errPCA  <- kEstimate(na.omit(pca_data_C),method="svd",evalPcs=1:((dim(pca_data_C)[2])-1),allVariables=TRUE)
}

#na.omit (regular, not imputed PCA, dropping rows that do not have values in each column
autoplot(prcomp(na.omit(pca_data_C)))

autoplot(prcomp(na.omit(pca_data_C)), data = na.omit(input_data), colour = input_cols[1,1], loadings = TRUE, loadings.label = TRUE)

}

```
