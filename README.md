# plotting_functions
Some functions for using R to plot data prettily


#PCA plot, displaying stage/overlapping experiment as colours/symbols
```
# This function expects a pathway and name of a csv file with qPCR count data in it. In the same folder, appended by ".cols", it expects
#you on the first line of this plain text file to give the name of the column you would like to display by colour in the PCA plot,
#followed by the names of the columns with the count data.

PCA_plot_qPCR <- function(csv_file) {

if(missing(csv_file)) {
stop('Please give the full path and name of your csv file with qPCR counts e.g. PCA_plot_qPCR("/Users/alanaalexander/Downloads/Practice data.csv")')
}

#loading required libraries
library(ggplot2)
library(ggfortify)
library(GGally)
library(pcaMethods)

#reading in data and column information
input_data <- read.csv(csv_file)
input_cols <- as.matrix(read.table(paste(csv_file,".cols",sep=""),sep="\t"))
input_cols <- gsub("\\s", ".", input_cols)

#extracting columns with numerical data and centering the values
pca_data <- input_data[ , input_cols[2:(dim(input_cols)[1]),1] ]

#For loop for looking at different scaling methods
for (scale_type in c("none","pareto","vector","uv")) {
pca_data_C <- prep(pca_data, scale=scale_type, center=TRUE)

#If control structure for imputing missing data if necessary
if(!(dim(pca_data_C)[1]==dim(na.omit(pca_data_C))[1])) {

#Getting error estimates for each imputation method - original PCA ("svd") is calculated by dropping all rows with NA
errPCA  <- kEstimate(na.omit(pca_data_C),method="svd",evalPcs=1:((dim(pca_data_C)[2])-1),allVariables=TRUE)
errPPCA <- kEstimate(pca_data_C,method="ppca",evalPcs=1:((dim(pca_data_C)[2])-1),allVariables=TRUE)
errBPCA <- kEstimate(pca_data_C,method="bpca",evalPcs=1:((dim(pca_data_C)[2])-1),allVariables=TRUE)
errSVDI <- kEstimate(pca_data_C,method="svdImpute",evalPcs=1:((dim(pca_data_C)[2])-1),allVariables=TRUE)
errNipals <- kEstimate(pca_data_C,method="nipals",evalPcs=1:((dim(pca_data_C)[2])-1),allVariables=TRUE)
errNLPCA <- kEstimate(pca_data_C,method="nlpca",evalPcs=1:((dim(pca_data_C)[2])-1),maxSteps=300,allVariables=TRUE)

#creating a matrix with the error estiamtes for each imputation method
header_row <- c("methods","PCA","PPCA","BPCA","SVDI","Nipals","NLPCA")
row_names <- c("bestNPcs","eError_1","eError_2","eErr_SS")
error_matrix <- matrix(NA,nrow=5,ncol=7)
error_matrix[1,] <- header_row
error_matrix[2:(dim(error_matrix)[1]),1] <- row_names
error_matrix[2:4,2] <- c(errPCA$bestNPcs,errPCA$eError[1],errPCA$eError[2])
error_matrix[2:4,3] <- c(errPPCA$bestNPcs,errPPCA$eError[1],errPPCA$eError[2])
error_matrix[2:4,4] <- c(errBPCA$bestNPcs,errBPCA$eError[1],errBPCA$eError[2])
error_matrix[2:4,5] <- c(errSVDI$bestNPcs,errSVDI$eError[1],errSVDI$eError[2])
error_matrix[2:4,6] <- c(errNipals$bestNPcs,errNipals$eError[1],errNipals$eError[2])
error_matrix[2:4,7] <- c(errNLPCA$bestNPcs,errNLPCA$eError[1],errNLPCA$eError[2])

for (i in 2:7) {
   error_matrix[5,i] <- as.numeric(error_matrix[3,i])^2+as.numeric(error_matrix[4,i])^2
}

#Writing out the matrix - you can use this to choose the imputation method that you prefer
dir_create_name <- paste(getwd(),"/",scale_type,sep="")
write.table(error_matrix,paste(dir_create_name,"/imputation_method_error_matrix.txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Loop for imputing missing data values for the differing numbers of principle components suggested by the kestimate step
for (npcs in min(as.numeric(na.omit(error_matrix[2,2:7]))):max(as.numeric(na.omit(error_matrix[2,2:7])))) {

#creating a directory for each optimal number of principle components
dir_create_name <- paste(getwd(),"/",scale_type,"/Imputed_Data_NPcs_",npcs,sep="")
dir.create(dir_create_name)

#imputing missing data
resPCA <- completeObs(pca(na.omit(pca_data_C), method="svd", center=FALSE, nPcs=npcs))
resPPCA <- completeObs(pca(pca_data_C, method="ppca", center=FALSE, nPcs=npcs))
resBPCA <- completeObs(pca(pca_data_C, method="bpca", center=FALSE, nPcs=npcs))
resSVDI <- completeObs(pca(pca_data_C, method="svdImpute", center=FALSE, nPcs=npcs))
resNipals <- completeObs(pca(pca_data_C, method="nipals", center=FALSE, nPcs=npcs))
resNLPCA <- completeObs(pca(pca_data_C, method="nlpca", center=FALSE, nPcs=npcs, maxSteps=300))

#taking note of missing rows for plotting
missing_rows <- attr((na.omit(pca_data_C)),"na.action")
missing_rows <- missing_rows[1:length(missing_rows)]

#Running the PCA on each of the imputed datasets
pcaPCA <- prcomp(resPCA)
pcaPPCA <- prcomp(resPPCA)
pcaBPCA <- prcomp(resBPCA)
pcaSVDI <- prcomp(resSVDI)
pcaNipals <- prcomp(resNipals)
pcaNLPCA <- prcomp(resNLPCA)

#Looking at the sdev of principle components for each imputational method and plotting this
eigenvalue_df <- data.frame(Eigenvalues=rep(c(1:length(pcaPCA$sdev)),6),Method=c(rep("PCA",(length(pcaPCA$sdev))),rep("PPCA",(length(pcaPCA$sdev))),rep("BPCA",(length(pcaPCA$sdev))),rep("SVDI",(length(pcaPCA$sdev))),rep("Nipals",(length(pcaPCA$sdev))),rep("NLPCA",(length(pcaPCA$sdev)))),sdev=c(pcaPCA$sdev,pcaPPCA$sdev,pcaBPCA$sdev,pcaSVDI$sdev,pcaNipals$sdev,pcaNLPCA$sdev))

ggplot(data=eigenvalue_df,aes(x=Eigenvalues,y=sdev,group=Method))+geom_line(aes(color=Method))+geom_point(aes(color=Method))+ggplot2::labs(title="Eigenvalue structure as obtained with different imputation methods",x="Eigenvalue", y="Standard deviation of PC")

ggsave("eigenvalue_structure.pdf", plot = last_plot(), device = NULL, path = dir_create_name)

#Looking at the loadings of each imputation method and plotting them against PCA
loadings_df <- data.frame(Method=c(rep("PCA",3),rep("PPCA",3),rep("BPCA",3),rep("SVDI",3),rep("Nipals",3),rep("NLPCA",3)),PC1=rbind(pcaPCA$rotation,pcaPPCA$rotation,pcaBPCA$rotation,pcaSVDI$rotation,pcaNipals$rotation,pcaNLPCA$rotation)[,1],PC1=rbind(pcaPCA$rotation,pcaPPCA$rotation,pcaBPCA$rotation,pcaSVDI$rotation,pcaNipals$rotation,pcaNLPCA$rotation)[,1],PC1=rbind(pcaPCA$rotation,pcaPPCA$rotation,pcaBPCA$rotation,pcaSVDI$rotation,pcaNipals$rotation,pcaNLPCA$rotation)[,3])



}

} else {
print("your dataset is complete, so no imputation has been performed")
errPCA  <- kEstimate(na.omit(pca_data_C),method="svd",evalPcs=1:((dim(pca_data_C)[2])-1),allVariables=TRUE)


}

#na.omit (regular, not imputed PCA, dropping rows that do not have values in each column
autoplot(prcomp(na.omit(pca_data_C)))

autoplot(prcomp(na.omit(pca_data_C)), data = na.omit(input_data), colour = input_cols[1,1], loadings = TRUE, loadings.label = TRUE)

}

```
