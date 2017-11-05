# plotting_functions
Some functions for using R to plot data prettily


# PCA plot of qPCR data, displaying stage/experiment as colours/symbols
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

input_data <- read.csv(csv_file)
input_cols <- as.matrix(read.table(paste(csv_file,".cols",sep=""),sep="\t"))
input_cols <- gsub("\\s", ".", input_cols)

autoplot(prcomp(input_data[ , input_cols[2:(dim(input_cols)[1]),1] ]))

autoplot(prcomp(input_data[ , input_cols[2:(dim(input_cols)[1]),1] ]), data = input_data, colour = input_cols[1,1], loadings = TRUE, loadings.label = TRUE)

}

```
