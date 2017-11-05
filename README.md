# plotting_functions
Some functions for using R to plot data prettily


# PCA plot of qPCR data, displaying stage/experiment as colours/symbols
```
PCA_plot_qPCR <- function(csv_file) {

if(missing(csv_file)) {
stop('Please give the full path and name of your csv file with qPCR counts e.g. PCA_plot_qPCR("/Users/alanaalexander/Downloads/Practice data.csv")')
}

library(ggplot2)
library(ggfortify)

input_data <- read.csv(csv_file)

colour_column <- NA

print("Please type the name of column that you would like to use colour to represent:")
print("You may choose out of the following columns.")
print(colnames(input_data))
colour_column <- readLines(con = stdin(), n=1)
if(is.null(input_data[ , colour_column ])) {
   print("Sorry, that does not seem to be a valid column name. Please check your spelling and try again")
   while(is.null(input_data[ , colour_column ])) {
      print("Please type the name of column that you would like to use colour to represent:")
      print("You may choose out of the following columns.")
      print(colnames(input_data))
      colour_column <- readLines(con = stdin(), n=1)
      if(is.null(input_data[ , colour_column ])) {
         print("Sorry, that does not seem to be a valid column name. Please check your spelling and try again")
      }    
   }
}

print("Please type the name of column that you would like to use colour to represent:")
print("You may choose out of the following columns.")
colour_column <- readLines(con = stdin(), n=1)







```
