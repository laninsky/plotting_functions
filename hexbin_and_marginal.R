# Defining the directory/file we are using
working_dir <- "/Users/alanaalexander/Downloads"
filename <- "Adephaga_11.1Kv1_probe_stats.txt"

# Loading in libraries we need  
library(readr)
library(ggplot2)
library(hexbin)
library(ggExtra)
library(grid)

# Setting the variables to our bits and pieces  
setwd("/Users/alanaalexander/Downloads")  
temp <- read_table2(filename)

# Making a transparent theme (so that we don't have a whole bunch of stuff plotted for the geom_points() option
transparent_theme <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_blank(),
  panel.background = element_rect(fill = "transparent",colour = NA),
  plot.background = element_rect(fill = "transparent",colour = NA))

# Making our hexagon binned data plot
GC_Tm_hex <- ggplot(data=temp,mapping=aes(x=`GC%`,y=Tm)) + 
  geom_hex(bins=10) + 
  theme_bw() +
  theme(aspect.ratio=1) +
  theme(legend.position="bottom")

# Pulling out some of the characteristics of our data. B/c of our very shonky overlaying of the hexbin on top of the other graph
# Needed to draw a line off the bottom left of the plot to give us enough space for our hexbin y-axis and legend to show up by
ymin = min(temp$Tm)
ymax = max(temp$Tm)
yrange = ymax-ymin
xmin = min(temp$`GC%`)
xmax = max(temp$`GC%`)
xrange = xmax-xmin

# Controls how much space there is on the left and bottom of the plot
multiplier <- 0.25

# The "dummy" geom _points layer (with the invisible line to make more space off the left and bottom of the plot)
GC_Tm_points <- ggplot(data=temp,mapping=aes(x=`GC%`,y=Tm)) + 
  geom_point(color="transparent") +
  transparent_theme +
  annotate("segment", x=(xmin-multiplier*xrange),xend=xmax,y=(ymin-multiplier*yrange),yend=ymax,colour="transparent") +
  theme(aspect.ratio=1)

# The boxplot marginal plot based off the geom_points plot
baseplot <- ggMarginal(GC_Tm_points  +
                         theme(aspect.ratio=1) +
                         theme(
                           panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(),
                           panel.background = element_rect(fill = "transparent",colour = NA),
                           plot.background = element_rect(fill = "transparent",colour = NA)
                         ) +
                         theme(legend.position="left"),
                         type="boxplot") 
)

# Printing the boxplots
print(baseplot,bg="transparent")

# Defining where we are going to shove the hexagon graph over top of the boxplots
vp <- viewport(x=0.453,y=0.41,width=0.77,height=0.77)

# If you want more space/bigger graph you could use these options
#multiplier <- 1/3
#vp <- viewport(x=0.465,y=0.42,width=0.74,height=0.74)

# Chucking the hexagon graph over top using the definitions we made above for vp
print(GC_Tm_hex,vp = vp)
