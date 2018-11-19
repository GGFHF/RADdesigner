#####Install and load required packages
#These packages were tested in the R3.5.1. version
##STAGE 1###
install.packages("tidyverse")
install.packages("devtools",dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)
install.packages("gridExtra")
install.packages("githubinstall")
library(devtools)         #We need to load devtools before using install_github()
install_github("kassambara/easyGgplot2")

##STAGE 2###
install.packages("vcfR", dependencies=TRUE)
install.packages("ggmap", dependencies=TRUE)
install.packages("StAMPP", dependencies=TRUE)
install.packages("gplots", dependencies=TRUE)
install.packages("adegraphics", dependencies=TRUE)
install.packages("ggmap", dependencies=TRUE)
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
install_github("sdray/ade4")
install_github("thibautjombart/adegenet")

##STAGE 3###
#No packages needed

##STAGE 4###
install_github("apache/spark/R/pkg")
install.packages('sqldf')
install.packages('gridExtra')
install.packages('grid')
install.packages('rlist')


