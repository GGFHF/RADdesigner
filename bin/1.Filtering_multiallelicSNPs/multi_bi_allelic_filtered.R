

#################################################################
#                                                              #    
#               Multi-allelic SNPs filtering                   #
#                           R script                           #
#                                                              #
#################################################################
                     
#-------------------------------------------------------------------------------

#This custom script has been developed by:

#       GI Sistemas Naturales e Historia Forestal (formerly known as GI Genetica, Fisiolog?a e Historia Forestal)
#       Dpto. Sistemas y Recursos Naturales
#       ETSI Montes, Forestal y del Medio Natural
#       Universidad Politecnica de Madrid
#       https://github.com/ggfhf/
#
#   Licence: GNU General Public Licence Version 3

#-------------------------------------------------------------------------------

#This script contains a custom workflow to filter multi-allelic SNPs from vcf files.
#-------------------------------------------------------------------------------
#Load required packages
library(tidyverse)
library(easyGgplot2)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(githubinstall)
#-------------------------------------------------------------------------------                           
#####SET PATH TO YOUR WORKING DIRECTORY (where you have stored the directory 1.Filtering_multiallelicSNPs) 
dir<-("/PATH_TO_WORKING_DIRECTORY/bin/1.Filtering_multiallelicSNPs")
setwd(dir)
#####SET PATH TO YOUR INPUT FILES
dirInput<-("/PATH_TO_WORKING_DIRECTORY/Input_files") ##Same for all the scripts
#####LOAD A .csv FILE FORMAT WITH THE CHARACTERISTICS OF EACH COMBINATION (see Ch_combinations.csv example in the Input_files directory)
Ch_names <- read.csv(paste(dirInput, "Ch_combinations.csv", sep = "/"))
#####LOAD A .csv FILE FORMAT WITH YOUR SAMPLES LISTED (see Samples_names.txt example in the Input_files directory)
S_names <- read.table(paste(dirInput,"Samples_names.txt", sep = "/"))
#####List .vcf files inside the input vcf files directory (INPUT_VCF_DATA). 
filelist <- list.files((paste(dirInput,"INPUT_VCF_DATA", sep = "/")), pattern="*.vcf", full.names=TRUE)

#####Create an object with the names of each combination merging the columns of the characteristics data frame
C_names<-paste(Ch_names$Read_type,Ch_names$Enzyme,Ch_names$Methodology,Ch_names$Read_number,Ch_names$parameters, sep = "_")

#####Create a transposed object of the samples names
S_names_t<- as.character(t(S_names))

#####Read the .vcf into a list file 
ldf <- lapply(filelist, read.table)

#####Define the name of each list element with the object defined as combination_names (C_names)
names(ldf)<-C_names

#####Defined the columns names 
vcfColNames<-c("#CHROM", "POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT")

for (i in 1:length(ldf)){
  colnames(ldf[[i]])<-c(vcfColNames,S_names_t) #####First nine columns are shared for all .vcf (vcfColNames). Next columns are the sample names (S_names_t)
}

#####Separation of .vcf files into multi-allelic SNPs and bi-allelic SNPs   

for (i in 1:length(ldf)){
nam<-paste(paste(names(ldf)[i], sep = ""), "multi_allelic", sep = "_")
 assign(nam,ldf[[i]] %>% filter(str_detect(ALT,",")))
 write.table(assign(nam,ldf[[i]] %>% filter(str_detect(ALT,","))), paste(nam,"vcf",sep = "."), 
              row.names = FALSE, sep = "	", quote = FALSE)
}
  
for (i in 1:length(ldf)){
  nam<-paste(paste(names(ldf)[i], sep = ""), "bi_allelic", sep = "_")
  assign(nam,ldf[[i]] %>% filter(!str_detect(ALT,",")))
  write.table(assign(nam,ldf[[i]] %>% filter(!str_detect(ALT,","))), paste(nam,"vcf",sep = "."), 
              row.names = FALSE, sep = "	", quote = FALSE)
}


print("finished first R phase")
