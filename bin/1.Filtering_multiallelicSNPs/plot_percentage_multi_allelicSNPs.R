

#################################################################
#                                                              #    
#               Multi-allelic SNPs PLOT                        #
#                           R script                           #
#                                                              #
#################################################################

#-------------------------------------------------------------------------------

#This custom script has been developed by:

#       GI Sistemas Naturales e Historia Forestal (formerly known as GI Genetica, Fisiolog?a e Historia Forestal)
#       Dpto. Sistemas y Recursos Naturales
#       ETSI Montes, Forestal y del Medio Natural
#       Universidad Politécnica de Madrid
#       https://github.com/ggfhf/
#
#   Licence: GNU General Public Licence Version 3

#-------------------------------------------------------------------------------

#This script contains a custom workflow to plot the percentage of multi-allelic SNPs from vcf filtered files.
#Also this script plot the number of loci/SNPs post multi-allelic filtered
#-------------------------------------------------------------------------------
#Load required packages  
library(tidyverse)
library(easyGgplot2)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(githubinstall)
#####SET PATH TO YOUR WORKING DIRECTORY (where you have stored the directory 1.Filtering_multiallelicSNPs) 
dir<-("/PATH_TO_WORKING_DIRECTORY/bin/1.Filtering_multiallelicSNPs")
setwd(dir)
#####SET PATH TO YOUR INPUT FILES
dirInput<-("/PATH_TO_WORKING_DIRECTORY/Input_files") ##Same for all the scripts
#####LOAD A .csv FILE FORMAT WITH THE CHARACTERISTICS OF EACH COMBINATION (see Ch_combinations.csv example in the Input_files directory)
Ch_names <- read.csv(paste(dirInput, "Ch_combinations.csv", sep = "/"))
#####LOAD A .csv FILE FORMAT WITH YOUR SAMPLES LISTED (see Samples_names.txt example in the Input_files directory)
S_names <- read.table(paste(dirInput,"Samples_names.txt", sep = "/"))
#####LIST .vcf files INSIDE THE INPUT VCF FILES DIRECTORY (INPUT_VCF_DATA).
filelist <- list.files((paste(dirInput,"INPUT_VCF_DATA", sep = "/")), pattern="*.vcf", full.names=TRUE)

#####Create an object with the names of each combination merging the columns of the characteristics data frame
C_names<-paste(Ch_names$Read_type,Ch_names$Enzyme,Ch_names$Methodology,Ch_names$Read_number,Ch_names$parameters, sep = "_")

#####Create a transposed object of the samples names
S_names_t<- as.character(t(S_names))

############## Number of total SNPs ############ 

#####Read total .vcf files into a list file 
ldf_total <- lapply(filelist, read.table)

#####Define the name of each list element with the object defined as combination_names (C_names)
names(ldf_total)<-C_names

#####Define the columns names 
vcfColNames<-c("#CHROM", "POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT")

for (i in 1:length(ldf_total)){
  colnames(ldf_total[[i]])<-c(vcfColNames,S_names_t) #####First nine columns are shared for all .vcf (vcfColNames). Next columns are the sample names (S_names_t)
}

#####Quantify NtotalSNPs
N_totalSNPs<-c()
for (i in 1:length(ldf_total)){
  N_totalSNPs[i]<-nrow(ldf_total[[i]])
}

#####merge characteristics info with Ntotal_SNPS
N_totalSNPs_CH<-cbind(C_names,N_totalSNPs, Ch_names)

#####Filtered by number of reads and by methodology
RN<-unique(Ch_names$Read_number)
Meth<-unique(Ch_names$Methodology)

All_possible_comb<-expand.grid(RN,Meth)

colnames(All_possible_comb)<-c("Read_number", "Methodology")

N_totalSNPs_CH_filter<-list()
for (i in 1:nrow(All_possible_comb)) {
  N_totalSNPs_CH_filter[[i]]<-dplyr::filter(N_totalSNPs_CH, Read_number==as.character(All_possible_comb[i,1]),
                                     Methodology==as.character(All_possible_comb[i,2]))
}

############## Number of multi-allelic SNPs ############ 

#####List multi-allelic .vcf files inside the Multi_allelic directory
filelist_multiA <- list.files((paste(dir,"Multi_allelic", sep = "/")), pattern="*.vcf", full.names=TRUE)

#####Read multi-allelic .vcf files into a list file
ldf_multiA <- lapply(filelist_multiA, read.table)

#####Define the name of each list element with the object defined as combination_names (C_names)
names(ldf_multiA)<-C_names

#####Define the columns names 
vcfColNames<-c("#CHROM", "POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT")

for (i in 1:length(ldf_multiA)){
  colnames(ldf_multiA[[i]])<-c(vcfColNames,S_names_t) #####First nine columns are shared for all .vcf (vcfColNames). Next columns are the sample names (S_names_t)
}

#####Quantify N multi-allelic SNPs
N_multiASNPs<-c()
for (i in 1:length(ldf_multiA)){
  N_multiASNPs[i]<-nrow(ldf_multiA[[i]])
}

#####merge characteristics info with N multi-allelic SNPs
N_multiASNPs_CH<-cbind(C_names,N_multiASNPs, Ch_names)

#####Filtered by number of reads and by methodology
N_multiASNPs_CH_filter<-list()
for (i in 1:nrow(All_possible_comb)) {
  N_multiASNPs_CH_filter[[i]]<-dplyr::filter(N_multiASNPs_CH, Read_number==as.character(All_possible_comb[i,1]), 
                                     Methodology==as.character(All_possible_comb[i,2]))
}


############## Percentage of multi-allelic SNPs ############ 

perc_multiA<-c()
for (i in 1:length(ldf_total)){
  perc_multiA[i]<-((N_totalSNPs_CH$N_totalSNPs[i]-N_multiASNPs_CH$N_multiASNPs[i])-
                     N_totalSNPs_CH$N_totalSNPs[i])/N_totalSNPs_CH$N_totalSNPs[i]*100
}

#####merge characteristics info with Ntotal_SNPS
perc_multiA_CH<-cbind(C_names,perc_multiA, Ch_names)

#####Filtered by number of reads and by methodology
perc_multiA_CH_filter<-list()
for (i in 1:nrow(All_possible_comb)) {
  perc_multiA_CH_filter[[i]]<-dplyr::filter(perc_multiA_CH, Read_number==as.character(All_possible_comb[i,1]), 
                                      Methodology==as.character(All_possible_comb[i,2]))
}



#####PLOT % MULTI-ALLELIC SNP#######
plotperc_multiA<-list()
for (i in 1:nrow(All_possible_comb)) { 
  plotperc_multiA[[i]]<-ggplot(data = perc_multiA_CH_filter[[i]])+
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "darkgrey"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "darkgrey"))+
    geom_bar(mapping = aes(x=parameters, y=abs(perc_multiA), fill = Read_type), position = "dodge", stat = "identity")  +
    theme(axis.text.x = element_text( angle = 20, size = 9),axis.text.y  = element_text( angle = 0, size = 10))+
    ylab("% multi-allelic SNPs filtred") +xlab("parameters")+
    guides(fill=guide_legend("Read type"))+
    ylim(0,15)
  
}


####Create the directory to store the plots
dir.create(path = "PERCENTAGEmultiAllelic_PLOTS")

mypath<-list()
for (i in 1:length(plotperc_multiA)){
  
  mypath[[i]] <- file.path("PERCENTAGEmultiAllelic_PLOTS",paste("Percentage_multiAllelic_", 
                    paste(All_possible_comb[i,1],All_possible_comb[i,2], sep = "_"), ".jpg", sep = ""))
  
  jpeg(file=mypath[[i]], width = 25, height = 10, units = "cm", res = 400)
  print(plotperc_multiA[[i]])
  dev.off()
}

############## Number of bi-allelic SNPs ############ 

#####List multi-allelic .vcf files inside the Bi_allelic directory
filelist_biA <- list.files((paste(dir,"Bi_allelic", sep = "/")), pattern="*.vcf", full.names=TRUE)

#####Read bi-allelic .vcf files into a list file
ldf_biA <- lapply(filelist_biA, read.table)

#####Define the name of each list element with the object defined as combination_names (C_names)
names(ldf_biA)<-C_names

#####Define the columns names 
vcfColNames<-c("CHROM", "POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT")

for (i in 1:length(ldf_biA)){
  colnames(ldf_biA[[i]])<-c(vcfColNames,S_names_t) #####First nine columns are shared for all .vcf (vcfColNames). Next columns are the sample names (S_names_t)
}

#####Quantify N bi-allelic SNPs
N_biASNPs<-c()
for (i in 1:length(ldf_biA)){
  N_biASNPs[i]<-nrow(ldf_biA[[i]])
}

#####merge characteristics info with N bi-allelic SNPs
N_biASNPs_CH<-cbind(C_names,N_biASNPs, Ch_names)

#####Filtered by number of reads and by methodology
N_biASNPs_CH_filter<-list()
for (i in 1:nrow(All_possible_comb)) {
  N_biASNPs_CH_filter[[i]]<-dplyr::filter(N_biASNPs_CH, Read_number==as.character(All_possible_comb[i,1]), 
                                      Methodology==as.character(All_possible_comb[i,2]))
}

############## Number of loci ############ 
#Get only information from CHROM (locus) column
ldf_biA_loci<-ldf_biA
for (i in 1:length(ldf_biA)){
  ldf_biA_loci[[i]]<-subset(ldf_biA_loci[[i]], !duplicated(ldf_biA_loci[[i]]$CHROM)) 
}


#####Quantify N bi-allelic loci
N_biAloci<-c()
for (i in 1:length(ldf_biA_loci)){
  N_biAloci[i]<-nrow(ldf_biA_loci[[i]])
}

#####merge characteristics info with N bi-allelic SNPs
N_biAloci_CH<-cbind(C_names,N_biAloci, Ch_names)

#####Filtered by number of reads and by methodology
N_biAloci_CH_filter<-list()
for (i in 1:nrow(All_possible_comb)) {
  N_biAloci_CH_filter[[i]]<-dplyr::filter(N_biAloci_CH, Read_number==as.character(All_possible_comb[i,1]), 
                                   Methodology==as.character(All_possible_comb[i,2]))
}


#####Plot NUMBER loci/SNPS#######
plot_biA_SNPs<-list()
plot_biA_loci<-list()
for (i in 1:nrow(All_possible_comb)) { 
  plot_biA_loci[[i]]<-ggplot(data = N_biAloci_CH_filter[[i]])+
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "darkgrey"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "darkgrey"))+
    geom_bar(mapping = aes(x=parameters, y=N_biAloci, fill = Read_type), position = "dodge", stat = "identity")  +
    theme(axis.text.x = element_text( angle = 20, size = 9),axis.text.y  = element_text( angle = 0, size = 10))+
    ylab("Number of loci") +xlab("parameters")+
    guides(fill=guide_legend("Read type"))+
    ylim(0,50000)
  
  plot_biA_SNPs[[i]]<-ggplot(data = N_biASNPs_CH_filter[[i]])+
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "darkgrey"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "darkgrey"))+
    geom_bar(mapping = aes(x=parameters, y=N_biASNPs, fill = Read_type), position = "dodge", stat = "identity")  +
    theme(axis.text.x = element_text( angle = 20, size = 9),axis.text.y  = element_text( angle = 0, size = 10))+
    ylab("Number of SNPs") +xlab("parameters")+
    guides(fill=guide_legend("Read type"))+
    ylim(0,210000)
  
}


####Create the directory to store the plots
dir.create(path = "NlocibiA_PLOTS")
dir.create(path = "NSNPsbiA_PLOTS")

mypath<-list()
for (i in 1:length(plot_biA_loci)){
  
  mypath[[i]] <- file.path("NlocibiA_PLOTS",paste("Nloci_biAllelic_", 
                    paste(All_possible_comb[i,1],All_possible_comb[i,2], sep = "_"), ".jpg", sep = ""))
  
  jpeg(file=mypath[[i]], width = 25, height = 10, units = "cm", res = 400)
  print(plot_biA_loci[[i]])
  dev.off()
}


mypath<-list()
for (i in 1:length(plot_biA_SNPs)){
  
  mypath[[i]] <- file.path("NSNPsbiA_PLOTS",paste("NSNPs_biAllelic_",
                    paste(All_possible_comb[i,1],All_possible_comb[i,2], sep = "_"), ".jpg", sep = ""))
  
  jpeg(file=mypath[[i]], width = 25, height = 10, units = "cm", res = 400)
  print(plot_biA_SNPs[[i]])
  dev.off()
}


print("finished second R phase")
