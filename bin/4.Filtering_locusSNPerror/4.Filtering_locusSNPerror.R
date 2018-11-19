

#################################################################
#                                                              #    
#           Filter Locus-error-rate & SNP-error-rate           #
#                           R script                           #
#                                                              #
#################################################################

#-------------------------------------------------------------------------------

#This custom script has been developed by:

#       GI Sistemas Naturales e Historia Forestal (formerly known as GI Genetica, Fisiologia e Historia Forestal)
#       Dpto. Sistemas y Recursos Naturales
#       ETSI Montes, Forestal y del Medio Natural
#       Universidad Politécnica de Madrid
#       https://github.com/ggfhf/
#
#   Licence: GNU General Public Licence Version 3


#-------------------------------------------------------------------------------

#This script contains a custom workflow to filter those loci/SNPs that did not match between replicates to eliminate locus/SNP errors. 
#-------------------------------------------------------------------------------
#Load required packages
library(tidyverse)
library(SparkR)
library(sqldf)
library(rlist)


#####SET PATH TO YOUR WORKING DIRECTORY (where you have stored the directory 4.Filtering_locusSNPerror) 
dir<-("/PATH_TO_WORKING_DIRECTORY/bin/4.Filtering_locusSNPerror")
setwd(dir)
#####SET PATH TO YOUR INPUT FILES
dirInput<-("/PATH_TO_WORKING_DIRECTORY/Input_files") ##Same for all the scripts
#####LOAD A .csv FILE FORMAT WITH THE CHARACTERISTICS OF EACH COMBINATION (see Ch_combinations.csv example in the Input_files directory)
Ch_names <- read.csv(paste(dirInput, "Ch_combinations.csv", sep = "/"))
#####LOAD A .txt FILE FORMAT WITH YOUR SAMPLES LISTED (see Samples_names.txt example in the Input_files directory)
S_names <- read.table(paste(dirInput,"Samples_names.txt", sep = "/"))
#####SET PATH TO vcf converted and prepared txt files (which are in the 3.Prepare_vcf directory)
dirInput_txt<-("/PATH_TO_WORKING_DIRECTORY/bin/3.Prepare_vcf/OUTPUT_txt")

#####Create an object with the names of each combination merging the columns of the characteristics data frame
C_names<-paste(Ch_names$Read_type,Ch_names$Enzyme,Ch_names$Methodology,Ch_names$Read_number,Ch_names$parameters, sep = "_")

#####Create a transposed object of the samples names
S_names_t<- as.character(t(S_names))

#####List .txt files inside the input directory (INPUT_sqldf)
filelist <- list.files(dirInput_txt, pattern="*.txt", full.names=TRUE)

#####Read the .CSV into a list file 
ldf_Chrom <- lapply(filelist, read.table)

#####Define the name of each list element with the object defined as combination_names (C_names)
names(ldf_Chrom)<-C_names

#####Link CHROM and POS columns to keep SNPs and defined the columns names
vcfContrast<-c("CHROM","POS")

ldf<-ldf_Chrom
for (i in 1:length(ldf_Chrom)) {
  colnames(ldf_Chrom[[i]])<-c(vcfContrast, S_names_t)
  colnames(ldf[[i]])<-c(vcfContrast,S_names_t)
  ldf[[i]]<-unite(ldf[[i]],CHROM_POS, CHROM, POS)
  }

#####Delete the Chrom_pos column to filtered the Locus/SNP error rates
ldf_compare<-ldf
for (i in 1:length(ldf_compare)) {
  ldf_compare[[i]]<-ldf_compare[[i]][,-1]
}

#####Compare each pair of individual+duplicate and create data frames of TRUE-FALSE
TF_record<-list()
for (i in 1: length(ldf_compare)){
  TF_record[[i]] <- ldf_compare[[i]][,seq(1,ncol(ldf_compare[[i]]),2)]==ldf_compare[[i]][,seq(2,ncol(ldf_compare[[i]]),2)]
}

#####Define the name of each list element with the object defined as combination_names (C_names)
names(TF_record)<-C_names

#####Quantify overall concordance for all samples
concordance<-list()
for (i in 1: length(ldf_compare)){
  concordance[[i]]<-colSums(ldf_compare[[i]][,seq(1,ncol(ldf_compare[[i]]),2)]==ldf_compare[[i]][,seq(2,ncol(ldf_compare[[i]]),2)])/nrow(ldf_compare[[i]])
}

#####Define the name of each list element with the object defined as combination_names (C_names)
names(concordance)<-C_names 

######Filter locus/SNP error 
ldf_filtered_SNP_error_total <- vector("list", length(ldf_compare))
for (i in 1:length(ldf)){
for (j in 1:nrow(TF_record[[i]])){
  if (all(TF_record[[i]][j,])==TRUE)
  {
    ldf_filtered_SNP_error_total[[i]] <- rbind(ldf_filtered_SNP_error_total[[i]][,], ldf[[i]][j,])
  }
  else {}
  }
}

#####Separate Chrom_pos column

ldf_filtered_Locus_error_total<-ldf_filtered_SNP_error_total
for (i in 1:length(ldf_filtered_SNP_error_total)){
  ldf_filtered_Locus_error_total[[i]]<-separate(ldf_filtered_Locus_error_total[[i]], CHROM_POS, sep = "_", remove = T, into = c("CHROM_POS_1", "CHROM_POS_2", "CHROM_POS_3"))
  }

##Get only information from uniques CHROM (locus)
for (i in 1:length(ldf_filtered_Locus_error_total)){
  ldf_filtered_Locus_error_total[[i]]<-subset(ldf_filtered_Locus_error_total[[i]], !duplicated(ldf_filtered_Locus_error_total[[i]]$CHROM_POS_2)) 
  ldf_filtered_Locus_error_total[[i]]<-unite(ldf_filtered_Locus_error_total[[i]],CHROM_POS, CHROM_POS_1, CHROM_POS_2, CHROM_POS_3)
  }

####Define the name of each list element with the object defined as combination_names (C_names)
names(ldf_filtered_SNP_error_total)<-C_names
names(ldf_filtered_Locus_error_total)<-C_names

####Final number of SNPs after the filtering
final_NUM_SNPs<-list()
for (i in 1:length(ldf_filtered_SNP_error_total)){
  final_NUM_SNPs[[i]]<-nrow(ldf_filtered_SNP_error_total[[i]])
}

####Final number of Loci after the filtering
final_NUM_Loci<-list()
for (i in 1:length(ldf_filtered_Locus_error_total)){
  final_NUM_Loci[[i]]<-nrow(ldf_filtered_Locus_error_total[[i]])
}

####Define the name of each list element with the object defined as combination_names (C_names)
names(final_NUM_SNPs) <- C_names
names(final_NUM_Loci) <- C_names

#Create directory to store the outputs
dir.create(path = "OUTPUT_FINAL_FILTERED/")
dir.out.final.filt<-paste(dir,"/OUTPUT_FINAL_FILTERED/", sep = "")
  
####Write out the final number of Loci/SNPs
final_NUM_SNPSs_Loci<-list()
for (i in 1: length(final_NUM_SNPs)){
final_NUM_SNPSs_Loci[[i]]<-do.call(rbind, Map(data.frame, A=final_NUM_SNPs[[i]], B=final_NUM_Loci[[i]]))
colnames(final_NUM_SNPSs_Loci[[i]])<-c("N.SNPs","N.Loci")
write.csv( x = final_NUM_SNPSs_Loci[[i]], file = paste(dir.out.final.filt,C_names[i],"_final_number_SNPS_Loci.csv"), row.names = F)
}

#Create directory to store the final filtered vcf
dir.create(path = paste(dir.out.final.filt, "FINAL_FILTERED_VCF_SNPs/", sep = ""))
dir.vcf.filt.SNPs<-paste(dir.out.final.filt,"FINAL_FILTERED_VCF_SNPs/", sep = "")

dir.create(path = paste( dir.out.final.filt, "FINAL_FILTERED_VCF_Loci/", sep = ""))
dir.vcf.filt.Loci<-paste(dir.out.final.filt,"FINAL_FILTERED_VCF_Loci/", sep = "")

####Write out the final filtered vcf
for (i in 1: length(ldf_filtered_SNP_error_total)){
write.table( x = ldf_filtered_SNP_error_total[[i]], file = paste(dir.vcf.filt.SNPs,C_names[i],"_final_SNPs_filtered.vcf", sep = ""),quote = F, row.names = F )
write.table( x = ldf_filtered_Locus_error_total[[i]], file = paste(dir.vcf.filt.Loci,C_names[i],"_final_Loci_filtered.vcf", sep = ""),quote = F, row.names = F )
  
  }
