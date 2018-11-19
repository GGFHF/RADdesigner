

#################################################################
#                                                              #    
#      Quantify depth, Locus-error-rate & SNP-error-rate       #
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

#This script contains a custom workflow for quantify and plot the read-depth and both locus/SNP error rates. 
#Also this script quantify and plot the dendrograms of similarity for each sample and for each combination
#-------------------------------------------------------------------------------
#Load required packages  
library(vcfR)
library(StAMPP)
library(reshape2)
library(gplots)
library(adegraphics)
library(ggmap)
library(adegenet)
library(pegas)
library(phyloseq)
library(ade4)
#####SET PATH TO YOUR WORKING DIRECTORY (where you have stored the directory 2.Depth_locusSNPerror) 
dir<-("/PATH_TO_WORKING_DIRECTORY/bin/2.Depth_locusSNPerror")
setwd(dir)
#####SET PATH TO YOUR INPUT FILES
dirInput<-("/PATH_TO_WORKING_DIRECTORY/Input_files") ##Same for all the scripts
#####SET PATH TO vcf filtered files (which are in the 1.Filtering_multiallelicSNPs directory )
dirInput_vcf<-("/PATH_TO_WORKING_DIRECTORY/bin/1.Filtering_multiallelicSNPs/FILTERED_OUTPUT")

#####LOAD A .csv FILE FORMAT WITH THE CHARACTERISTICS OF EACH COMBINATION (see Ch_combinations.csv example in the Input_files directory)
Ch_names <- read.csv(paste(dirInput, "Ch_combinations.csv", sep = "/"))
#####LOAD A .txt FILE FORMAT WITH YOUR SAMPLES LISTED (see Samples_names.txt example in the Input_files directory)
S_names <- read.table(paste(dirInput,"Samples_names.txt", sep = "/"))

#####Create an object with the names of each combination merging the columns of the characteristics data frame
C_names<-paste(Ch_names$Read_type,Ch_names$Enzyme,Ch_names$Methodology,Ch_names$Read_number,Ch_names$parameters, sep = "_")

#####Create a transposed object of the samples names
S_names_t<- as.character(t(S_names))

####Quantify number of pairs of individual+duplicate
number_pairs<-length(S_names_t)/2

#####Create an object from Ch_names characteristics x number of pairs
Ch_names_error<-do.call("rbind", replicate(number_pairs, Ch_names, simplify = FALSE))

#####Create a name for each pair of samples
dup <- grep("*_d*",S_names_t,value=TRUE) # get the replicates (ending with _d or _id)
samps <- match(sub("*_d*","",dup),S_names_t) # match against its sample (ie names w/o _d or _id)
samps<- S_names_t[samps]
pairs<-paste(samps,dup, sep = "_") 

#####List .vcf files inside the input directory (inputs)
filelist <- list.files(dirInput_vcf, pattern="filtered*", full.names=TRUE)

#####Read the .vcf into a list file 
ldf <- lapply(filelist, read.vcfR)

#####Define the name of each list element with the object defined as combination_names (C_names)
names(ldf)<-C_names

#####################Quantifying Read-depth###############
dp<-list()
for (i in 1:length(ldf)){
  dp[[i]] <- extract.gt(ldf[[i]], element='DP', as.numeric=TRUE)
}

#####Define the name of each list element with the object defined as combination_names (C_names)
names(dp)<-C_names

#####Reorganize read-depth values to plot it
dpf<-list()
for (i in 1:length(ldf)){
  dpf[[i]] <- melt(dp[[i]], varnames=c("Index", "Sample"), value.name = "Depth", na.rm=TRUE)
  dpf[[i]] <- dpf[[i]][ dpf[[i]]$Depth > 0,] ##Deleted zero values
}

#####Define the name of each list element with the object defined as combination_names (C_names)
names(dpf)<-C_names


#####Plot read-depth for each sample for the whole set of combinations
#Create directory to store the plots
dir.create(path = "DEPTH_PLOTS")

mypath<-list()
for (i in 1:length(dpf)){
  
  mypath[[i]] <- file.path("DEPTH_PLOTS",paste("Depth_", C_names[i], ".jpg", sep = ""))
  
  jpeg(file=mypath[[i]], width = 25, height = 10, units = "cm", res = 400)
  print(ggplot(dpf[[i]], aes(x=Sample, y=Depth)) + geom_boxplot(fill="#C0C0C0") + theme_bw() +
          ylab("Read Depth (DP)") + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
                                          axis.text.y = element_text(angle = 0, size = 10))+
          ggtitle(C_names[i])+ ylim(0,800))
  dev.off()
}


#####Average of read-depth for each pair of individuals for the whole dataset combinations

average_perpair<-function(M){
  average<-matrix(NA,nrow = nrow(M),ncol = ncol(M))
  for (i in 1:nrow(M)){
    for (j in 1:ncol(M)){
      if(j%%2!=0) {
        average[i,j] <- (M[i,j]+M[i,j+1])/2
      }    
    }  
  }
  average<-average[, colSums(is.na(average)) != nrow(average)]
  average
}

dp_average_pair<-lapply(dp, average_perpair)

#####Average of read-depth for each pair and for the whole locus of individuals for the whole dataset combinations

average_perlocus<-function(M){
  average<-matrix(NA,nrow = 1,ncol = ncol(M))
  for (i in 1:nrow(M)){
    for (j in 1:ncol(M)){
      average[,j] <- (sum(M[,j]))/nrow(M)
    }  
  }
  colnames(M)<-pairs
  average
}  

dp_average_pair_locus<-lapply(dp_average_pair, average_perlocus)


#####Set the names of the pair samples
for (i in 1:length(dp_average_pair_locus)){
  colnames(dp_average_pair_locus[[i]])<-pairs
}

dp_average_pair_locus_V<-lapply(dp_average_pair_locus, as.vector)


dp_average_pair_locus_df<-do.call("rbind", dp_average_pair_locus_V)
colnames(dp_average_pair_locus_df)<-pairs

#####Transform into a vector
dp_average_pair_locus_vector<-c(dp_average_pair_locus_df)

#####Transform from .vcf format into genlight objects
genlight<- lapply(ldf, vcfR2genlight)

#####Dendograms for all samples for each combination

#####Set the function to create dendograms
dendrofunction<-function(x) nj(dist(as.matrix(x)))

#####define names trees from parameter name file
dendro_names <- C_names

#####Removes NAs from genlight objects
genlight_b<-list()
for(i in 1:length(genlight)){
  #To find the culplrits:
  toRemove <- is.na(glMean(genlight[[i]], alleleAsUnit = FALSE)) # TRUE where NA
  which(toRemove) # position of entirely non-typed loci
  genlight_b[[i]] <- genlight[[i]][, !toRemove]
}

#####Apply the dendrogram function in all the genlight objects
dendro_all<-lapply(genlight_b, dendrofunction)


#####Calculate the boots for each combination
myBoots<-list()
for (i in 1:length(genlight_b)){
  myBoots[[i]] <- boot.phylo(dendro_all[[i]], genlight_b[[i]], dendrofunction, mc.cores = 4)
}

for (i in 1:length(genlight_b)){
  dendro_all[[i]]$node.label<-myBoots[[i]]  
}

#####Plot the dendrograms in .jpeg format individual files in the DENDRO_PLOTS
#Create directory to store the plots
dir.create(path = "DENDRO_PLOTS")

for (i in 1:length(genlight_b)){
  
  mypath <- file.path("DENDRO_PLOTS",paste("TREE_", dendro_names[i], ".jpg", sep = ""))
  
  jpeg(file=mypath)
  mytitle = (dendro_names[i])
  plot(dendro_all[[i]], type="phylogram", edge.w =1, show.node.label=T, main = mytitle, cex =1.2)
  add.scale.bar(cex = 1.2, font = 5, col = "black")
  dev.off()
  
}

#####Locus/SNP error
#####Quantify of locus-error-rate and SNP-error-rate considering read-depth
#####Load Mastretta-Yanes et al (2014) functions
source("Mastretta_functions/LociAllele_error_pyrad.R")
source("Mastretta_functions/SNPs_error.R")


Locus_error_rate<-data.frame()
for(i in 1:length(genlight)){
  x<-LociAllele_error_pyrad(genlight[[i]], names(genlight)[[i]])
  Locus_error_rate<-rbind(Locus_error_rate, data.frame(x))
}

SNP_error_rate<-data.frame()
for(i in 1:length(genlight)){
  x<-SNP_error(genlight[[i]], names(genlight)[[i]])
  SNP_error_rate<- rbind(SNP_error_rate, data.frame(x))
}

#####Add parameter names column from Locus-error-rate dataframe to SNP-error-rate 
SNP_error_rate<-cbind(Locus_error_rate[,1],SNP_error_rate)

#####Set column names 
colnames(SNP_error_rate)<-c("parameter","pair", "SNP.error.rate")

#####Change "-" in the pair names column to "_" 
changeguionLAE<-gsub("-", "_", Locus_error_rate$pair)
Locus_error_rate$pair<-changeguionLAE

changeguionSE<-gsub("-","_", SNP_error_rate$pair)
SNP_error_rate$pair<-changeguionSE

##############################LOCUS_ERROR_RATE##########################

#####Order by name of pairs the data frame
Locus_error_rate<- plyr::arrange(Locus_error_rate,pair)

#####merge characteristics info with locus-error-rate
Locus_error_rate<-cbind(Locus_error_rate, Ch_names_error)

#####Add average read-depth into Locus-error-rate and SNP-error-rate objects
Locus_error_rate[,18]<-dp_average_pair_locus_vector

#####Set column names
colnames(Locus_error_rate)<-c( "parameter", "pair", "nloci", "nMissLoc", "MissTotProp", "shareMissLoc", "loci.mismatches",
                               "unshareMissLoc", "loci.error.rate", "n.loci.woNA", "allele.mismatches",
                               "allele.error.rate", "Read_type","Enzyme","Methodology", "Read_number", "parameters", 
                               "dp_average_pair_locus")


#####Transform loci.error.rate into numeric argument to operate with the locus-error-rate values
Locus_error_rate$loci.error.rate<-as.numeric(Locus_error_rate$loci.error.rate)

#####Divide locus-error-rate values for each pair by the average of the read-depth 
Locus_error_rate[,19]<-Locus_error_rate$loci.error.rate/Locus_error_rate$dp_average_pair_locus

#####Set column names
colnames(Locus_error_rate)<-c( "parameter", "pair", "nloci", "nMissLoc", "MissTotProp", "shareMissLoc", "loci.mismatches",
                               "unshareMissLoc", "loci.error.rate", "n.loci.woNA", "allele.mismatches",
                               "allele.error.rate", "Read_type","Enzyme","Methodology", "Read_number", "parameters", 
                               "dp_average_pair_locus", "loci.error.rate_depth")


#####Filter by number of reads and by methodology
RN<-unique(Ch_names_error$Read_number)
Meth<-unique(Ch_names_error$Methodology)

All_possible_comb<-expand.grid(RN,Meth)

colnames(All_possible_comb)<-c("Read_number", "Methodology")

LocusError<-list()
for (i in 1:nrow(All_possible_comb)) {
  LocusError[[i]]<-dplyr::filter(Locus_error_rate, Read_number==as.character(All_possible_comb[i,1]), 
                          Methodology==as.character(All_possible_comb[i,2]))
}

##############################SNP_ERROR_RATE##########################
#####order by name of pairs the data frame
SNP_error_rate<-dplyr::arrange(SNP_error_rate,pair)

#####merge characteristics info with SNP-error-rate
SNP_error_rate<-cbind(SNP_error_rate, Ch_names_error)

#####Add average read-depth into Locus-error-rate and SNP-error-rate objects
SNP_error_rate[,9]<-dp_average_pair_locus_vector

#####Set column names
colnames(SNP_error_rate)<-c( "parameter","pair","SNP.error.rate" ,"Read_type", "Enzyme", "Methodology", "Read_number","parameters",   
                             "dp_average_pair_locus")

#####Transform SNP.error.rate into numeric argument to operate with the SNP-error-rate values
SNP_error_rate$SNP.error.rate<-as.character(SNP_error_rate$SNP.error.rate)
SNP_error_rate$SNP.error.rate<-as.numeric(SNP_error_rate$SNP.error.rate)

#####Divide locus-error-rate values for each pair by the average of the read-depth 
SNP_error_rate[,10]<-SNP_error_rate$SNP.error.rate/SNP_error_rate$dp_average_pair_locus

#####Set column names
colnames(SNP_error_rate)<-c( "parameter","pair","SNP.error.rate" ,"Read_type", "Enzyme", "Methodology", "Read_number","parameters",   
                             "dp_average_pair_locus","SNP.error.rate_depth")


#####Filter by number of reads and by methodology
colnames(All_possible_comb)<-c("Read_number", "Methodology")

SNPsError<-list()
for (i in 1:nrow(All_possible_comb)) {
  SNPsError[[i]]<-dplyr::filter(SNP_error_rate, Read_number==as.character(All_possible_comb[i,1]), 
                         Methodology==as.character(All_possible_comb[i,2]))
  
}

#####PLOT LOCUS/SNP ERROR RATE#######
plotlocuserror<-list()
plotSNPerror<-list()
for (i in 1:nrow(All_possible_comb)) { 
  plotlocuserror[[i]]<-ggplot(data = LocusError[[i]])+
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "darkgrey"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "darkgrey"))+
    geom_bar(mapping = aes(x=parameters, y=loci.error.rate_depth, fill = Read_type), position = "dodge", stat = "identity")  +
    theme(axis.text.x = element_text( angle = 20, size = 9),axis.text.y  = element_text( angle = 0, size = 10))+
    guides(fill=guide_legend("Read type"))+
    ylab("Locus error rate")+
    ylim(0,0.04)
  
  plotSNPerror[[i]]<-ggplot(data = SNPsError[[i]])+
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "darkgrey"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "darkgrey"))+
    geom_bar(mapping = aes(x=parameters, y=SNP.error.rate_depth, fill = Read_type), position = "dodge", stat = "identity")  +
    theme(axis.text.x = element_text( angle = 20, size = 9),axis.text.y  = element_text( angle = 0, size = 10))+
    guides(fill=guide_legend("Read type"))+
    ylab("SNP error rate")+
    ylim(0,0.002)
  
}


####Create the directory to store the plots
dir.create(path = "locusError_PLOTS")
dir.create(path = "SNPsError_PLOTS")

mypath<-list()
for (i in 1:length(plotlocuserror)){
  
  mypath[[i]] <- file.path("locusError_PLOTS",paste("LocusError_",
                      paste(All_possible_comb[i,1],All_possible_comb[i,2], sep = "_"), ".jpg", sep = ""))
  
  jpeg(file=mypath[[i]], width = 25, height = 10, units = "cm", res = 400)
  print(plotlocuserror[[i]])
  dev.off()
}


mypath<-list()
for (i in 1:length(plotSNPerror)){
  
  mypath[[i]] <- file.path("SNPsError_PLOTS",paste("SNPError_",
                      paste(All_possible_comb[i,1],All_possible_comb[i,2], sep = "_"), ".jpg", sep = ""))
  
  jpeg(file=mypath[[i]], width = 25, height = 10, units = "cm", res = 400)
  print(plotSNPerror[[i]])
  dev.off()
}


#####Filtered by number of reads, type of reads and by methodology to plot errors by pair of samples
RN<-unique(Ch_names_error$Read_number)
RT<-unique(Ch_names_error$Read_type)
Meth<-unique(Ch_names_error$Methodology)

All_possible_comb_b<-expand.grid(RN,RT,Meth)

colnames(All_possible_comb_b)<-c("Read_number","Read_type","Methodology")

LocusError_b<-list()
for (i in 1:nrow(All_possible_comb_b)) {
  LocusError_b[[i]]<-dplyr::filter(Locus_error_rate, Read_number==as.character(All_possible_comb_b[i,1]), 
                            Read_type==as.character(All_possible_comb_b[i,2]),
                            Methodology==as.character(All_possible_comb_b[i,3]))
  
}

SNPsError_b<-list()
for (i in 1:nrow(All_possible_comb_b)) {
  SNPsError_b[[i]]<-dplyr::filter(SNP_error_rate, Read_number==as.character(All_possible_comb_b[i,1]), 
                           Read_type==as.character(All_possible_comb_b[i,2]),
                           Methodology==as.character(All_possible_comb_b[i,3]))
  
}

#####Store the plots of the error rates per sample for all the combinations to compare in more detail 

plotlocuserror_b<-list()
plotSNPerror_b<-list()
for (i in 1:nrow(All_possible_comb_b)) { 
  plotlocuserror_b[[i]]<-ggplot(data = LocusError_b[[i]])+
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "darkgrey"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "darkgrey"))+
    geom_bar(mapping = aes(x=parameters, y=loci.error.rate_depth, fill = pair), position = "dodge", stat = "identity")  +
    theme(axis.text.x = element_text( angle = 0, size = 12),axis.text.y  = element_text( angle = 0, size = 10))+
    guides(fill=guide_legend("Samples"))+
    ylab("Locus error rate")+xlab("parameters")+
    ylim(0,0.004)
  
  plotSNPerror_b[[i]]<-ggplot(data = SNPsError_b[[i]])+
    theme(panel.background = element_rect(fill = "white", colour = "white", size = 0.5, linetype = "solid"), 
          panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "darkgrey"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "darkgrey"))+
    geom_bar(mapping = aes(x=parameters, y=SNP.error.rate_depth, fill = pair), position = "dodge", stat = "identity")  +
    theme(axis.text.x = element_text( angle = 0, size = 12),axis.text.y  = element_text( angle = 0, size = 10))+
    guides(fill=guide_legend("Samples"))+
    ylab("SNP error rate")+xlab("parameters")+
    ylim(0,0.00045)
}  


####Create the directory to store the plots
dir.create(path = "locusError_b_PLOTS")
dir.create(path = "SNPsError_b_PLOTS")

mypath<-list()
for (i in 1:length(plotlocuserror_b)){
  
  mypath[[i]] <- file.path("locusError_b_PLOTS",paste("LocusError_b_", 
                  paste(All_possible_comb_b[i,1],All_possible_comb_b[i,2],All_possible_comb_b[i,3] ,sep = "_"), ".jpg", sep = ""))
  
  jpeg(file=mypath[[i]], width = 25, height = 10, units = "cm", res = 400)
  print(plotlocuserror_b[[i]])
  dev.off()
}


mypath<-list()
for (i in 1:length(plotSNPerror_b)){
  
  mypath[[i]] <- file.path("SNPsError_b_PLOTS",paste("SNPError_b_", 
                      paste(All_possible_comb_b[i,1],All_possible_comb_b[i,2],All_possible_comb_b[i,3], sep = "_"), ".jpg", sep = ""))
  
  jpeg(file=mypath[[i]], width = 25, height = 10, units = "cm", res = 400)
  print(plotSNPerror_b[[i]])
  dev.off()
}
