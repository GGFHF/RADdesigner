#!/bin/bash

#################################################################
#                                                              #    
#               Multi-allelic SNPs filtering                   #
#                         bash script                          #
#                                                              #
#################################################################

#-------------------------------------------------------------------------------

#This custom script has been developed by:

#       GI Sistemas Naturales e Historia Forestal(formerly known as GI Genética, Fisiología e Historia Forestal)
#       Dpto. Sistemas y Recursos Naturales
#       ETSI Montes, Forestal y del Medio Natural
#       Universidad Politécnica de Madrid
#       https://github.com/ggfhf/
#
#   Licence: GNU General Public Licence Version 3


#-------------------------------------------------------------------------------

#This script contains a custom workflow to filter milti-allelicSNPs from vcf files.

#-------------------------------------------------------------------------------                     
#####SET PATH TO YOUR WORKING DIRECTORY (where you have stored the directory 1.Filtering_multiallelicSNPs) 
dir="/PATH_TO_WORKING_DIRECTORY/bin/1.Filtering_multiallelicSNPs"  #Select the path to your working directory for the first script
#####SET PATH TO YOUR INPUT FILES
dirInput="/PATH_TO_WORKING_DIRECTORY/Input_files" ##Same for all the scripts

#####CREATE IpyRAD_output_vcf (where are stored the vcf files) DIRECTORY INSIDE DIRECTORY Input_files 
#####Set path to vcf input files
dirvcf=$dirInput/INPUT_VCF_DATA

#####Create and set path to the bi-allelic output directory
mkdir $dir/FILTERED_OUTPUT
outdir=$dir/FILTERED_OUTPUT

#####Create and set path to the multi-allelic output directory
mkdir $dir/Multi_allelic
multi_allelic=$dir/Multi_allelic

#####Create and set path to the bi-allelic output directory
mkdir $dir/Bi_allelic
bi_allelic=$dir/Bi_allelic

####Move to IpyRAD .vcf output directory
cd $dirvcf

###Save the header (first ten rows) of the .vcf files in new .txt files called header_*name_of_the_file.txt 
echo "Saving the header of vcf files in header*.txt files"

ls *.vcf > IpyRAD_output_vcf.txt

while read FILE; do

	head $FILE > $dir/header"$FILE".txt
	if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi

done < IpyRAD_output_vcf.txt


####Move to initial working directory
cd $dir

#Installing required R packages

#echo "Installing required R packages"

#Rscript multi_bi_allelic_filtered.R             #Installing required R packages


###Open the .vcf files without header in R studio software. 
#multi_bi_allelic_filtered.R Rscript separates the multi-allelic SNPs and bi-allelic SNPs in multi_allelic and bi_allelic files, respectively.

echo "Calling the multi_bi_allelic_filtered.R Rscript to separate the multi- and bi-allelic SNPs."

Rscript multi_bi_allelic_filtered.R             #Calling the multi_bi_allelic_filtered.R script

###Paste each header to the bi_allelic*.vcf files
echo "Pasting each header to the bi_allelic*.vcf files"


ls *bi_allelic.vcf > bi_allelic_files.txt

ls header*.txt > header_files.txt


while read FILE1; do
	while read FILE2; do
		cat $FILE2 $FILE1 > filtered_"$FILE1"
		if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi
	done < header_files.txt
done < bi_allelic_files.txt


###Move filtered files to FILTERED_OUTPUT directory
echo "Moving the filtered files to the output directory (FILTERED_OUTPUT)"

mv filtered*.vcf $outdir

###Move the ouput bi_allelic files to Bi_allelic directory
echo "Moving the ouput directory for Bi-allelic vcf files"

mv *bi_allelic.vcf $bi_allelic

###Move the ouput multi_allelic files to Multi_allelic directory
echo "Moving the ouput directory for Multi_allelic vcf files"

mv *multi_allelic.vcf $multi_allelic

###Remove of the header*.txt files an the list of bi-allelic files
echo "Removing the header*.txt files, bi_allelic_files.txt file and IpyRAD_output_vcf.txt file "
rm *header*.txt
rm *bi_allelic_files.txt
rm $dirvcf/IpyRAD_output_vcf.txt

###Open the .vcf files filtered in R studio software. 
##This R script contains a custom workflow to plot the percentage of multi-allelic SNPs from vcf filtered files, and to plot the number of loci/SNPs post multi-allelic filtered

echo "Calling the plot_percentage_multi_allelicSNPs.R Rscript to plot the percentage of multi-allelic SNPs from vcf filtered files and to plot the number of loci/SNPs post multi-allelic filtered."

Rscript plot_percentage_multi_allelicSNPs.R             #Calling the plot_percentage_multi_allelicSNPs.R script

#-------------------------------------------------------------------------------

# End
echo '*********************END*****************************'
exit 0

#-------------------------------------------------------------------------------
