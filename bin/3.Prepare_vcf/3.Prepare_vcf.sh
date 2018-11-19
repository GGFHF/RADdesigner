#!/bin/bash


#################################################################
#                                                              #    
#         VCF-CONTRAST tool to organize the .vcf files         #
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

#This script contains a custom bash script to organize the .vcf files for better files management in the next step.

#-------------------------------------------------------------------------------                     

#####Set path to the software VCFTOOLS
dirVCFTOOLS=/PATH_TO_VCFTOOLS/perl
export PERL5LIB=/PATH_TO_VCFTOOLS/perl
#####SET PATH TO YOUR WORKING DIRECTORY (where you have stored the directory 3.Prepare_vcf) 
dir=/PATH_TO_WORKING_DIRECTORY/bin/3.Prepare_vcf  #Select the path to your working directory for the third script
#####SET PATH TO vcf filtered files (which are in the 1.Filtering_multiallelicSNPs directory)
dirInput_vcf=/PATH_TO_WORKING_DIRECTORY/bin/1.Filtering_multiallelicSNPs/FILTERED_OUTPUT
#####SET PATH TO YOUR INPUT FILES
dirInput=/PATH_TO_WORKING_DIRECTORY/Input_files 
#####LOAD A .txt FILE FORMAT WITH YOUR SAMPLES LISTED (see Samples_names.txt example in the Input_files directory)
SAMPLE_FILE=$dirInput/Samples_names.txt

#####Create and set path to store the output files
mkdir OUTPUT_txt 
DIRoutput=$dir/OUTPUT_txt

#####Create the individual samples list from the samples names file

INDIVIDUALS=""

while read LINE; do
	if [ "$INDIVIDUALS" == "" ]; then 
		INDIVIDUALS=$LINE
	else
		INDIVIDUALS=$INDIVIDUALS","$LINE
	fi
done < $SAMPLE_FILE

#####Move to input directory
cd $dirInput_vcf

#####Convert and prepare the .vcf files to use it in the next step (4.Filtering_locusSNPerror)
ls *.vcf > vcf_files.txt

while read FILE; do
	$dirVCFTOOLS/vcf-contrast +"$INDIVIDUALS" -"$INDIVIDUALS" $FILE | $dirVCFTOOLS/vcf-query -f '%CHROM\t%POS[\t%GTR]\n' > $DIRoutput/vcf_contrast_"$FILE".txt
		if [ $? -ne 0 ]; then echo 'Script ended with errors.'; exit 1; fi
done < vcf_files.txt


#-------------------------------------------------------------------------------

# End
echo '*********************END*****************************'
exit 0

#-------------------------------------------------------------------------------





