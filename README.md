README

### `/bin/`
Contains the workflow scripts for the sections *Data post-filtering and error rates* of the manuscript "RADdesigner: a workflow to select the optimal sequencing methodology in genotyping experiments on woody plant species".

The scripts in `/bin` should be run in the order they are numbered. R functions used by some of these scripts are not numbered.
The parameters that need to be modified in each script are indicated at the top of the script. Comments with capital letters
indicate the command lines that need to be modified to include path to specific applications or to the user's working directory.


Scripts content:

* `0.R_dependendies_installation.R` R scritp to install and call the packages required to run subsequent R scripts.
* `VCFtools_installation_instructions.txt` Text file with the command line instructions to install VCFtools.
* `1.Filtering_multiallelicSNPs.sh` Bash script to filter multi-allelicSNPs from .vcf files, inside this script we call two R scripts:
* `multi_bi_allelic_filtered.R` R script to separate multi-allelic from bi-allelic SNPs in each vcf.file
* `plot_percentage_multi_allelicSNPs.R` R script to plot the percentage of multi-allelic SNPs filtered for each combination. Also this script plots the number of filtered loci/SNPs after filtering.
* `2.Depth_locusSNPerror.R` R script that contains a custom workflow to quantify and plot read-depth and both locus/SNP error rates. Also this script quantifies and plots the dendrograms of similarity for each sample and for each combination.
* `LociAllele_error_pyrad.R` Mastretta-Yanes et al. (2014), R function to quantify locus-error-rate.
* `SNPs_error.R` Mastretta-Yanes et al. (2014), R function to quantify SNP-error-rate.
* `3.Prepare_vcf.sh`contains a custom bash script to organize the .vcf files to ease file management in the next step.
* `4.Filtering_locusSNPerror.R`quantifies and filter the locus/SNP error.

These scripts use the vcf files data recovered from the IpyRAD software. 
Each script uses the output of the previous script as an input.

### `/Input_files/`
Contains the required information files to run the scripts placed in `/bin` directory.

* `INPUT_VCF_DATA` is a directory to stage the input .vcf files. Be sure to store your duplicate input vcf data using the string "_d". 
* `Ch_combinations.csv` is a file with the characteristics of the combinations.
* `Samples_names.txt`is a file with the names of your individuals and the duplicates placed in a column. 
