#!/bin/bash
echo -e "\n"
echo "#########################################################################################################################################"
echo "#   Wally-Slurm Script move fastqc file into folder ::: Romain Savary 12.07.2019 #"
echo "#########################################################################################################################################"
echo -e "\n"
echo "#########################################################################################"
echo "################### RNAseq project 02 ::: 02_RSS-long_vs_CBASS-short ####################"
echo "#########################################################################################"
echo -e "\n"
echo -e "\n"
echo "RUN in lane folder, lane per lane, start with lane 1"

######## move fastqc html file in folder : FastQC_before_trimmo

a=0;for Lane in $(ls -d */);do cd $Lane; 
	for folder in $(ls -d */); do cd $folder; 
	mkdir  FastQC_before_trimmo_$folder; mv *html* FastQC_before_trimmo_$folder/; mv *zip FastQC_before_trimmo_$folder/; 
	cd ..; done
    cd ..; done
