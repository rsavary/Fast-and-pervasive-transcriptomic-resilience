#!/bin/bash
echo -e "\n"
echo "#########################################################################################################################################"
echo "#   Wally-Slurm Script for reodering files into folder before reads cleaning *fastq.gz (Trimmomatic) *fq and quality control (fastQC)::: Romain Savary 12.07.2019 #"
echo "#########################################################################################################################################"
echo -e "\n"
echo "#########################################################################################"
echo "################### RNAseq project 02 ::: 02_RSS-long_vs_CBASS-short ####################"
echo "#########################################################################################"
echo -e "\n"
echo -e "\n"
echo "RUN in lane folder, lane per lane, start with lane 1"

####First from sequencing files move file from within version and have only Lane*/*gz and all gz files
#for i in $(ls -d Lane*);do cd $i; pwd; mv version*/*gz . ; cd ..; done

########
echo -e "\n"
echo "Change name of fastq files to have a simple name with Sample name _ lane _ Reads1-2 _ file"
echo -e "\n"
for i in $(ls *gz); do j=$(echo $i | cut -d "_" -f4,6,7,8); echo $i "=>" $j; mv $i $j; lane=$(echo $i | cut "_" -f6); done

######## move all fastq file in the RNAseq folder
echo "--> Number of fastq files";
ls | wc -l ;
echo "-->Path of the fastq files";
ls *fastq.gz ;


########## Make one folder per sample and move the samples in their respective folder
ARRAY=(); printf "%0.s-" {1..10};echo; echo "-> Filenames";for filename in $(find . -maxdepth 1 -type f | sed s,^./,,); 
do echo $filename; ARRAY+=($( echo $filename | cut -d'_' -f1)); done;  printf "%0.s-" {1..10};echo; echo "-> Folders"; 
sorted_unique_ARRAY=$(echo "${ARRAY[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '); for i in ${sorted_unique_ARRAY[*]}; 
do echo $i; mkdir $i; for files in $(find . -type f -maxdepth 1 -name $i"*" | sed s,^./,,); do echo $files; mv $files $i'/'; done; done
