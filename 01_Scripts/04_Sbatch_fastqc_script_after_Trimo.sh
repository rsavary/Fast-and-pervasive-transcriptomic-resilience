#!/bin/bash
echo -e "\n"
echo "#########################################################################################################################################"
echo "#   Wally-Slurm Script for reads cleaning after Trimmomatic ::: Romain Savary 12.07.2019 #"
echo "#########################################################################################################################################"
echo -e "\n"
echo "#########################################################################################"
echo "################### RNAseq project 02 ::: 02_RSS-long_vs_CBASS-short ####################"
echo "#########################################################################################"
echo -e "\n"
echo -e "\n"
echo "RUN in lane folder, lane per lane, start with lane 1"


###### submit jobs of the script Test_script_trimo.sh
a=0;for folder in $(ls -d */);
do cd $folder; sbatch /scratch/wally/PRTNR/EPFL/ENAC/ameibom/redseacoral/Wally_SCRIPTS/04_fastqc_script_after_Trimo.sh ; cd ..;
done
