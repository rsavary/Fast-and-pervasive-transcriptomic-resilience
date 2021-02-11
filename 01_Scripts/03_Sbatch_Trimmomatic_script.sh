#!/bin/bash

echo -e "\n"
echo "#########################################################################################################################################"
echo "#   Wally-Slurm Script submit reads cleaning to cluster (Trimmomatic) *fq ::: Romain Savary 12.07.2019 #"
echo "#########################################################################################################################################"
echo -e "\n"
echo "#########################################################################################"
echo "################### RNAseq project 02 ::: 02_RSS-long_vs_CBASS-short ####################"
echo "#########################################################################################"
echo -e "\n"
echo -e "\n"
echo "RUN in lane folder, lane per lane, start with lane 1"



a=0;for folder in $(ls -d */); do cd $folder; cp -s /software/UHTS/Analysis/trimmomatic/0.36/adapters/TruSeq3-PE-2.fa .;
	sbatch /scratch/wally/PRTNR/EPFL/ENAC/ameibom/redseacoral/Wally_SCRIPTS/03_Trimmomatic.sh;
	cd .. ;done

