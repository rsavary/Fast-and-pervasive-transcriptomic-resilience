#!/bin/bash

#SBATCH --account=ameibom_redseacoral
#SBATCH --export=NONE
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 3:00:00
#SBATCH --job-name=fastqc
#SBATCH --export=None
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user r.savary@epfl.ch
#SBATCH --mail-type FAIL,END

module load Bioinformatics/Software/vital-it
module add UHTS/Quality_control/fastqc/0.11.7

fastqc *fq -f fastq

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
