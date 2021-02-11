#!/bin/bash

# Slurm options:
#SBATCH --account=ameibom_redseacoral
#SBATCH --export=NONE
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 3:00:00
#SBATCH --job-name=Kallisto-index
#SBATCH --export=None
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user r.savary@epfl.ch
#SBATCH --mail-type END, FAIL

###### generate index of Spist transcripts
module load Bioinformatics/Software/vital-it
module add UHTS/Analysis/kallisto/0.44.0
kallisto index --make-unique -i /scratch/wally/PRTNR/EPFL/ENAC/ameibom/redseacoral/GENOME/03_combined_transcript_file_Spist_Smicro/Smicro_Spist_transcript_kallisto_index /scratch/wally/PRTNR/EPFL/ENAC/ameibom/redseacoral/GENOME/03_combined_transcript_file_Spist_Smicro/Smic_Spist_combined_CDS_file.fa

