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
#SBATCH --job-name=Trimmomatic
#SBATCH --export=None
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user r.savary@epfl.ch
#SBATCH --mail-type END, FAIL

module load Bioinformatics/Software/vital-it
module add UHTS/Analysis/trimmomatic/0.36

file=$(ls *_R1_*.fastq.gz)
Basename1=$(echo $file | cut -d'.' -f1 | cut -d'_' -f1,2);
Basename2=$(echo $file | cut -d'.' -f1 | cut -d'_' -f4); 
Secondpair=$(ls $Basename1'_R2_'$Basename2*.fastq.gz); echo $file' - '$Secondpair;


trimmomatic PE -phred33 -threads $SLURM_CPUS_PER_TASK $file $Secondpair $Basename1'_R1_'$Basename2'.T1.fq' $Basename1'_R1_'$Basename2'.T1.fq.unpaired' $Basename1'_R2_'$Basename2'.T1.fq' $Basename1'_R2_'$Basename2'.T1.fq.unpaired' ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40
