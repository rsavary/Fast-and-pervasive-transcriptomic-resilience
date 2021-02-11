#!/bin/bash

#SBATCH --account=ameibom_redseacoral
#SBATCH --export=NONE
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 1:30:00
#SBATCH --job-name=Kallisto-Smicro-Spist
#SBATCH --export=None
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user r.savary@epfl.ch
#SBATCH --mail-type FAIL,END


module load Bioinformatics/Software/vital-it
module add UHTS/Analysis/kallisto/0.44.0

folder=$(pwd | cut -d "/" -f11)
q=$(pwd)

kallisto quant --index=/scratch/wally/PRTNR/EPFL/ENAC/ameibom/redseacoral/GENOME/03_combined_transcript_file_Spist_Smicro/Smicro_Spist_transcript_kallisto_index --output-dir=$q/results$folder --threads=$SLURM_CPUS_PER_TASK --plaintext ../T1*/*R1*.fq ../T1*/*R2*.fq
