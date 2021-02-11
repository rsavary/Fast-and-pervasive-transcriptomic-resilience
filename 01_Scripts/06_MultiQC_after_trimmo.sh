!#/bin/bash
echo -e "\n"
echo "#########################################################################################################################################"
echo "#   Combine with MultiQC reads quality check results after trimmomatic ::: Romain Savary 18.05.2019 #"
echo "#########################################################################################################################################"
echo -e "\n"
echo "#########################################################################################"
echo "################### RNAseq project 02 ::: 02_RSS-long_vs_CBASS-short ####################"
echo "#########################################################################################"
echo -e "\n"
echo -e "\n"
echo "RUN in top folde, CBASS vs RSS"

module load Bioinformatics/Software/vital-it
module add UHTS/Analysis/MultiQC/1.7;

multiqc */*/*after*/*fastqc.zip -d -o multiqc_after_trimo
