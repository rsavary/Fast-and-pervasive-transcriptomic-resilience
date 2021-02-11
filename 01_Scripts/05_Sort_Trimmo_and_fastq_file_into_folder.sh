!#/bin/bash
echo -e "\n"
echo "#########################################################################################################################################"
echo "#   Sort clean reads into folder and reads cleaning results into folder ::: Romain Savary 17.04.2019 #"
echo "#########################################################################################################################################"
echo -e "\n"
echo "#########################################################################################"
echo "################### RNAseq project 02 ::: 02_RSS-long_vs_CBASS-short ####################"
echo "#########################################################################################"
echo -e "\n"
echo -e "\n"
echo "RUN in lane folder, lane per lane, start with lane 1"



######## sort file into folder
for folder in $(ls -d */); do echo "-> "$folder; cd $folder; mkdir T1_after_trimmo; mkdir FastQC_after_trimmo_$folder; mv *html FastQC_after_trimmo_$folder/; rm *unpaired*; mv *.T1.fq T1_after_trimmo/; mkdir gz; mv *gz gz/; rm Tru*; mv *zip FastQC_after_trimmo_$folder/;  cd ..; done
