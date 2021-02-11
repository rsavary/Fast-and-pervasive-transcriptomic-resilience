
echo -e "\n"
echo "##############################################################################################################################"
echo "#   Script for pseudoaligning clean read with Kallisto (after Trimmomatic) *fq on the combined genes models of Smic and Spis  ::: Romain Savary 25.07.2019 #"
echo "##############################################################################################################################"
echo -e "\n"
echo "#########################################################################################"
echo "################### RNAseq project 02 ::: 02_RSS-long_vs_CBASS-short ####################"
echo "#########################################################################################"
echo -e "\n"
#####
echo "To run lane per lane inside the Lane* folder"
echo -e "\n"

rm -r */T2*

for folder in $(ls -d */); do echo '-> '$folder; cd $folder; foldername=$(echo $folder | sed 's/.$//'); runDir="T2_Kallisto_"$folder; mkdir $runDir; cd $runDir; q=$(pwd); sbatch /scratch/wally/PRTNR/EPFL/ENAC/ameibom/redseacoral/Wally_SCRIPTS/08_Kallisto_Smicro-Spist.sh ; cd ../../; done
