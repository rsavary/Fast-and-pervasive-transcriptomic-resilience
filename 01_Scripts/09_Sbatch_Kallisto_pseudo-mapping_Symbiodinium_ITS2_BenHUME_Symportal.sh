
echo -e "\n"
echo "##############################################################################################################################"
echo "#   Script for pseudoaligning clean read with Kallisto (after Trimmomatic) *fq the Smic ITS data base of SymPortal 2019 ::: Romain Savary 25.07.2019 #"
echo "##############################################################################################################################"
echo -e "\n"
echo "#########################################################################################"
echo "################### RNAseq project 02 ::: 02_RSS-long_vs_CBASS-short ####################"
echo "#########################################################################################"
echo -e "\n"
#####
echo "To run lane per lane inside the Lane* folder"
echo -e "\n"

rm -r */T7*

for folder in $(ls -d */); do echo '-> '$folder; cd $folder; foldername=$(echo $folder | sed 's/.$//'); runDir="T7_Kallisto_1192_Symportal_its2_"$folder; mkdir $runDir; cd $runDir; q=$(pwd); sbatch /scratch/wally/PRTNR/EPFL/ENAC/ameibom/redseacoral/Wally_SCRIPTS/13_Kallisto_pseudo-mapping_Symbiodinium_ITS2_BenHUME_Symportal ; cd ../../; done
