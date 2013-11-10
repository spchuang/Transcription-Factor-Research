#this script reads in the all.footprint file and split them to different files according to to each chromosome

clear
#declare the string of chromosomes
chromsomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)
celltypes=(AG10803 AoAF CD20+ GM06990 GM12865 H7-hESC HAEpiC HA-h HCF HCM HCPEpiC HEEpiC HepG2 HFF HIPEpiC HMF HMVEC-dBl-Ad HPAF HPdLF HPF HRCEpiC HSMM HVMF K562 NB4 NH-A NHDF-Ad NHDF-neo NHLF SAEC SKMC Th1)

cd ../data
echo "start splitting the footprint file..."
mkdir footprintsData;
cd footprintsData;

cat all.footprints
for (( i = 0; i< 32; i++ ))
   do
   	  echo "dumpling fp for ${celltypes[$i]}" ;
      mkdir ${celltypes[i]}
      cd ${celltypes[i]}
      for (( j = 0; j< 23; j++ ))
         do
            echo "processing  ${chromsomes[$i]}" ;
            cat "../../all.footprints" | grep "${chromsomes[$j]}$(echo -ne \\t)" | grep "${celltypes[$i]}$(echo -ne \\t)" > "${chromsomes[$j]}";
         done
      cd ..
   done
