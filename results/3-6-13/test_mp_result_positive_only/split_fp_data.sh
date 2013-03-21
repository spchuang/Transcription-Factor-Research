#this script reads in the all.footprint file and split them to different files according to to each chromosome

clear
#declare the string of chromosomes
chromsomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)

echo "start splitting the footprint file..."
mkdir chr.footprints;
cd chr.footprints;

cat all.footprints
for (( i = 0; i<= 22; i++ ))
   do
      echo "dumpling fp for ${chromsomes[$i]}" ;
      cat "../all.footprints" | grep "${chromsomes[$i]}$(echo -ne \\t)" > "${chromsomes[$i]}.footprints";
   done
