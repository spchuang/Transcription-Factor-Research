#this script reads in the K562Sig_filter.wig 
clear
#declare the string of chromosomes
chromsomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)

echo "start splitting the signal file..."
mkdir motif_non_consbased;
cd motif_non_consbased;

#cat hg19_instances-thresh8-0.4.txt
for (( i = 0; i<= 22; i++ ))
   do
      echo "dumpling signal for ${chromsomes[$i]}" ;
      cat "../0.0" | grep -v _C |  grep "${chromsomes[$i]}$(echo -e ' ')" > "${chromsomes[$i]}.motif";
   done
