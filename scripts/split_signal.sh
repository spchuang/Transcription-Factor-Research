#this script reads in the signal data and then split into different chromosomes
clear

#all the cell types for signal data
celltypes=(AG10803 CD20+ GM12865 HAEpiC HCF HCPEpiC HepG2 HIPEpiC HMVEC-dBl-Ad HPdLF HRCEpiC HVMF NB4 NHDF-Ad NHLF SKMC Th1 AoAF GM06990 H7-hESC HA-h HCM HEEpiC HFF HMF HPAF HPF HSMM K562 NH-A NHDF-neo SAEC)

#declare the string of chromosomes
chromsomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)

echo "start splitting the signal file..."
for (( j = 0; j<= 33; j++ ))
	do
	    echo "Start on Cell Type: ${celltypes[$j]}" ;
		mkdir ${celltypes[$j]}_split;
		cd ${celltypes[$j]}_split;
		
		for (( i = 0; i<= 22; i++ ))
		   do
		      echo "dumpling signal for ${chromsomes[$i]}" ;
		      cat "../${celltypes[$j]}" | grep "${chromsomes[$i]}$(echo -ne \\t)" > "${chromsomes[$i]}.${celltypes[$j]}";
		   done
		   
		cd ..
	
	done