#include "readMotifMaster.cpp"

int main(){
	//step 1: read all the motif data for a given cell type
	//step 2: average the signals/cons for the same motif
	//setp3: output average signals for each motif 
	string cellTypes[] = {"AG10803", "AoAF", "CD20+", "GM06990", "GM12865","H7-hESC","HAEpiC","HA-h","HCF","HCM","HCPEpiC","HEEpiC","HepG2","HFF","HIPEpiC","HMF","HMVEC-dBl-Ad","HPAF","HPdLF","HPF","HRCEpiC","HSMM","HVMF","K562","NB4","NH-A","NHDF-Ad","NHDF-neo","NHLF","SAEC","SKMC","Th1"};
	string chromosomes[] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                           "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                           "chr21", "chr22","chrX"};
	int totalChr = 23;     //23
    int totalCell = 1; //31       
    string rootDir = "TEST"; 
    
    ReadMotifMaster readMaster(rootDir);
    
    for(int cell=0; cell<totalCell; cell++){
    	readMaster.readMotif(cellTypes[cell], chromosomes, totalChr);
		readMaster.averageMotifData(cellTypes[cell]);
		readMaster.printMotifData(cellTypes[cell]);
    }
    			
}