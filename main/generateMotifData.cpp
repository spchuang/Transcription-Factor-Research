#include "dataUtility.cpp"
/*
	Question: given a list of motifs (start,end), a list of fp of any celltype(start,end), filter otu the motifs that don't overlap with any fp
	step 1. aggregate all the motifs first lol 
		-> celltype/motif1,motif2,motif3... (non aggregated)
	if say, we can't hold all the motif instances and their signals of one cell in one run, then get the data for each chromosome then append each motif signal in the file. So open 456 file descriptors?
*/

int main(){
	
	string cellTypes[] = {"AG10803", "AoAF", "CD20+", "GM06990", "GM12865","H7-hESC","HAEpiC","HA-h","HCF","HCM","HCPEpiC","HEEpiC","HepG2","HFF","HIPEpiC","HMF","HMVEC-dBl-Ad","HPAF","HPdLF","HPF","HRCEpiC","HSMM","HVMF","K562","NB4","NH-A","NHDF-Ad","NHDF-neo","NHLF","SAEC","SKMC","Th1"};
	string chromosomes[] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                           "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                           "chr21", "chr22","chrX"};
    int totalChr = 23;    
    int totalCell = 31;  //31     
    string rootDir = "noFilterMotif";            
    
	MotifMaster motif_master;
	motif_master.setReadCells(cellTypes, 1);
	
	FootprintMaster fp_master;
	//foot print reads each chr for ALL cell types
	fp_master.setReadCells(cellTypes, totalCell);

	
	cout << "[DEBUG]LETS START THIS SHIT" <<endl;
	for(int chr=0; chr<totalChr; chr++){
		cout <<"[DEBUG]for motif chr " << chr<<endl;
		motif_master.clearSegments();
		motif_master.setReadChromosomes(chromosomes, chr, chr+1);
		motif_master.printProperties();
		motif_master.startReadingData();
		
		cout <<"[DEBUG]for fp chr " << chr<<endl;
		fp_master.setReadChromosomes(chromosomes,chr, chr+1);
		fp_master.startReadingData();
		fp_master.sortSegmentsBySegStart();
		
		
		//filter out motifs that are not included in any footprint
		
		cout <<"[DEBUG]Filter motif instances to those that includes at least one fp"<<endl;
		vector<segment>* fp_segs = fp_master.getSegment();
		vector<segment>* mt_segs = motif_master.getSegment();
		
		int mt_size = mt_segs->size();
		int fp_size = fp_segs->size();
		int j=0;
		/*
		vector<segment> new_mt;
		for(int i=0; i<mt_size; i++){
			bool containsFootprint = false;
			while(j<fp_size){
				//using segStart cuz we're matching the actual segment coverage and not the frame (which is user defined)
				int fp_start = (*fp_segs)[j].segStart;
				int fp_end   = fp_start + (*fp_segs)[j].segLength;
				int mt_start = (*mt_segs)[i].segStart;
				int mt_end   = mt_start + (*mt_segs)[i].segLength;
				
				if(fp_start > mt_end && fp_end >mt_end){
					break;
				}
				if(
		          ((fp_end >= mt_start+1) && (fp_end<=mt_end)) || 
		           ((fp_start >= mt_start) && (fp_start < mt_end))
		        ){
		        	//includes
					containsFootprint = true;
					break;
		        }
		        j++;
			}
			
			if(containsFootprint){
				new_mt.push_back((*mt_segs)[i]);
			}
		}
		
		motif_master.clearSegments();
		for(int i=0; i<new_mt.size(); i++){	
			mt_segs->push_back(new_mt[i]);
		}
		new_mt.clear();*/
		motif_master.sortSegmentsByFrameStart();

		int new_mt_size = motif_master.getSegSize();
		
		//cout <<"new size: " <<new_mt.size() <<endl;
		//cout <<"fp size: "<<fp_segs->size()<<endl;
		
		cout <<"[DEBUG]new motif size: "<<motif_master.getSegSize()<< " / " <<mt_size<<endl;
		
		
		string motif_list = rootDir+"/"+chromosomes[chr];
		ofstream motifListFile(motif_list.c_str(), ios_base::trunc );
		for(int i=0;i<mt_segs->size(); i++){
			motifListFile<< (*mt_segs)[i].name << " : " <<(*mt_segs)[i].frameStart << ". seg at "<<(*mt_segs)[i].segStart<<endl;
		}
		
		for(int cell = 0; cell<totalCell; cell++){
			cout <<"[DEBUG]start reading signal for motif at chr " << chromosomes[chr]<< " , " <<cellTypes[cell]<<endl;
			motif_master.readSignalData(cellTypes[cell], chromosomes[chr], 0);
			motif_master.flipNegativeSegments();	
			
			string dirname = rootDir+"/motif_signals/"+cellTypes[cell];
			cout <<"[DEBUG]output direcotry:" + dirname <<endl;
			struct stat st = {0};
			if (stat(dirname.c_str(), &st) == -1) {
				if (mkdir(dirname.c_str(), 0700) != 0){
					cout <<"PROBLEM FUCK.."<<endl;
				}
	            
				
			}else{
				//cout <<"PROBLEM AND SHIT?"<<endl;
			}
			
			string filename = dirname+"/"+chromosomes[chr];
			ofstream outputFile(filename.c_str(), ios_base::trunc );
			
			for(int i=0; i< new_mt_size; i++){

			    outputFile << (*mt_segs)[i].name << "  ";
			    outputFile << (*mt_segs)[i].segStart << " " <<(*mt_segs)[i].segStart+(*mt_segs)[i].segLength<<endl;
			             
			    for(int j=0; j<WINDOW_SIZE; j++){
			      outputFile <<(*mt_segs)[i].signal[j] << " " << (*mt_segs)[i].con_level[j] <<endl;
			    }
			} 
		}
	
		fp_master.clearSegments();
		
	}
	
	
	
	return 0;
}