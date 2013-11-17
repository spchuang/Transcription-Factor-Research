#include "dataUtility.cpp"

struct averageMotifSegment{
   float con_level[WINDOW_SIZE];
   float signal[WINDOW_SIZE];       
   char* name;
   int count;
};

class ReadMotifMaster{
	public:
		ReadMotifMaster(string dir);
		
		void readMotif(string cell, string* chr, int n);
		void averageMotifData(string cell);
		void printMotifData(string cell);
		
	private:
		vector<segment> segs;
		vector<averageMotifSegment> unique_motifs;
		string rootDir;
};

ReadMotifMaster::ReadMotifMaster(string dir){

	rootDir = dir;
}

void ReadMotifMaster::readMotif(string cell, string* chr, int n){
	cout <<"[DEBUG]Reaing motifs for " << cell <<endl;
	segs.clear();
	for(int i=0; i<n; i++){
		cout <<"[DEBUG]read chr " <<chr[i]<<endl;

		string motiFile = rootDir+"/motif_signals/"+cell+"/"+chr[i];

		ifstream motifInFile(motiFile.c_str());
		
		string motif_name;
		int s, e;
		
		while(motifInFile >> motif_name >> s >>e){
			//cout <<motif_name << " : " <<s <<" , " <<e<<endl;
			segment mtt;
			mtt.name= new char [motif_name.size()+1];
			strcpy(mtt.name, motif_name.c_str());
			mtt.segStart  = s;
			mtt.segLength = e-s;
			mtt.startIndex =0;
			int signal;
			string cons_value;
			int i=0;
			while(motifInFile >> signal >> cons_value){
				mtt.signal[i] = signal;
				mtt.con_level[i] = atof(cons_value.c_str());
				i++;
				if(i == WINDOW_SIZE) break;
			}
			//print to debug
			//cout <<mtt.name << " : " <<mtt.segStart<< " , " <<mtt.segLength<<endl;
			for(i=0; i< FRAME_SIZE; i++){
				//cout << mtt.signal[i] << " " << mtt.con_level[i] <<endl;
			}
			segs.push_back(mtt);
		}
			
	}
	cout <<"[DEBUG]read total of motif: " <<segs.size()<<endl; 	
	
	
}

bool sortByCount(const averageMotifSegment &a, const averageMotifSegment &b)
{
    return a.count < b.count;
}


void ReadMotifMaster::averageMotifData(string cell){
	cout <<"[DEBUG]Average all motifs " <<endl;
	//find the number of unique motifs first?
	unique_motifs.clear();


	for(int i=0; i<segs.size(); i++){
		//cout <<(*mt)[i].motifName<<endl;
	    bool exists = false;
	    
	    for(int j=0; j<unique_motifs.size(); j++){
	    	if(strcmp(segs[i].name, unique_motifs[j].name) == 0){
				exists = true;
			}
		}
    
		//create this many centroids
		if(!exists){
			//cout << unique_motifs.size() <<" : " << segs[i].name <<endl;
			size_t c=0;
			while(segs[i].name[c] != 0)
				c++;
			averageMotifSegment mtt;
			for(int j=0; j< WINDOW_SIZE; j++){
				mtt.signal[j] = 0;
				mtt.con_level[j] = 0;
			}
			mtt.name = new char [c];
			mtt.count = 0;
			strcpy(mtt.name, segs[i].name);
			unique_motifs.push_back(mtt);
		}
	}
	cout <<"[DEBUG]total unique motifs: " <<unique_motifs.size()<<endl;
	
	cout<<"[DEBUG]aggregate signal and conservation signals"<<endl;
	//add the sig and con to centroids
	for(int i=0; i<unique_motifs.size(); i++){
		for(int k=0; k<segs.size(); k++){
			if(strcmp(segs[k].name, unique_motifs[i].name) == 0){
				int startIndex      = segs[k].startIndex;
				for(int j=0; j<WINDOW_SIZE; j++){
					unique_motifs[i].signal[j] += segs[k].signal[j+startIndex];
					unique_motifs[i].con_level[j] += segs[k].con_level[j+startIndex];
				}
				unique_motifs[i].count++;
			}
		}
	
	}
	cout<<"[DEBUG]average the aggregation" <<endl;
	//average the signals and reposition the centroid
	int total_count =0;
	for(int i=0; i<unique_motifs.size(); i++){
		for(int j=0; j<WINDOW_SIZE; j++){
			unique_motifs[i].signal[j] = (float)(unique_motifs[i].signal[j] / unique_motifs[i].count);
			unique_motifs[i].con_level[j] = (float)(unique_motifs[i].con_level[j] / unique_motifs[i].count);
		}
		total_count+=unique_motifs[i].count;
	}
	cout <<"TOTAL COUNT:  "<<total_count<<endl;

	
	cout <<"[DEBUG]sort"<<endl;
	sort (unique_motifs.begin(), unique_motifs.end(), sortByCount); 
	
	

}

void ReadMotifMaster::printMotifData(string cell){
	cout <<"[DEBUG]output result" <<endl;
	struct stat st = {0};
	
	//create directory if it doesn't exists
	string dirname = rootDir+"/average/"+cell;
	
	cout <<"[DEBUG]output direcotry:" + dirname <<endl;
	if (stat(dirname.c_str(), &st) == -1) {
		if (mkdir(dirname.c_str(), 0700) != 0){
			cout <<"PROBLEM FUCK.."<<endl;
		}
	}
	
	string filename = dirname+"/fpsig";
	ofstream outputFile(filename.c_str(), ios_base::trunc );
	filename = dirname+"/consSig";
	ofstream consFile(filename.c_str(), ios_base::trunc );
	filename = dirname+"/motifOrder";
	ofstream motifFile(filename.c_str(), ios_base::trunc );
	for(int i=0; i<unique_motifs.size(); i++){
		outputFile <<"\"[DEBUG]index " << i << " " <<unique_motifs[i].name<<"\n";
		consFile <<"\"[DEBUG]index " << i << " " <<unique_motifs[i].name<<"\n";
		
		for(int j=0; j<WINDOW_SIZE; j++){
			outputFile <<j <<" "<<unique_motifs[i].signal[j] << endl;
			consFile <<j <<" "<<unique_motifs[i].con_level[j] << endl;
		}
		motifFile << i << " : " << unique_motifs[i].name << "(" << unique_motifs[i].count << ")"<<endl;
	} 
	
}

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