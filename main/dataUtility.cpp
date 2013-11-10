#ifndef DATA_UTILITY_CPP           
#define DATA_UTILITY_CPP 

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <list>
#include <limits.h>
#include <cstring>
#include <cmath>
#include <float.h>

using namespace std;

const int WINDOW_SIZE = 30;
const int FRAME_SIZE  = 60;
      
struct segment{
   float con_data[FRAME_SIZE];
   int 	 signal[FRAME_SIZE];       //holding the signal level for each basepair in the frame
   int   segStart;
   int   frameStart;                 //the start DNA location for the frame (including the actual footprint sequence)
   int   length;                   //used for grabbing signal data
   char* chr;                    //the chromosome
   char* cell;
   char* name;
   int   segLength;                 //length of the footprint (from the fp data)
   int   startIndex;             //the offset for the fixed window in the frame
   bool flip;
};

/*
	Footprint Master aka the ultimate class to handle DNA footprints
*/
class FootprintMaster
{
  public:                    
    FootprintMaster();     
    ~FootprintMaster();                 
	void setReadSignal(bool read);
	void setReadCons(bool read);
	void setReadCells(string* cs, int n);
	void setReadChromosomes(string* chrs, int n);
	void setDataDir(string dir);
	void printProperties();
	void startReadingData();


 private:         
 	int readFootprints(string cell, string chr);
 	void readSignalData(string cell, string chr, int startIndex);
 	void readConsData(string cell, string chr, int startIndex);         
  	
  	bool readSignal;
  	bool readCons;
  	string dataDir;
  	double minFOS;
 	vector<segment> fps;
    vector<string> cells;
    vector<string> chromosomes;
};


FootprintMaster::FootprintMaster(){
	readSignal = false;
	readCons   = false;
	dataDir	   = "../data";
	minFOS	   = 0;
	fps.clear();
	cells.clear();
	chromosomes.clear();
}
FootprintMaster::~FootprintMaster(){
	cout << "master is dead. " <<endl;
}

void FootprintMaster::setReadSignal(bool read){
	readSignal = read;
}
void FootprintMaster::setReadCons(bool read){
	readCons = read;
	
}
void FootprintMaster::setReadCells(string* cs, int n){
	for(int i=0; i<n; i++){
		cells.push_back(cs[i]);
	}
	
}
void FootprintMaster::setReadChromosomes(string* chrs, int n){
	for(int i=0; i<n; i++){
		chromosomes.push_back(chrs[i]);
	}
}

void FootprintMaster::setDataDir(string dir){
	dataDir = dir;
}

void FootprintMaster::printProperties(){
	cout << "read signal: " << readSignal <<endl;
	cout << "read conservation data: " << readCons <<endl;
	cout << "reading data from directory: '" << dataDir <<"'"<<endl;
	cout << "read cells: " << endl;
	for(int i=0; i<cells.size(); i++){
		cout << cells[i] << " ";
	}
	cout << endl;
	cout << "read chromosomes: " << endl;
	for(int i=0; i<chromosomes.size(); i++){
		cout << chromosomes[i] << " ";
	}
	cout << endl;
}

void FootprintMaster::startReadingData(){
	for(int i=0; i<cells.size(); i++){
		cout <<"[DEBUG]reading footprints for " << cells[i] <<endl;
		int count=0;
		
		for(int j=0; j<chromosomes.size(); j++){
			cout <<"[DEBUG]reading " << chromosomes[j] <<endl;
			int startIndex = fps.size();
			int chr_count=readFootprints(cells[i], chromosomes[j]);
			cout <<"[DEBUG]"<<chromosomes[j] << " has " << chr_count << " footprints" << endl;
			if(readSignal){
				readSignalData(cells[i], chromosomes[j], startIndex);
			} 
			if(readCons){
				readConsData(cells[i], chromosomes[j], startIndex);
			}
			
			count+=chr_count;
		}
		cout <<"[DEBUG]Read " << count << " footprints for " << cells[i] << endl;
	}

	
	cout << "Total number of footprints: " << fps.size() <<endl;
	//if(readSignal) readSignalData();
	//if(readCons) readConsData();
}

int FootprintMaster::readFootprints(string cell, string chr){
	/*
		footprint folder structure:
		footprintsData/
			/celltype
				/chr1
				/chr2
				...
			...
		This file structure is created using script/split_fp_data.sh	
		fp data format:
			
			chr		start*	end*	cell		FOS (lower the better)
			chr1	10150	10157	HMVEC-LLy	0.942857
			
			*location is 0-based
		for more: ftp://ftp.ebi.ac.uk/pub/databases/ensembl/encode/supplementary/integration_data_jan2011/byDataType/footprints/README
		returns the number of footprints read for this cell/chr
	*/
	int s, e;   //start sequence, end sequence
    double fos; //footprint occupancy score
    string c, ctype; //chromosome number, cell type
	string FPfile = dataDir + "/footprintsData/"+cell+"/"+chr;
	ifstream FPinfile(FPfile.c_str());
	int count =0;
	while(FPinfile >> c >> s >> e >> ctype >> fos){
     	
		if(fos >= minFOS){
			segment fpp;
			fpp.frameStart = s - (FRAME_SIZE-(e-s))/2;
			fpp.length   = 0;
			fpp.segStart  = s;
			fpp.chr      = new char [chr.size()+1];
			fpp.segLength = e-s;
			fpp.flip = false;
			fpp.startIndex = 15; //window size of 30bp aligned in the center
			strcpy(fpp.chr, chr.c_str());
			fps.push_back(fpp);
			count ++;
			//cout << count << ": " << s << " " << e<<endl;
		}
	}
	FPinfile.close();
	return count;
		
	
}

void FootprintMaster::readSignalData(string cell, string chr, int startIndex){
	/*
		signal folder structure:
		signals/
			/celltype
				/chr1.celltype
				/chr2.celltype
				...
			...
		This file structure is created using script/split_signal.sh	
		signal data format:
			
			chr		start*	end*	signal
			chr10	60823	60824	1
	*/
	
	int s, e, signal, START;
	string c; 
	bool first = false;
	vector<segment>::iterator it_start;
    vector<segment>::iterator it;
    it_start = fps.begin()+startIndex;
    int size = fps.size();
	int totalToBeProcessed = size - startIndex;
	//read through the signal data, if the window contains footprint sequence
	//first iteration, find the start of sequence
	string sFile = "../data/signals/"+cell+"_split/"+chr+"."+cell;
	ifstream singalinfile(sFile.c_str());

	while(singalinfile >> c >> s >> e >> signal){
		
		if(!first){
            first = true;
            START = s+1;
         }else{
            //skip useless zero signals...
            int firstFPSeq = (it_start->frameStart)+(it_start->length);
            if(START < firstFPSeq && firstFPSeq < s){
               START = firstFPSeq;
            }else if(START<firstFPSeq && s < firstFPSeq){
               START = s+1;
            }
            //insert zero signal
            while(START < s){
               //iterate through the temp footprint vector and insert the signalFile  
               //it = temp_f.begin()
               for(it=it_start; it<fps.end(); it++){
               		
                  //break if the seq is bigger
                  if( (it->frameStart+it->length) > START)
                     break; 

                  //if it's the sequence to enter
                  if( (it->frameStart+ it->length) == START){
                 // cout<<it->frameStart << " : " <<it->length<<" : " <<0<<endl;
                     it->signal[it->length] = 0; //ADD THE EMPTY SIGNAL
                     it->length++;
                     //if the point is full, push it to vector f 
                     if(it->length >= FRAME_SIZE){
                        it->length = 0;
                        totalToBeProcessed--;
                        //cout << totalToBeProcessed<<endl;
                        it++;
						if(totalToBeProcessed == 0)
							return;
                     }
                  }
                  
               }
               START++;            
            }
            START = s+1;
         }
         //iterate through the temp footprint vector and insert the signalFile  
         //for(it=temp_f.begin(); it<temp_f.end(); it++){
         for(it=it_start; it<fps.end(); it++){
            //break if the seq is bigger
            if( (it->frameStart+it->length) > s)
               break; 
            
            //if it's the sequence to enter
            if( (it->frameStart+ it->length) == s){
            
              	// cout<<it->frameStart << " : " <<it->length<<" : " <<signal<<endl;
               it->signal[it->length] = signal;
               it->length++;
               //if the point is full, push it to vector f 
               if(it->length >= FRAME_SIZE){
                 it->length = 0;
                 totalToBeProcessed--;
                 //cout << totalToBeProcessed<<endl;
                 it++;
                 if(totalToBeProcessed == 0)
                 	return;
          
               }
            }
            
         }     
   
	}
	
}

void FootprintMaster::readConsData(string cell, string chr, int startIndex){
	
	
}

int main(){
	cout << "okay" <<endl;
	string cellTypes[] = {"AG10803", "AoAF", "CD20+", "GM06990", "GM12865","H7-hESC","HAEpiC","HA-h","HCF","HCM","HCPEpiC","HEEpiC","HepG2","HFF","HIPEpiC","HMF","HMVEC-dBl-Ad","HPAF","HPdLF","HPF","HRCEpiC","HSMM","HVMF","K562","NB4","NH-A","NHDF-Ad","NHDF-neo","NHLF","SAEC","SKMC","Th1"};
	string chromosomes[] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                           "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                           "chr21", "chr22","chrX"};
	FootprintMaster test;
	test.setReadSignal(true);
	test.setReadCells(cellTypes, 1);
	test.setReadChromosomes(chromosomes, 1);
	test.printProperties();
	test.startReadingData();
}

#endif 


