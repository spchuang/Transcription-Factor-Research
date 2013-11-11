#ifndef DATA_UTILITY_CPP           
#define DATA_UTILITY_CPP 

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <list>
#include <limits.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <float.h>

using namespace std;

const int WINDOW_SIZE = 30;
const int FRAME_SIZE  = 60;
      
struct segment{
   float con_level[FRAME_SIZE];
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
   int segIndex;
};

/*
	Data Master aka the base (abstract) class to handle reading of Dnase signals and conservation data
*/
class BaseDataMaster
{
  public:                    
    BaseDataMaster();     
    ~BaseDataMaster();                 
	void setReadSignal(bool read);
	void setReadCons(bool read);
	void setReadCells(string* cs, int n);
	void setReadChromosomes(string* chrs, int n);
	void setDataDir(string dir);
	void printProperties();
	virtual void startReadingData() = 0; 	//all children classes have to define this


 protected:         
 	
 	void readSignalData(string cell, string chr, int startIndex);
 	void readConsData(string chr, int startIndex);         
  	
  	bool readSignal;
  	bool readCons;
  	string dataDir;
 
 	vector<segment> segs;
    vector<string> cells;
    vector<string> chromosomes;
};




BaseDataMaster::BaseDataMaster(){
	readSignal = false;
	readCons   = false;
	dataDir	   = "../data";
	segs.clear();
	cells.clear();
	chromosomes.clear();
}
BaseDataMaster::~BaseDataMaster(){
	cout << "master is dead. " <<endl;
}

void BaseDataMaster::setReadSignal(bool read){
	readSignal = read;
}
void BaseDataMaster::setReadCons(bool read){
	readCons = read;
	
}
void BaseDataMaster::setReadCells(string* cs, int n){
	for(int i=0; i<n; i++){
		cells.push_back(cs[i]);
	}
	
}
void BaseDataMaster::setReadChromosomes(string* chrs, int n){
	for(int i=0; i<n; i++){
		chromosomes.push_back(chrs[i]);
	}
}

void BaseDataMaster::setDataDir(string dir){
	dataDir = dir;
}

void BaseDataMaster::printProperties(){
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




void BaseDataMaster::readSignalData(string cell, string chr, int startIndex){
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
    it_start = segs.begin()+startIndex;
    int size = segs.size();
	int totalToBeProcessed = size - startIndex;
	//read through the signal data, if the window contains footprint sequence
	//first iteration, find the start of sequence
	string sFile = dataDir+"/signals/"+cell+"_split/"+chr+"."+cell;
	ifstream singalinfile(sFile.c_str());
	
	while(singalinfile >> c >> s >> e >> signal){
		//cout <<"SIGNAL DATA: "<< c << " " << s << "  " <<signal;
		if(!first){
            first = true;
            START = s+1;
         }else{
            //skip useless zero signals...
            int firstsegseq = (it_start->frameStart)+(it_start->length);
            if(START < firstsegseq && firstsegseq < s){
               START = firstsegseq;
            }else if(START<firstsegseq && s < firstsegseq){
               START = s+1;
            }
            //insert zero signal
            while(START < s){
            	
               //iterate through the temp footprint vector and insert the signalFile  
               //cout << endl<<"Zero DATA: "<<c << " " << START << "  0" ;
               if( (it_start->frameStart + it_start->length) >s)
               	break;
               for(it=it_start; it<segs.end(); it++){
               	  //cout<<endl<<(it->frameStart+it->length);
                  //break if the seq is bigger
                  if( (it->frameStart+it->length) > START)
                     break; 
					
                  //if it's the sequence to enter
                  if( (it->frameStart+ it->length) == START){
                  	
                     //cout<<endl<<it->frameStart << " : " <<it->length<<" ("<<(it->frameStart+it->length)<<"): 0";
                     it->signal[it->length] = 0; //ADD THE EMPTY SIGNAL
                     it->length++;
                     //if the point is full, push it to vector f 
                     if(it->length >= FRAME_SIZE){
                     	//cout<<endl<<"break!"<<endl;
                        it->length = 0;
                        totalToBeProcessed--;
                        //cout << totalToBeProcessed<<endl;
                        it_start++;
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
         for(it=it_start; it<segs.end(); it++){
            //break if the seq is bigger
            if( (it->frameStart+it->length) > s)
               break; 
            //if it's the sequence to enter
            if( (it->frameStart+ it->length) == s){
            
               //cout<<endl<<it->frameStart << " : " <<it->length<<" ("<<(it->frameStart+it->length)<<"): " <<signal;
               it->signal[it->length] = signal;
               it->length++;
               //if the point is full, push it to vector f 
               if(it->length >= FRAME_SIZE){
               //cout<<endl<<"break!"<<endl;
                 it->length = 0;
                 totalToBeProcessed--;
                 //cout << totalToBeProcessed<<endl;
                 it_start++;
                 if(totalToBeProcessed == 0)
                 	return;
          
               }
            }
            
         }     
		 //cin.get(); 
	}
	
}

void BaseDataMaster::readConsData(string chr, int startIndex){
	/*
		*conservation data is not cell specific
		conservation data folder structure:
		cons_data/
			/chr1.phyloP46way.placental
			/chr2.phyloP46way.placental
			...
			...
	
		conservation data format:
			fixedStep chrom=chr1 start=10918 step=1
			0.064
			0.056
			0.064
			...
			so the data is (10918:0.064), (10919:0.056)...
		TODO: add offset for the 5 special cell types?
	*/
	string junk, start_text, value;
    int start;
    int size = segs.size();
	int totalToBeProcessed = size - startIndex;

    
	string cons_file = dataDir+"/cons_data/"+chr+".phyloP46way.placental.wigFix";
	ifstream Cinfile(cons_file.c_str());
	
    Cinfile >> junk>> junk >> start_text >> junk;
    start = atoi(start_text.substr(6).c_str());
    vector<segment>::iterator it;
    vector<segment>::iterator it_start;
    it_start = segs.begin()+startIndex;
    
    while(Cinfile >> value){
      //cout <<"start: " << start << ", cons: " << value <<endl;
      if(value == "fixedStep"){
        //cout <<"hit new start: " << start << endl;
        
        Cinfile >> junk >> start_text >> junk;
        //cout <<"test2"<<endl;
        int new_start = atoi(start_text.substr(6).c_str());

        // <<"new start is ... " << new_start <<endl;
        while(start<new_start){
        	//don't waste time
		  if(it_start->frameStart+it_start->length > new_start)
              break;
              
          for(it=it_start; it<segs.end(); it++){
            if(it->frameStart+it->length > start)
              break;
            if(it->frameStart+it->length == start){
              //add 0 for non specified conservation level
              it->con_level[it->length] = 0;
              //cout << "start seq: " << it->startSeq << " , length: " << it->length << ", cons: " << it->con_level[it->length]<<endl;
              it->length++;
              if(it->length >= FRAME_SIZE){
				  it->length = 0;
                  it_start++;
                  totalToBeProcessed--;
                   //cout <<totalToBeProcessed<<endl;
                  if(totalToBeProcessed == 0)
                 	return;
               }
            }
           
          }
          start++;
       } 
        //cout <<"new start: " << start <<endl;
        continue;
      }
      for(it=it_start; it<segs.end(); it++){
         if(it->frameStart+it->length > start)
            break;
         if(it->frameStart+it->length == start){
             
            it->con_level[it->length] = atof(value.c_str());
            //cout << "start seq: " << it->startSeq << " , length: " << it->length << ", cons: " << it->con_level[it->length]<<endl;
            it->length++;
            if(it->length >= FRAME_SIZE){
   
	              it->length = 0;
	              it_start++;
	              totalToBeProcessed--;
	              //cout <<totalToBeProcessed<<endl;
	              if(totalToBeProcessed == 0)
	             	return;
            }
         }
      }
      start++;
      //cin.get();
   }
   Cinfile.close();
   return;


	
}

/*
	Footprint Master aka the ultimate class to handle DNA footprints
*/

class FootprintMaster: public BaseDataMaster
{
   public:
   	  FootprintMaster();
      void startReadingData();
   private: 
   	  int readFootprints(string cell, string chr);
   	  double minFOS;
};

FootprintMaster::FootprintMaster(){
	minFOS	   = 0;
}

void FootprintMaster::startReadingData(){
	for(int i=0; i<cells.size(); i++){
		cout <<"[DEBUG]reading footprints for " << cells[i] <<endl;
		int count=0;
		
		for(int j=0; j<chromosomes.size(); j++){
			cout <<"[DEBUG]reading " << chromosomes[j] <<endl;
			int startIndex = segs.size();
			int chr_count=readFootprints(cells[i], chromosomes[j]);
			cout <<"[DEBUG]"<<chromosomes[j] << " has " << chr_count << " footprints" << endl;
			if(readSignal){
				readSignalData(cells[i], chromosomes[j], startIndex);
			} 
			cout <<"[DEBUG]Done reading signals" << endl;
			if(readCons){
				readConsData(chromosomes[j], startIndex);
			}
			cout <<"[DEBUG]Done reading cons" << endl;
			
			count+=chr_count;
		}
		cout <<"[DEBUG]Read " << count << " footprints for " << cells[i] << endl;
	}
	
	for(int i=segs.size()-10; i<segs.size()-5; i++){
		cout <<"Signals for footprint " << i<<endl;
		for(int f=0; f<FRAME_SIZE; f++){
			cout <<f<<":"<<segs[i].signal[f]<<" , "<<segs[i].con_level[f]<<endl;
		}
	}

	
	cout << "Total number of footprints: " << segs.size() <<endl;
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
		fp data format (BED format):
			
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
			segs.push_back(fpp);
			count ++;
			//cout << count << ": " << s << " " << e<<endl;
		}
	}
	FPinfile.close();
	return count;
		
	
}

/*
	Motif Master aka the ultimate class to handle DNA footprints
*/

class MotifMaster: public BaseDataMaster
{
   public:
   	  MotifMaster();
      void startReadingData();
   private: 
   	  int readMotifData(string chromosome);
   	  double score;
};


MotifMaster::MotifMaster(){
	
}

void MotifMaster::startReadingData(){
	int count=0;
	
	for(int j=0; j<chromosomes.size(); j++){
		cout <<"[DEBUG]reading motif instances for " << chromosomes[j] <<endl;
		int startIndex = segs.size();
		
		int chr_count = readMotifData(chromosomes[j]);
		
		
		count+=chr_count;
	}
	cout <<"[DEBUG]Read " << count << " motifs " <<endl;


	
}

int MotifMaster::readMotifData(string chromosome){
	
	/*
		Motif is not cell specific so each all cells share the same motif. However, with the same motif there are cell specific Dnase signals
		motif data folder structure:
		motif_non_consbased/
			/chr1.motif
			/chr2.motif
			...

		motif data format:
			
			motif_name			 chr  start*    end*	  flip   idk idk
			ERalpha-a_disc4_8mer chr1 2387761   2387774   -      1   0.013341
						
			*location is 1-indexed, end inclusive
		for more: http://www.broadinstitute.org/~pouyak/motif-disc/human/
		returns the number of motifs read for this cell/chr
		Quesiton: what does 1-indexed exactly mean...
	*/
   string motifFIle = dataDir+"/motif_non_consbased/"+chromosome+".motif";

   int s, e;   //start sequence, end sequence
   string chr, motif_name, strand; //chromosome number, motif name

   ifstream MTinFile(motifFIle.c_str());
   string useless1, useless2;
   //NOTE: motif instance interval is inclusive
   int count=0;
   while(MTinFile >> motif_name >> chr >> s >> e >> strand >> useless1 >>useless2){

       e=e+1;
       segment mtt;
       mtt.frameStart = s - (FRAME_SIZE-(e-s))/2;
       mtt.length    = 0;
       mtt.segIndex  = -1;
       mtt.segStart  = s;
       mtt.chr       = new char [chromosome.size()+1];
       mtt.segLength = e-s;
       if(strand == "+")
          mtt.flip = false;
       else
          mtt.flip = true;
      
       mtt.startIndex = 15; //window size of 30bp aligned in the center
       strcpy(mtt.chr, chromosome.c_str());
       mtt.name= new char [motif_name.size()+1];
       strcpy(mtt.name, motif_name.c_str());
       segs.push_back(mtt);
       count++;
      // cout << mtt.fpStart << " , " <<mtt.motifName<< " , " <<mtt.fpLength<< " , " <<mtt.chr<< " , " <<strand<<endl;
  
   }
   cout <<"[DEBUG] motif instance size for this chromosome.." << segs.size()<<endl;
   
   MTinFile.close();
   return count;
}


int main(){
	cout << "okay" <<endl;
	string cellTypes[] = {"AG10803", "AoAF", "CD20+", "GM06990", "GM12865","H7-hESC","HAEpiC","HA-h","HCF","HCM","HCPEpiC","HEEpiC","HepG2","HFF","HIPEpiC","HMF","HMVEC-dBl-Ad","HPAF","HPdLF","HPF","HRCEpiC","HSMM","HVMF","K562","NB4","NH-A","NHDF-Ad","NHDF-neo","NHLF","SAEC","SKMC","Th1"};
	string chromosomes[] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                           "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                           "chr21", "chr22","chrX"};
	/*FootprintMaster test;
	test.setReadSignal(true);
	test.setReadCons(true);
	test.setReadCells(cellTypes, 1);
	test.setReadChromosomes(chromosomes, 1);
	test.printProperties();
	test.startReadingData();*/
	MotifMaster test2;
	test2.setReadSignal(true);
	test2.setReadCons(true);
	test2.setReadCells(cellTypes, 1);
	test2.setReadChromosomes(chromosomes, 1);
	test2.printProperties();
	test2.startReadingData();
}

#endif 


