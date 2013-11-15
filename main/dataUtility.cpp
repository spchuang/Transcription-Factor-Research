#ifndef DATA_UTILITY_CPP           
#define DATA_UTILITY_CPP 

#include <iostream>
#include <algorithm>    // std::sort
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
const int FRAME_SIZE  = 30;
      
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
	void setReadChromosomes(string* chrs, int start, int end);
	void setDataDir(string dir);
	void printProperties();
	virtual void startReadingData() = 0; 	//all children classes have to define this
	vector<segment>* getSegment();
	void clearSegments();
	int getSegSize();
	void sortSegmentsByStart();
	void readSignalData(string cell, string chr, int startIndex);
	void readConsData(string chr, int startIndex);   
	void saveSegs(string dir);
 protected:         
 	
 	
 	      
  	
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
void BaseDataMaster::setReadChromosomes(string* chrs, int start, int end){
	for(int i=start; i<end; i++){
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

vector<segment>* BaseDataMaster::getSegment(){
	return &segs;
	
}

void BaseDataMaster::clearSegments(){
	segs.clear();
}

int BaseDataMaster::getSegSize(){
	return segs.size();
}

//sorting function
bool compareBySegStart(const segment &a, const segment &b)
{
    return a.segStart < b.segStart;
}

void BaseDataMaster::sortSegmentsByStart(){
	
	sort(segs.begin(), segs.end(), compareBySegStart);
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
	int cccc = 0;
	while(singalinfile >> c >> s >> e >> signal){
				//cout <<"SIGNAL DATA: "<< c << " " << s << "  " <<signal <<endl;
		if(!first){
            first = true;
            START = s+1;
         }else{
            
            //insert zero signal
            while(START < s){
            	
               
               //skip useless zero signals...
	            int firstsegseq = (it_start->frameStart)+(it_start->length);
	            if(START < firstsegseq && firstsegseq < s){
	               START = firstsegseq;
	            }else if(s < firstsegseq){
	               break;
	               
	            }
               	
               for(it=it_start; it<segs.end(); it++){
               //cccc++;
               	int current = (it->frameStart+it->length);
               	  //cout<<endl<<(it->frameStart+it->length);
                  //break if the seq is bigger
                  if( current > START)
                     break; 
					
                  //if it's the sequence to enter
                  if( current == START){
                  	//cout << cccc <<endl;
                    // cout<<endl<<it->name << it->frameStart << " : " <<it->length<<" ("<<(it->frameStart+it->length)<<"): 0";
                     it->signal[it->length] = 0; //ADD THE EMPTY SIGNAL
                     it->length++;
                     //if the point is full, push it to vector f 
                     if(it->length >= FRAME_SIZE){
                     	//cout<<endl<<"break!"<<endl;
                        it->length = 0;
                        totalToBeProcessed--;
                        if(totalToBeProcessed % 1000 ==0) cout << totalToBeProcessed <<endl;

                        //cout << totalToBeProcessed<<endl;
                        ++it_start;
						if(totalToBeProcessed == 0){
							return;
						}
                     }
                    // cin.get(); 
                  }
                  
               }
               START++;   
                        
            }
            START = s+1;
         }
         
         //iterate through the temp footprint vector and insert the signalFile  
         //for(it=temp_f.begin(); it<temp_f.end(); it++){
         for(it=it_start; it<segs.end(); it++){
        // cccc++;
        	int current = (it->frameStart+it->length);
            //break if the seq is bigger
            if( current > s)
               break; 
            //if it's the sequence to enter
            if( current == s){
            
               //cout<<endl<<it->frameStart << " : " <<it->length<<" ("<<(it->frameStart+it->length)<<"): " <<signal;
               it->signal[it->length] = signal;
               it->length++;
               //if the point is full, push it to vector f 
               if(it->length >= FRAME_SIZE){
               //cout<<endl<<"break!"<<endl;
                 it->length = 0;
                 totalToBeProcessed--;
                 //cout << totalToBeProcessed<<endl;
                 if(totalToBeProcessed % 1000 ==0) cout << totalToBeProcessed <<endl;

                 ++it_start;
                 if(totalToBeProcessed == 0){
                 	return;
                 }
          
               }
     
            }
            
         }     
		 
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
              cout << "start seq: " << it->frameStart << " , length: " << it->length << ", cons: " << it->con_level[it->length]<<endl;
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
            cout << "start seq: " << it->frameStart << " , length: " << it->length << ", cons: " << it->con_level[it->length]<<endl;
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
      cin.get();
   }
   Cinfile.close();
   return;


	
}

void BaseDataMaster::saveSegs(string dir ){
	
	

	
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
			//cout <<"[DEBUG]reading " << chromosomes[j] <<endl;
			int startIndex = segs.size();
			int chr_count=readFootprints(cells[i], chromosomes[j]);
			cout <<"[DEBUG]"<<chromosomes[j] << " has " << chr_count << " footprints" << endl;
			if(readSignal){
				readSignalData(cells[i], chromosomes[j], startIndex);
				cout <<"[DEBUG]Done reading signals" << endl;

			} 
			if(readCons){
				readConsData(chromosomes[j], startIndex);
				cout <<"[DEBUG]Done reading cons" << endl;
			}
			
			
			count+=chr_count;
		}
		cout <<"[DEBUG]Read " << count << " footprints for " << cells[i] << endl;
	}
	/*
	for(int i=segs.size()-10; i<segs.size()-5; i++){
		cout <<"Signals for footprint " << i<<endl;
		for(int f=0; f<FRAME_SIZE; f++){
			cout <<f<<":"<<segs[i].signal[f]<<" , "<<segs[i].con_level[f]<<endl;
		}
	}*/

	
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
			fpp.startIndex = (FRAME_SIZE-WINDOW_SIZE)/2; //window size of 30bp aligned in the center
			strcpy(fpp.chr, chr.c_str());
			for(int i=0; i<FRAME_SIZE; i++){
				fpp.signal[i] = 0;
				fpp.con_level[i] = 0;
		    }
			
			
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
   	  void writeMotifData(string cell, string chromosome);
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
		cout <<"[DEBUG]"<<chromosomes[j] << " has " << chr_count << " motifs" << endl;
		//Filter with footprints
		
		if(readSignal){
			readSignalData(cells[0], chromosomes[j], startIndex);
			cout <<"[DEBUG]Done reading signals" << endl;

		} 
		
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
		NOTE: The same motifs have the same segment length (HOLDS TRUE TO ALL MOTIF INSTANCES)
	*/
   string motifFIle = dataDir+"/motif_non_consbased/"+chromosome+".motif";

   int s, e;   //start sequence, end sequence
   string chr, motif_name, strand; //chromosome number, motif name

   ifstream MTinFile(motifFIle.c_str());
   string useless1, useless2;
   //NOTE: motif instance interval is inclusive
   int count=0;
   while(MTinFile >> motif_name >> chr >> s >> e >> strand >> useless1 >>useless2){

       e=e-1;
       s=s-1;
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
      
       mtt.startIndex = (FRAME_SIZE-WINDOW_SIZE)/2; //window size of 30bp aligned in the center
       strcpy(mtt.chr, chromosome.c_str());
       mtt.name= new char [motif_name.size()+1];
       strcpy(mtt.name, motif_name.c_str());
       for(int i=0; i<FRAME_SIZE; i++){
	      mtt.signal[i] = 0;
	      mtt.con_level[i] = 0;
      }
       
       
       segs.push_back(mtt);
       
       count++;
      // cout << mtt.segStart << " , " <<mtt.motifName<< " , " <<mtt.fpLength<< " , " <<mtt.chr<< " , " <<strand<<endl;
      
      //initialize 
      
   }
   
   MTinFile.close();
   return count;
}

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
                           
	MotifMaster motif_master;
	//motif_master.setReadSignal(true);
	FootprintMaster fp_master;
	//fp_master.setReadSignal(true);
	motif_master.setReadCells(cellTypes, 1);
	//foot print reads each chr for ALL cell types
	fp_master.setReadCells(cellTypes, 31);
	
	
	cout << "LETS START THIS SHIT" <<endl;
	for(int chr=0; chr<1; chr++){
		cout <<"for motif chr " << chr<<endl;
		motif_master.printProperties();
		motif_master.setReadChromosomes(chromosomes, chr, chr+1);
		motif_master.startReadingData();
		
		cout <<"for fp chr " << chr<<endl;
		fp_master.setReadChromosomes(chromosomes,chr, chr+1);
		fp_master.startReadingData();
		
		
		//cout <<"sort"<<endl;
		fp_master.sortSegmentsByStart();
		//filter out motifs that are not included in any footprint
		vector<segment>* fp_segs = fp_master.getSegment();
		vector<segment>* mt_segs = motif_master.getSegment();
		int mt_size = mt_segs->size();
		int fp_size = fp_segs->size();
		int j=0;
		vector<segment> new_mt;
		for(int i=0; i<mt_size; i++){
			bool containsFootprint = false;
			while(j<fp_size){
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
		        	//cout << fp_start<< " : " <<fp_end << " ,  "<< mt_start <<" : " <<mt_end<<endl;
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
		new_mt.clear();
		int new_mt_size = motif_master.getSegSize();
		
		//cout <<"new size: " <<new_mt.size() <<endl;
		//cout <<"fp size: "<<fp_segs->size()<<endl;
		cout <<"new motif size: "<<motif_master.getSegSize()<< " / " <<mt_size<<endl;
		
		for(int cell = 0; cell<31; cell++){
			cout <<"start reading signal for motif at chr " << chromosomes[chr]<< " , " <<cellTypes[cell]<<endl;
			motif_master.readSignalData(cellTypes[cell], chromosomes[chr], 0);
			/*
			string dirname = "TEST/motif_signals/"+cellTypes[cell];
			cout <<dirname<<endl;
			cout <<"[DEBUG]output direcotry:" + dirname <<endl;
			if (stat(dirname.c_str(), &st) == -1) {
				mkdir(dirname.c_str(), 0700);
			}
			string filename = dirname+"/"+chromosomes[chr];
			cout <<filename<<endl;
			ofstream outputFile(filename.c_str(), ios_base::trunc );
			*/
			for(int i=0; i< new_mt_size; i++){
			    //outputFile <<"[DEBUG]index " << i <<" centroid... \n";
			    cout << (*mt_segs)[i].name << " : (" <<(*mt_segs)[i].segStart << " , " <<(*mt_segs)[i].segStart+(*mt_segs)[i].segLength<<endl;
			             
			    for(int j=0; j<WINDOW_SIZE; j++){
			      cout <<(*mt_segs)[i].signal[j] << endl;
			    }
		
			} 
			
		}
		
		
		
		
		fp_master.clearSegments();
		
	}
	
	
	
	return 0;
}

#endif 


