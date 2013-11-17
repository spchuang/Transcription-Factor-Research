#include <limits.h>
#include <cstring>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <float.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>
#include <sstream>
#include "readfp.cpp"

#include "template_generate.cpp"
/******************************************************/
using namespace std;
int total_motif=0;
int real_motif = 0;

//footprint data and the signals
//This represents a single point in the domain that will be used for clustering
struct MotifSignalFrame{
   float con_level[FRAME_SIZE];
   int signal[FRAME_SIZE];       //holding the signal level for each basepair in the frame
   int mtStart;
   int startSeq;                 //the start DNA location for the frame (including the actual footprint sequence)
   int length;                   //used for grabbing signal data
   char* chr;                    //the chromosome
   char* motifName;              //name of the motif
   int mtLength;                 //length of the motif (from the fp data)
   int mtStartIndex;             //the offset for the fixed window in the frame
   int clusterAssigned;          //the cluster this point is assigned to
   bool flip;
   int motifIndex;
   
   //also contain fp information that it contains?  actualy, no it doesn't matter
};


void getMotifSignal(string cell, string chromosome, vector<MotifSignalFrame> &f,double score);
void getMotifGraph(vector<MotifSignalFrame>* m, string cellType);


bool sortByCount (centroid i,centroid j) { return (i.count>j.count); }
void mtMatch(vector<centroid> C, vector<MotifSignalFrame> *f, string cellType);
int main(){
	vector<centroid> C;
	
	string templateFile = "templateSignal";
	ifstream readTemplateFile(templateFile.c_str());
	readTemplateFromFile(C, readTemplateFile);
	
	//generateTemplate(C, 1, 8, 15);
  //generateTemplate(C, 1, 9 , 13);
	
  //reduceVector(C,100);
  
    vector<MotifSignalFrame>* mt = new vector<MotifSignalFrame>;
  string chromosomes[] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                           "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                           "chr21", "chr22","chrX"};
  //string cellTypes[] = {"AG10803", "AoAF", "CD20+", "GM06990", "GM12865","H7-hESC","HAEpiC","HA-h","HCF","HCM","HCPEpiC","HEEpiC","HepG2","HFF","HIPEpiC","HMF","HMVEC-dBl-Ad","HPAF","HPdLF","HPF","HRCEpiC","HSMM","HVMF","K562","NB4","NH-A","NHDF-Ad","NHDF-neo","NHLF","SAEC","SKMC","Th1"};
//we have AG10803, AoAF, CD20+, GM12865, HCF, HCM, HMVEC, HPAF, NG-A, NHDF-Ad

//the 5 special ones: K562, HepG2, GM06990, Sknshra, TH1
  /*string cellTypes[] = {"AG10803", "AoAF", "CD20+", "K562","HepG2","GM06990", "GM12865", "H7-hESC","HAEpiC","HA-h","HCF","HCM","HCPEpiC","HEEpiC","HFF","HIPEpiC","HMF","HMVEC-dBl-Ad","HPAF","HPdLF","HPF","HRCEpiC","HSMM","HVMF","NB4","NH-A","NHDF-Ad","NHDF-neo","NHLF","SAEC","SKMC","Th1"};*/
  string cellTypes[] = {"H7-hESC","HAEpiC","HA-h","HCPEpiC","HEEpiC","HFF","HIPEpiC","HMF","HPdLF","HPF","HRCEpiC","HSMM","HVMF","NB4","NHDF-neo","NHLF","SAEC","SKMC"};
  
  int max = 23;
  int maxCellNum = 28;
  
  //create directory folder
	struct stat st = {0};

  if (stat("../tmp/test_motif_pattern", & st) == -1) {
		if(mkdir("../tmp/test_motif_pattern", 0700) == -1){
			cout <<"[DEBUG]Error creating assignment folder" <<endl;
		}
		
	}
	
  for(int x=16; x<18; x++){
	  for(int i=0; i<max; i++){
	     getMotifSignal(cellTypes[x], chromosomes[i], (*mt), 0);
	
	  }
	  cout <<"post read signal size: " <<(*mt).size()<<endl;
	  
	  //reverse the ones that are flip
	  for(int i=0; i< (*mt).size(); i++){
	    if((*mt)[i].flip){
	      int zz=FRAME_SIZE-1; 
	      for(int z=0; z<FRAME_SIZE/2; z++){
	        int temp_s = (*mt)[i].signal[zz];
	        float temp_con = (*mt)[i].con_level[zz];
	        (*mt)[i].signal[zz] = (*mt)[i].signal[z];
	        (*mt)[i].con_level[zz] =(*mt)[i].con_level[z];
	        (*mt)[i].signal[z] = temp_s;
	        (*mt)[i].con_level[z] =temp_con;
	        
	        zz--;
	      }
	    }
	  }
	  cout <<"TOTAL: " <<total_motif<<endl;
	  cout <<"real: " <<real_motif<<endl;
	  cout <<"Percent: " << real_motif/total_motif<<endl;
      
      //mtMatch(C, mt, cellTypes[x]);
	  ///method 2: average the signal for the same motif instance
	  getMotifGraph(mt, cellTypes[x]);
	  mt->clear();


  }
  
 /* for(int i=0; i<mt->size(); i++){
         cout <<"\"[DEBUG]index " << i <<" centroid... \"\n";
         
         cout <<"start at: " <<(*mt)[i].startSeq<<endl;
         cout <<"\"[DEBUG]motif: " << (*mt)[i].motifName <<"\n";
         cout <<"\"[DEBUG]motif actual start: " << (*mt)[i].mtStart <<"\n";
         for(int j=0; j<WINDOW_SIZE; j++){
            cout <<j+(*mt)[i].mtStartIndex +(*mt)[i].startSeq<<" : "<<(*mt)[i].signal[j + (*mt)[i].mtStartIndex] <<
            " , " << (*mt)[i].con_level[j + (*mt)[i].mtStartIndex]<< endl;
          }
       }
*/

}
void mtMatch(vector<centroid> C, vector<MotifSignalFrame> *f, string cellType){
 
  //get list of motif namespace
  cout<<"incoming size: " <<(*f).size()<<endl;
  //char *motifs[500];
  vector<centroid> mt_list;
  int size=0;
  for(int i=0; i<(*f).size(); i++){
    //cout <<(*mt)[i].motifName<<endl;
    bool exists = false;
    for(int j=0; j<size; j++){
      if(strcmp((*f)[i].motifName, mt_list[j].motifName) == 0){
        exists = true;
      }
    }
    
    //create this many centroids
    if(!exists){
      cout << size <<" : " << (*f)[i].motifName <<endl;
      size_t c=0;
      while((*f)[i].motifName[c] != 0)
        c++;
      centroid s;
      for(int j=0; j< WINDOW_SIZE; j++){
        s.signal[j] = 0;
        s.cons_signal[j] = 0;
      }
      s.motifName = new char [c];
      s.count = 0;
      strcpy(s.motifName, (*f)[i].motifName);
      mt_list.push_back(s);
      size++;
    }
  }
  
  cout <<"number of motif: " << mt_list.size() <<endl;
  
  //increment count for the motifs
  for(int i=0; i<mt_list.size(); i++){
    for(int k=0; k<(*f).size(); k++){
      if(strcmp((*f)[k].motifName, mt_list[i].motifName) == 0){
        mt_list[i].count++;
      }
    }
    
  }
  //sort the motifs
  sort (mt_list.begin(), mt_list.end(), sortByCount); 
  
  
  //assign motif numbers to each motif
  for(int i=0; i<(*f).size(); i++){
  	for(int j=0; j<(mt_list.size()); j++){
  		if(strcmp((*f)[i].motifName, mt_list[j].motifName)==0){
  			(*f)[i].motifIndex = j;
  			break;
  		}
  	}
  }
  
  int fpSize = f->size();
  //match each footprint to a template structure

  cout <<"[DEBUG]Total footprints " <<fpSize <<endl;
  cout <<"[DEBUG]Assigning footprints..." <<endl;
  for(int i=0; i<fpSize; i++){
    //cout << "counting for " << i << " / "<< fpSize <<endl;
    int shiftTo = (*f)[i].mtStartIndex;
    float max = FLT_MIN;
    for(int j=0; j< C.size(); j++){
    
      //OFFSET range
      int si = (*f)[i].mtStartIndex - MAX_OFFSET;
      if(si<0) 
        si=0;
      int ei = (*f)[i].mtStartIndex + MAX_OFFSET;
      if(ei>= (FRAME_SIZE-WINDOW_SIZE))
        ei = FRAME_SIZE-WINDOW_SIZE-1;

      for(int offset= si ; offset <= ei; offset++){
        float d = correlation((*f)[i].signal, offset, C[j].signal);

        //if distance is smaller than prev min
        if(d > max){
          max = d;
          //change the offset
          shiftTo = offset;
          //assign this point to the centroid
          (*f)[i].clusterAssigned = j;
        }
      }
    }
    if(max <=0){
      (*f)[i].clusterAssigned = -1;
      
    }
    (*f)[i].mtStartIndex = shiftTo;
  }
  cout <<"[DEBUG]calculating centroid average " <<endl;
  //make a empty copy of temporary centroids
  vector<centroid> temp_C;
  vector<int> ptCount;
  for(int i=0; i<C.size(); i++){
    centroid s;
    for(int j=0 ; j<WINDOW_SIZE; j++){
      s.signal[j] = 0;
      s.cons_signal[j] = 0;
      //s.cons_CI[j] =0;
    }
    temp_C.push_back(s);
    ptCount.push_back(0);
  }
  cout <<"[DEBUG]average centroids" <<endl;
      
  //step 2. reposition each centroid to the average of all the points assigned to iterate
  for(int i=0; i<fpSize; i++){
    int clusterIndex = (*f)[i].clusterAssigned;
    if(clusterIndex == -1)
      continue;
    int fpIndex      = (*f)[i].mtStartIndex;
    ptCount[clusterIndex] = ptCount[clusterIndex]+1;
        
    for(int j=0; j<WINDOW_SIZE; j++){
        temp_C[clusterIndex].signal[j] += (*f)[i].signal[j+fpIndex];
        temp_C[clusterIndex].cons_signal[j] += (*f)[i].con_level[j+fpIndex];
    }
  }
  
  struct stat st = {0};
  
  //create directory if it doesn't exists
  string dirname = "../tmp/test_motif_pattern/"+cellType;
  cout <<"[DEBUG]output direcotry: " << dirname <<endl;
  if (stat(dirname.c_str(), &st) == -1) {
    mkdir(dirname.c_str(), 0700);
  }
  
  string motifFileName = dirname+"/motifOrder.txt";
  ofstream motifFile(motifFileName.c_str(), ios_base::trunc );
  
  for(int i=0; i<mt_list.size(); i++){
    motifFile << mt_list[i].motifName << "\t" << mt_list[i].count <<endl;
  } 
  
  //two output
  //1. For each template, how many percentage from each motif instance
  //template -> Motif
  //create a folder for template,
  cout <<"------------------------"<<endl;
  cout <<"step 1:  For each template, how many percentage from each motif instance " <<endl;
  vector< vector<int> > template_mt;
  vector<int> t_total;
  
  for(int i=0; i<C.size();i++){
  	vector<int> m;
  	for(int j=0; j<mt_list.size(); j++){
  		m.push_back(0);
  	}
  	template_mt.push_back(m);
  	t_total.push_back(0);
  }
  for(int i=0; i<(*f).size(); i++)
  {
  	(template_mt[(*f)[i].clusterAssigned])[(*f)[i].motifIndex]++;
  	t_total[(*f)[i].clusterAssigned]++;
  }
  string dirname1 = dirname+"/mt_perc_in_temp";
  cout <<"[DEBUG]output direcotry:" + dirname1 <<endl;
  if (stat(dirname1.c_str(), &st) == -1) {
    mkdir(dirname1.c_str(), 0700);
  }
  
  //output
  for(int i=0; i<C.size();i++){
  
     assert( t_total[i] == ptCount[i]);
     ostringstream convert;   // stream used for the conversion

	convert << i;      // insert the textual representation of 'Number' in the characters in the stream

     string filename = dirname1 + "/" + string("template_")+ convert.str()+ string(".txt");
     cout <<"[DEBUG]out put file:" + filename <<endl;
     ofstream outputFile(filename.c_str(), ios_base::trunc );
  	 
     cout.precision(2);
  	 assert(template_mt[i].size() == mt_list.size());
  	 for(int j=0; j<template_mt[i].size(); j++){
  	 	outputFile << 100.0*(float)template_mt[i][j]/(float)t_total[i]<<endl;
  	 }
  }
  
 
  //2. for each motif instance, how many percentage is assigned to each template
  //motif -> template
  cout <<"------------------------"<<endl;
  cout <<"step 2: for each motif instance, how many percentage is assigned to each template " <<endl;
  
  vector< vector<int> >mt_template;
  vector<int> mt_total;
  
  for(int i=0; i<mt_list.size();i++){
     vector<int> t;
     for(int j=0; j<C.size(); j++){
     	t.push_back(0);
     }
     mt_template.push_back(t);
     mt_total.push_back(0);
  }
  
  for(int i=0; i<(*f).size(); i++)
  {
  	(mt_template[(*f)[i].motifIndex])[(*f)[i].clusterAssigned]++;
  	mt_total[(*f)[i].motifIndex]++;
  }
  
  
  string dirname2 = dirname+"/tmp_perc_in_mt";
  cout <<"[DEBUG]output direcotry:" + dirname2 <<endl;
  if (stat(dirname2.c_str(), &st) == -1) {
    mkdir(dirname2.c_str(), 0700);
  }
  
  //output
  for(int i=0; i<mt_list.size();i++){
     cout <<"motif: "<<i <<endl;
     cout << mt_total[i] << " against  " << mt_list[i].count<<endl;
//     assert( mt_total[i] == mt_list[i].count);
	ostringstream convert;   // stream used for the conversion

	convert << i;      // insert the textual representation of 'Number' in the characters in the stream

     string filename = dirname2 + string("/") + convert.str() +"-"+string(mt_list[i].motifName)+string(".txt");
     ofstream outputFile(filename.c_str(), ios_base::trunc );
  	 
     cout.precision(2);
  	 assert(mt_template[i].size() == C.size());
  	 for(int j=0; j<mt_template[i].size(); j++){
  	 	outputFile << 100.0*(float)mt_template[i][j]/(float)mt_list[i].count<<endl;
  	 }
  }
  
}



//read the base level conservational energy data for each cluster points
void getMotifConsLevel(string chr, vector<MotifSignalFrame> &temp_con_f, vector<MotifSignalFrame> &f)
{
   cout <<"[DEBUG]reading conservation data for " << chr <<endl;
   string cons_file = "../data/cons_data/"+chr+".phyloP46way.placental.wigFix";
   ifstream Cinfile(cons_file.c_str());
   string a;
   string ss;
   int start;
   Cinfile >> a>> a >> ss >> a;
   
   start = atoi(ss.substr(6).c_str());
   //the offset
   //start--;
   //cout << "cons start at " << start <<endl;
   string value;

   vector<MotifSignalFrame>::iterator it;
   while(Cinfile >> value){
      //cout <<"start: " << start << ", cons: " << value <<endl;
      if(value == "fixedStep"){
        //cout <<"hit new start: " << start << endl;
        
        Cinfile >> a >> ss >> a;
        int new_start = atoi(ss.substr(6).c_str());
        //the offset again
        //new_start--;
       
        while(start<new_start){

          for(it=temp_con_f.begin(); it<temp_con_f.end(); it++){
            if(it->chr == chr && it->startSeq+it->length > start)
              break;
            if(it->chr == chr && it->length+it->startSeq == start){
              //add 0 for non specified conservation level
              it->con_level[it->length] = 0;
              //cout << "start seq: " << it->startSeq << " , length: " << it->length << ", cons: " << it->con_level[it->length]<<endl;
              it->length++;
              if(it->length >= FRAME_SIZE){
                  //cout <<"push1: " <<f.size()<<endl;
                  f.push_back(*it);
                  temp_con_f.erase(it);
                  it--;
                  if(temp_con_f.size() ==0 ){
                    cout <<"before return: " <<f.size()<<endl;
                     return;
                  }
               }
            }
           
          }
          start++;
       } 
        //cout <<"new start: " << start <<endl;
        continue;
      }else{
      
      for(it=temp_con_f.begin(); it<temp_con_f.end(); it++){
         if(it->chr == chr && it->startSeq+it->length > start)
            break;
         if(it->chr == chr && it->length+it->startSeq == start){
            //if(start > 249160630)
            //  cout <<" save " << atof(value.c_str()) << "  in to "<< it->length+it->startSeq <<endl;
            it->con_level[it->length] = atof(value.c_str());
            //cout << "start seq: " << it->startSeq << " , length: " << it->length << ", cons: " << it->con_level[it->length]<<endl;
            it->length++;
            if(it->length >= FRAME_SIZE){
                  //cout <<"push2: " <<f.size()<<endl;
                  f.push_back(*it);
                  temp_con_f.erase(it);
                  it--;
                  if(temp_con_f.size() ==0 ){
                     cout <<"before return: " <<f.size()<<endl;
                     return;
                  }
               }
         }
      }
      }
      start++;
   }
   Cinfile.close();
   cout <<"reach the end: " <<f.size()<<endl;
   return;
}

//read the footprint file for foorprint sequences
void getMotifSignal(string cell, string chromosome, vector<MotifSignalFrame> &f,double score)
{
   string motifFIle = "../data/motif_non_consbased/"+chromosome+".motif";
   vector<MotifSignalFrame> temp_m;
   int s, e;   //start sequence, end sequence
   string chr, mname, strand; //chromosome number, motif name
   cout <<"[DEBUG] get motif instance for " << chromosome<<endl;

   ifstream MTinFile(motifFIle.c_str());
   string useless1, useless2;
   //NOTE: motif instance interval is inclusive
   int z=0;
   while(MTinFile >> mname >> chr >> s >> e >> strand >> useless1 >>useless2){
      //get a specific cell type     
      if(chr == chromosome){
           e=e+1;
           MotifSignalFrame mtt;
           mtt.startSeq = s - (FRAME_SIZE-(e-s))/2;
           mtt.length   = 0;
           mtt.motifIndex = -1;
           mtt.mtStart  = s;
           mtt.chr      = new char [chromosome.size()+1];
           mtt.mtLength = e-s;
           if(strand == "+")
              mtt.flip = false;
           else
              mtt.flip = true;
           mtt.clusterAssigned = -1;
           mtt.mtStartIndex = 15; //window size of 30bp aligned in the center
           strcpy(mtt.chr, chromosome.c_str());
           mtt.motifName= new char [mname.size()+1];
           strcpy(mtt.motifName, mname.c_str());
           temp_m.push_back(mtt);
           z++;
          // cout << mtt.fpStart << " , " <<mtt.motifName<< " , " <<mtt.fpLength<< " , " <<mtt.chr<< " , " <<strand<<endl;
      }
   }
   cout <<"[DEBUG] motif instance size for this chromosome.." << temp_m.size()<<endl;
   total_motif+=temp_m.size();
   MTinFile.close();
   string FPfile = "../data/chr.footprints/"+chromosome+".footprints";
   //string FPfile = "../data/signals/"+cell+"_split/"+chromosome+"."+cell;
   //temporarily holding the fp data
   vector<fpSignalFrame> temp_f;
   vector<MotifSignalFrame> temp_con_f;
   //variables to read in
   double fos; //footprint occupancy score
   string ctype; //cell type
   cout <<"[DEBUG] get fp for " << chromosome<<", "<<cell<<endl;

   ifstream FPinfile(FPfile.c_str());
   while(FPinfile >> chr >> s >> e >> ctype >> fos){
      
      //get a specific cell type     
      if(ctype == cell && chr == chromosome){
        if(fos >= score){
           fpSignalFrame fpp;
           fpp.startSeq = s - (FRAME_SIZE-(e-s))/2;
           fpp.length   = 0;
           fpp.fpStart  = s;
           fpp.chr      = new char [chromosome.size()+1];
           fpp.fpLength = e-s;
           fpp.clusterAssigned = -1;
           fpp.flip = false;
           fpp.fpStartIndex = 15; //window size of 30bp aligned in the center
           strcpy(fpp.chr, chromosome.c_str());
           temp_f.push_back(fpp);
         }
      }
   }
   cout <<"[DEBUG]fp size for this chromosome.." << temp_f.size()<<endl;
   FPinfile.close();
   //check if the motif instance include any fp?
   vector<MotifSignalFrame> real_m;
   vector<fpSignalFrame>::iterator fpi;
   vector<MotifSignalFrame>::iterator mi;
   fpi=temp_f.begin();
   for(mi=temp_m.begin(); mi<temp_m.end(); mi++){
      //cout << temp_m.size() << " , " <<temp_f.size()<<endl;
      while(fpi<temp_f.end()){
        //remove fp instance the fp is chrStart is smaller than mt chrStart
        //if head and tail of fp is less than head of mt
        if(fpi->fpStart < mi->mtStart && fpi->fpStart+fpi->fpLength <mi->mtStart){
          fpi++;

        }else if(
          (fpi->fpStart>=mi->mtStart)&&(fpi->fpStart<= mi->mtStart+mi->mtLength-1) || 
          (fpi->fpStart+fpi->fpLength-1>=mi->mtStart)&&(fpi->fpStart+fpi->fpLength-1<= mi->mtStart+mi->mtLength-1)
        ){
          // cout <<"fp start, end : " << fpi->fpStart << " , " <<fpi->fpStart+fpi->fpLength<<endl;
          // cout <<"mt start, end : " << mi->mtStart << ", " <<mi->mtStart+mi->mtLength<<endl;
          //motif instance include the footprint interval
          real_m.push_back(*mi);
          break;  
        
        }else{
          //this means the motif does not include any footprints...
          break;
        }        
      }         
      if(fpi>=temp_f.end())
        break;       
   }
   cout <<"[DEBUG]Total motif size that include fp: " <<real_m.size()<<endl;
   real_motif+=real_m.size();
   
   //READ FP signals
   int START,signal;
   bool first = false;

   string sFile = "../data/signals/"+cell+"_split/"+chromosome+"."+cell;
   ifstream singalinfile(sFile.c_str());
   
   /*for(int i=0; i<real_m.size();i++){
      cout << real_m[i].startSeq <<endl;
   
   }*/
  int remove =0;
   vector<MotifSignalFrame>::iterator it;
   vector<MotifSignalFrame>::iterator it_start;
   it_start = real_m.begin();
   //read through the signal data, if the window contains footprint sequence
   //first iteration, find the start of sequence
   while(singalinfile >> chr >> s >> e >> signal){
      if(chr == chromosome){

         if(!first){
            first = true;
            START = s+1;
         }else{
            //skip useless zero signals...
            int firstFPSeq = (it_start->startSeq)+(it_start->length);
            if(START < firstFPSeq && firstFPSeq < s)
               START = firstFPSeq;
            else if(START<firstFPSeq && s < firstFPSeq){
               START = s+1;
                if (DEBUG)
                  cout <<"jumped!" <<endl;
               }
               
               
            //insert zero signal
            while(START < s){
     
               //iterate through the temp footprint vector and insert the signalFile  
               for(it=it_start; it<real_m.end(); it++){
                  //break if the seq is bigger
                  if( (it->startSeq+it->length) > START)
                     break; 

                  //if it's the sequence to enter
                  if( (it->startSeq+ it->length) == START){
                     if (DEBUG)
                        cout << it->length << ", " << START<< " : " << 0 <<endl;
                     it->signal[it->length] = 0; //ADD THE EMPTY SIGNAL
                     it->length++;
                     //if the point is full, push it to vector f 
                     if(it->length >= FRAME_SIZE){
                        it->length = 0;
                         remove++;
                        if(remove%1000 == 0)
                           cout <<remove<<endl;
                
                        temp_con_f.push_back(*it);
                        real_m.erase(it);
                        //it--;
                           
                        it_start++;
                        if(remove >=real_m.size() ){
                           //read in conservational level for each point
                           getMotifConsLevel(chromosome, temp_con_f, f);
                           return;
                        }
                        
                     }
                  }
                  
               }
                START++;            
            }
            START = s+1;
         }
         //cout << START <<endl;
         //iterate through the temp footprint vector and insert the signalFile  
         for(it=it_start; it<real_m.end(); it++){
            //break if the seq is bigger
            if( (it->startSeq+it->length) > s)
               break; 
            
            //if it's the sequence to enter
            if( (it->startSeq+ it->length) == s){
              if (DEBUG)
                cout << it->length << ", " << s<< " : " << signal <<endl;

               //cout <<"[DEBUG]: "<<s<<" : " << signal <<endl;
               it->signal[it->length] = signal; //ADD THE EMPTY SIGNAL
               it->length++;
               //if the point is full, push it to vector f 
               if(it->length >= FRAME_SIZE){
                  it->length = 0;
                  temp_con_f.push_back(*it);
                  //real_m.erase(it);
                  //it--;
                  remove++;
                  if(remove%1000 == 0)
                    cout <<remove<<endl;
                  it_start++;
                  if(remove >=real_m.size() ){
                    //read in conservational level for each point
                    getMotifConsLevel(chromosome, temp_con_f, f);
                    return;
                  }

               }
            }
            
         }         

      }
      
   }
   int k=0;
  /* for(it=it_start; it<real_m.end(); it++){
    //print the real_m 
    cout << it_start->mtStart << " , " << it_start->startSeq <<endl;  
    k++;
   }*/
   cout << k  << " left overs..."<<endl;
   cout <<"reach the end when reading signal: " <<f.size()<<endl;
   getMotifConsLevel(chromosome, temp_con_f, f);
   singalinfile.close();
}
void getMotifGraph(vector<MotifSignalFrame>* mt, string cellType){

  cout<<"incoming size: " <<(*mt).size()<<endl;
  //char *motifs[500];
  vector<centroid> C;
  int size=0;
  for(int i=0; i<(*mt).size(); i++){
    //cout <<(*mt)[i].motifName<<endl;
    bool exists = false;
    for(int j=0; j<size; j++){
      if(strcmp((*mt)[i].motifName, C[j].motifName) == 0){
        exists = true;
      }
    }
    
    //create this many centroids
    if(!exists){
      cout << size <<" : " << (*mt)[i].motifName <<endl;
      size_t c=0;
      while((*mt)[i].motifName[c] != 0)
        c++;
      centroid s;
      for(int j=0; j< WINDOW_SIZE; j++){
        s.signal[j] = 0;
        s.cons_signal[j] = 0;
      }
      s.motifName = new char [c];
      s.count = 0;
      strcpy(s.motifName, (*mt)[i].motifName);
      C.push_back(s);
      size++;
    }
  }
  cout <<"number of motif: " << C.size() <<endl;
  
  cout<<"[DEBUG]aggregate signal and conservation signals"<<endl;
  //add the sig and con to centroids
  for(int i=0; i<C.size(); i++){
    for(int k=0; k<(*mt).size(); k++){
      if(strcmp((*mt)[k].motifName, C[i].motifName) == 0){
        int fpIndex      = (*mt)[k].mtStartIndex;
        for(int j=0; j<WINDOW_SIZE; j++){
          C[i].signal[j] += (*mt)[k].signal[j+fpIndex];
          C[i].cons_signal[j] += (*mt)[k].con_level[j+fpIndex];
        }
        C[i].count++;
      }
    }
    
  }
  cout<<"[DEBUG]average the aggregation" <<endl;
  //average the signals and reposition the centroid
  for(int i=0; i<C.size(); i++){
    for(int j=0; j<WINDOW_SIZE; j++){
      C[i].signal[j] = (float)(C[i].signal[j] / C[i].count);
      C[i].cons_signal[j] = (float)(C[i].cons_signal[j] / C[i].count);
    }
  }
  

  cout <<"[DEBUG]sort"<<endl;
  sort (C.begin(), C.end(), sortByCount); 
  cout <<C.size() << " , " << (*mt).size()<<endl;
  
  //print the signals out to data
  cout <<"[DEBUG]output result" <<endl;
  
  struct stat st = {0};
  
  //create directory if it doesn't exists
  //create directory if it doesn't exists
  string dirname = "../tmp/test_motif_pattern"+cellType;
  cout <<"[DEBUG]output direcotry:" + dirname <<endl;
  if (stat(dirname.c_str(), &st) == -1) {
    mkdir(dirname.c_str(), 0700);
  }
  string filename = dirname+"/fpsig";
  ofstream outputFile(filename.c_str(), ios_base::trunc );
  filename = dirname+"/consSig";
  ofstream consFile(filename.c_str(), ios_base::trunc );
  filename = dirname+"/motifOrder";
  ofstream motifFile(filename.c_str(), ios_base::trunc );
  motifFile << "Total" << (*mt).size() <<endl;
  for(int i=0; i<C.size(); i++){
    outputFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
    consFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
              
    for(int j=0; j<WINDOW_SIZE; j++){
      outputFile <<j <<" "<<C[i].signal[j] << endl;
      consFile <<j <<" "<<C[i].cons_signal[j] << endl;
    }
    motifFile << i << " : " << C[i].motifName << "(" << C[i].count << ")"<<endl;
  } 

    //calculate anticorrelation
  filename = dirname+"/cor-level";
  ofstream anticFile(filename.c_str(), ios_base::trunc );
  
  /*ofstream anticFile("tmp/template_assign/anticorrelation_lv.txt", ios_base::trunc );*/
  anticFile <<"correlation level (red against blue): " << endl;
  for(int i=0; i<C.size(); i++){
    float d = correlation(C[i].signal, 0, C[i].cons_signal);
    anticFile  << d << endl;
  }  
  
}

