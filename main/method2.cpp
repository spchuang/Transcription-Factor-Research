//this code matches footprints to predefined template profiles (a few ways to generate this)
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <limits.h>
#include <time.h>
#include <cmath>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "readfp.cpp"
#include "struct.h"
#include "template_generate.cpp"

using namespace std;



void fpMatch(vector<centroid> &C, vector<fpSignalFrame>* f, string cellType);
void print_frame(float Sig[]);

int main(){
	struct stat st = {0};
	

	string cellTypes[] = {"AG10803", "AoAF", "CD20+", "GM06990", "GM12865","H7-hESC","HAEpiC","HA-h","HCF","HCM","HCPEpiC","HEEpiC","HepG2","HFF","HIPEpiC","HMF","HMVEC-dBl-Ad","HPAF","HPdLF","HPF","HRCEpiC","HSMM","HVMF","K562","NB4","NH-A","NHDF-Ad","NHDF-neo","NHLF","SAEC","SKMC","Th1"};
  
  
  string chromosomes[] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                           "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                           "chr21", "chr22","chrX"};
  //23
  int maxChr = 23;
  //33
  int maxCellNum = 33;
  cout <<"\n------------------------------------------------------------\n\n" ;

  cout <<"[DEBUG]Create template signals"<<endl;                  
  vector<centroid> C;
  vector<fpSignalFrame>* fpF = new vector<fpSignalFrame>;
  //generateTemplate(C, 1, 6, 21);
  //generateTemplate(C, 1, 8, 15);
  generateTemplate(C, 1, 9 , 13);

  reduceVector(C,3);
  //print the centroid
  /*for(int i=0; i<C.size(); i++){
     cout <<"\"[DEBUG]index " << i <<" centroid... \"\n";
     for(int j=0; j<WINDOW_SIZE; j++){
        cout <<j <<" "<<C[i].signal[j] << endl;
      }

   }*/
  	//creat tmp folder
  	if (stat("../tmp", &st) == -1) {
		if(mkdir("../tmp", 0700) == -1){
			cout <<"[DEBUG]Error creating tmp folder" <<endl;
		}
		
	}
  	//create directory if it doesn't exists
  	
	if (stat("../tmp/template_assign_full", &st) == -1) {
		if(mkdir("../tmp/template_assign_full", 0700) == -1){
			cout <<"[DEBUG]Error creating assignment folder" <<endl;
		}
		
	}
	
  for(int i=0; i<maxCellNum; i++){
	 cout <<"\n------------------------------------------------------------\n\n" ;
  	cout <<"[DEBUG]Analyzing cell type: " << cellTypes[i] <<endl;
  	for(int j=22; j<maxChr; j++){
		getFootPrint(cellTypes[i], chromosomes[j], (*fpF), 0, true);
	
	}
	fpMatch(C, fpF, cellTypes[i]);
  
  }
  
  /*for(int i=0; i<fpF->size(); i++){
    cout << "foot print start: " << (*fpF)[i].fpStart <<endl;
    for(int j=0; j<FRAME_SIZE; j++){
      cout <<(*fpF)[i].startSeq+j  <<" "<<(*fpF)[i].signal[j] << "  ,  " << (*fpF)[i].con_level[j] << endl;
    }
    
  
  }*/
  
  
 
}

void print_frame(float Sig[])
{
   for(int i=0; i<WINDOW_SIZE; i++){
      cout <<" "<<(Sig[i]) ;
   }
   cout <<endl;
}




//predefine the fixed structure of the centroids
//2 peaks


void fpMatch(vector<centroid> &C, vector<fpSignalFrame>* f, string cellType){
 
  int fpSize = f->size();
  //match each footprint to a template structure

  cout <<"[DEBUG]Total footprints " <<fpSize <<endl;
  cout <<"[DEBUG]Assigning footprints..." <<endl;
  for(int i=0; i<fpSize; i++){
    //cout << "counting for " << i << " / "<< fpSize <<endl;
    int shiftTo = (*f)[i].fpStartIndex;
    bool doFlip = (*f)[i].flip;
    float max = FLT_MIN;
    for(int j=0; j< C.size(); j++){
      //flip
      bool testFlip = false;    
    
      //OFFSET range
      int si = (*f)[i].fpStartIndex - MAX_OFFSET;
      if(si<0) 
        si=0;
      int ei = (*f)[i].fpStartIndex + MAX_OFFSET;
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
          doFlip = testFlip;
        }
      }
      
      //try flip?
      if(allowFlip){
      
       testFlip = true;
        //OFFSET range
      
       int flipped_signal[FRAME_SIZE];       //holding the signal level for each basepair in the frame
       //get flipped signal from (*f)[i]
       int zz=FRAME_SIZE-1;
         for(int z=0; z<FRAME_SIZE; z++){
          flipped_signal[z] = (*f)[i].signal[zz];
          zz--;
        }
       si = (FRAME_SIZE - (*f)[i].fpStartIndex - (*f)[i].fpLength) - MAX_OFFSET;
       if(si<0) 
          si=0; 
        ei = (FRAME_SIZE - (*f)[i].fpStartIndex - (*f)[i].fpLength) + MAX_OFFSET;
        if(ei>= (FRAME_SIZE-WINDOW_SIZE))
          ei = FRAME_SIZE-WINDOW_SIZE-1;
    
       for(int offset= si ; offset <= ei; offset++){
          float d = correlation(flipped_signal, offset, C[j].signal);

          //if distance is smaller than prev min
          if(d > max){
           max = d;
           //change the offset
           shiftTo = offset;
           //assign this point to the centroid
            (*f)[i].clusterAssigned = j;
           doFlip = testFlip;
         }
        }
      }
           
  
  
    }
    if(max <=0){
      (*f)[i].clusterAssigned = -1;
    
    }
    (*f)[i].fpStartIndex = shiftTo;
    (*f)[i].flip         = doFlip;
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
    int fpIndex      = (*f)[i].fpStartIndex;
    ptCount[clusterIndex] = ptCount[clusterIndex]+1;
    int flip_fpsignal[FRAME_SIZE];
    float flip_consSignal[FRAME_SIZE];
    
    if((*f)[i].flip){
      //get flipped signal from (*f)[i]
      int zz=FRAME_SIZE-1;
      for(int z=0; z<FRAME_SIZE; z++){
        flip_fpsignal[z] = (*f)[i].signal[zz];
        flip_consSignal[z] = (*f)[i].con_level[zz];
        zz--;
      }
    }
    
    for(int j=0; j<WINDOW_SIZE; j++){
      if((*f)[i].flip){
        temp_C[clusterIndex].signal[j] += flip_fpsignal[j+fpIndex];
        temp_C[clusterIndex].cons_signal[j] += flip_consSignal[j+fpIndex];
      }else{
        temp_C[clusterIndex].signal[j] += (*f)[i].signal[j+fpIndex];
        temp_C[clusterIndex].cons_signal[j] += (*f)[i].con_level[j+fpIndex];
      }
    }
  
  }
  //remove centroid with empty assignments..
  cout << "[DEBUG]pre size: " << temp_C.size() <<endl;
  for(int i=0; i<temp_C.size(); i++){
    if(ptCount[i] == 0){
      //remove it
      temp_C.erase(temp_C.begin()+i);
      ptCount.erase(ptCount.begin()+i);
      C.erase(C.begin()+i);
      i--;
    }
  }
  cout << "[DEUBG]after size: " << temp_C.size() <<endl;
  //average the signals and reposition the centroid
  for(int i=0; i<temp_C.size(); i++){
    //cout <<"index " << i <<" point count: " << ptCount[i]<<endl;
    for(int j=0; j<WINDOW_SIZE; j++){
      temp_C[i].signal[j] = (float)(temp_C[i].signal[j] / ptCount[i]);
      temp_C[i].cons_signal[j] = (float)(temp_C[i].cons_signal[j] / ptCount[i]);
    }
  }
  
  //calculate confidence interval for the conservation level
  //calculate sum of (Xi - X bar)^2
 
    

  //print the signals out to data
  cout <<"[DEBUG]output result" <<endl;
  struct stat st = {0};
  //create directory if it doesn't exists
  string dirname = "../tmp/template_assign_full/"+cellType;
  if (stat(dirname.c_str(), &st) == -1) {
    mkdir(dirname.c_str(), 0700);
  }
  string filename = dirname+"/fpsig";
  ofstream outputFile(filename.c_str(), ios_base::trunc );
  filename = dirname+"/consSig";
  ofstream consFile(filename.c_str(), ios_base::trunc );
  filename = dirname+"/templateSignal";
  ofstream tempFile(filename.c_str(), ios_base::trunc );
  filename = dirname+"/"+cellType+".count";
  ofstream countFile(filename.c_str(), ios_base::trunc );
  cout.precision(10);
  for(int i=0; i<temp_C.size(); i++){
    outputFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
    tempFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
    consFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";

    for(int j=0; j<WINDOW_SIZE; j++){
      outputFile <<j <<" "<<temp_C[i].signal[j] << endl;
      consFile <<j <<" "<<temp_C[i].cons_signal[j] << endl;
      tempFile <<j <<" "<<C[i].signal[j] << endl;
    }
  } 
  //print centroid count
  cout.precision(2);
  for(int i=0; i<ptCount.size(); i++){
  	
    countFile <<i<<": " << 100.0*(float)ptCount[i]/(float)fpSize<<endl;
  }
/*
  //calculate anticorrelation
  ofstream anticFile("tmp/template_assign/anticorrelation_lv.txt", ios_base::trunc );
  anticFile <<"correlation level (red against blue): " << endl;
  for(int i=0; i<temp_C.size(); i++){
    float d = correlation(temp_C[i].signal, 0, temp_C[i].cons_signal);
    anticFile << i <<" : " << d << endl;
  } 
       */
   C = temp_C;
}
