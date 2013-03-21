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



void fpMatch(vector<centroid> &C, vector<fpSignalFrame>* f, bool getCons);
void print_frame(float Sig[]);

int main(){

  //parameters: window size, frame size, chromosomes, fos lower bound, 
struct stat st = {0};
  //create directory if it doesn't exists


  vector<centroid> C;
  generateTemplate(C, 1, 7, 17);
  //generateCentroids(C, 0);
  reduceVector(C,100);
  //vector<centroid> C2;
  //generateCentroids(C2, 0);
  /*for(int i=0; i<C.size(); i++){
    print_frame(C[i].signal);
 
  }
  exit(1);*/
  cout <<"template size: " << C.size()<<endl;
  vector<fpSignalFrame>* fpF = new vector<fpSignalFrame>;
  string chromosomes[] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                           "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                           "chr21", "chr22","chrX"};
  //cout <<"[DEBUG]start getting footprint data"<<endl;
  /*for(int i=0; i<C.size(); i++){
         cout <<"\"[DEBUG]index " << i <<" centroid... \"\n";
         for(int j=0; j<WINDOW_SIZE; j++){
            cout <<j <<" "<<C[i].signal[j] << endl;
          }

       } */
  //23
  int max = 23;
  for(int i=0; i<max; i++){
     getFootPrint("K562", "K562Sig_filter", chromosomes[i], (*fpF), 0, true);

  }
  /*for(int i=0; i<fpF->size(); i++){
    cout << "foot print start: " << (*fpF)[i].fpStart <<endl;
    for(int j=0; j<FRAME_SIZE; j++){
      cout <<(*fpF)[i].startSeq+j  <<" "<<(*fpF)[i].signal[j] << "  ,  " << (*fpF)[i].con_level[j] << endl;
    }
    
  
  }*/
  //exit(1);
  fpMatch(C, fpF, true);
  
  /*
  ofstream outputFile("tmp/test_data_cluster_new_result_allChr_factor");
   for(int i=0; i<C.size(); i++){
      outputFile << "[START] Centroid : " << i <<endl;
      for(int j=0; j< fpF->size(); j++){
      
            if((*fpF)[j].clusterAssigned == i){
               
               outputFile <<"-----------------------" <<endl;
               outputFile << "window start at: " << (*fpF)[j].startSeq <<endl;
               outputFile << "fpStartIndex: " << (*fpF)[j].fpStartIndex <<endl;
               outputFile << "contain footPrint:";

               outputFile << "footprint start: " << (*fpF)[j].fpStart <<endl;
               outputFile << "footprint length: " << (*fpF)[j].fpLength <<endl;
               for(int z=0; z<FRAME_SIZE; z++){
    
                  if(z>=(*fpF)[j].fpStartIndex && z<(*fpF)[j].fpStartIndex+WINDOW_SIZE)
                     outputFile <<"* (" << z-(*fpF)[j].fpStartIndex<<")";
                  if(z<10)
                     outputFile <<" ";
                  outputFile << (*fpF)[j].startSeq+z <<" : " << (*fpF)[j].signal[z] ;
                  if(z>=(*fpF)[j].fpStart-(*fpF)[j].startSeq && z<(*fpF)[j].fpStart+(*fpF)[j].fpLength-(*fpF)[j].startSeq)
                     outputFile<<" @ |";
                  else 
                     outputFile<<"   |";
                     
                  //print a pseudo graph
                  if(z>=(*fpF)[j].fpStartIndex && z<(*fpF)[j].fpStartIndex+WINDOW_SIZE)
                  {
                     for(int h=0; h<(*fpF)[j].signal[z]; h++){
                        outputFile <<"-" ; 
                     }
                  }
                  outputFile <<endl;
               }
               outputFile <<"-----------------------" <<endl<<endl;
            }
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


void fpMatch(vector<centroid> &C, vector<fpSignalFrame>* f, bool getCons ){
   cout <<"[DEBUG]Assignment fps " <<endl;
  int fpSize = f->size();
  //match each footprint to a template structure

  cout <<fpSize <<endl;
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
  cout << "pre size: " << temp_C.size() <<endl;
  for(int i=0; i<temp_C.size(); i++){
    if(ptCount[i] == 0){
      //remove it
      temp_C.erase(temp_C.begin()+i);
      ptCount.erase(ptCount.begin()+i);
      C.erase(C.begin()+i);
      i--;
    }
  }
  cout << "after size: " << temp_C.size() <<endl;
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
 /* for(int i=0; i<fpSize; i++){
    int clusterIndex = (*f)[i].clusterAssigned;
    int fpIndex      = (*f)[i].fpStartIndex;
    float flip_consSignal[FRAME_SIZE];
    
    if((*f)[i].flip){
      //get flipped signal from (*f)[i]
      int zz=FRAME_SIZE-1;
      for(int z=0; z<FRAME_SIZE; z++){
        
        flip_consSignal[z] = (*f)[i].con_level[zz];
        zz--;
      }
    }
    for(int j=0; j<WINDOW_SIZE; j++){
      if((*f)[i].flip){
        temp_C[clusterIndex].cons_CI[j] += (float)(flip_consSignal[j+fpIndex]-(float)temp_C[clusterIndex].cons_signal[j])*
                                           (float)(flip_consSignal[j+fpIndex]-(float)temp_C[clusterIndex].cons_signal[j]);
      }else{
        temp_C[clusterIndex].cons_CI[j] += (float)((*f)[i].con_level[j+fpIndex]-(float)temp_C[clusterIndex].cons_signal[j])*
                                           (float)((*f)[i].con_level[j+fpIndex]-(float)temp_C[clusterIndex].cons_signal[j]);
      }
    }
    
  }*/
  //finish the CI: caluclaute division by (N-1) and square root
  /*for(int i=0; i<temp_C.size(); i++){
    //cout <<"index " << i <<" point count: " << ptCount[i]<<endl;
    for(int j=0; j<WINDOW_SIZE; j++){
      temp_C[i].cons_CI[j] =  sqrt(temp_C[i].cons_CI[j]/(ptCount[i]-1));
    }
  }*/
    

  //print the signals out to data
  cout <<"[DEBUG]output result" <<endl;
  struct stat st = {0};
  //create directory if it doesn't exists
  if (stat("tmp/test_template", &st) == -1) {
    mkdir("tmp/test_template", 0700);
  }
  ofstream outputFile("tmp/test_template/testing_flip_fpsig", ios_base::trunc );
  ofstream consFile("tmp/test_template/testing_flip_consSig", ios_base::trunc );
  ofstream tempFile("tmp/test_template/testing_filp_temp", ios_base::trunc );
  ofstream anticFile("tmp/test_template/anticorrelation_lv.txt", ios_base::trunc );
  ofstream countFile("tmp/test_template/count.txt", ios_base::trunc );
  //ofstream CIFile("tmp/test_template/cons_CI.txt", ios_base::trunc );
  for(int i=0; i<temp_C.size(); i++){
    outputFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
    tempFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
    //CIFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";

    if(getCons){
      consFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
              
    }
    for(int j=0; j<WINDOW_SIZE; j++){
      outputFile <<j <<" "<<temp_C[i].signal[j] << endl;
      
      //CI range
      //CIFile <<j <<" "<<temp_C[i].cons_signal[j]-1.96*temp_C[i].cons_CI[j] 
      //       << " ~ " <<temp_C[i].cons_signal[j]+1.96*temp_C[i].cons_CI[j] << endl;
      
      if(getCons){
        consFile <<j <<" "<<temp_C[i].cons_signal[j] << endl;
              
      }
      tempFile <<j <<" "<<C[i].signal[j] << endl;
    }
  } 
  //print centroid count
  for(int i=0; i<ptCount.size(); i++){
    countFile <<i<<": " << ptCount[i]<<endl;
  }
  

  
  //calculate anticorrelation
  anticFile <<"correlation level (red against blue): " << endl;
  for(int i=0; i<temp_C.size(); i++){
    float d = correlation(temp_C[i].signal, 0, temp_C[i].cons_signal);
    anticFile << i <<" : " << d << endl;
  } 
       
       /*
       cout <<"conservational data for each cluster: " <<endl<<endl;
      //print conservational average for each cluster
      for(int i=0; i<C.size(); i++){
         centroid s;
         int assigned = 0;
         for(int j=0 ; j<WINDOW_SIZE; j++){
            s.signal[j] = 0;
         }
         
         for(int j=0; j< domainSize; j++){
            if((*f)[j].clusterAssigned == i){
               assigned ++;
               int fpIndex = (*f)[j].fpStartIndex;
               for(int z=0; z<WINDOW_SIZE; z++){
                  s.signal[z] += (*f)[j].con_level[z+fpIndex];
               }   
            }
         }
         for(int j=0; j<WINDOW_SIZE; j++){
            s.signal[j] = s.signal[j]/(float)assigned;
         
         }
         
         //print out this data
         cout <<"\"index " << i <<" centroid...conser \"\n";
         print_frame(s.signal);
         cout <<"\n\n\n";
      }*/
       
        C = temp_C;
}
