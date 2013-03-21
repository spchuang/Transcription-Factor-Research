#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <list>
#include <limits.h>
#include <cstring>
#include <cstdlib>
#include <time.h>
#include <cmath>
#include <float.h>
#include <sys/stat.h>
#include <unistd.h>
#include "readfp.cpp"
#include <thread>
#include "struct.h"
#define MAX_THREADS 4
using namespace std;
const int CLUSTER_NUM = 80;   
//cluster centroid



//return a centroid given a fp point, and the starting index of the fixed window
centroid getCentroid(vector<fpSignalFrame>* f, int index)
{
    centroid s;
    int startIndex = (*f)[index].fpStartIndex;
    for(int i=0 ; i<WINDOW_SIZE; i++){
        s.signal[i] = (*f)[index].signal[i+startIndex];
        s.cons_signal[i] = (*f)[index].con_level[i+startIndex];
    }
    return s;
}


void print_frame(float Sig[])
{
   for(int i=0; i<WINDOW_SIZE; i++){
      cout <<i <<" "<<(Sig[i]) << endl;
   }
   cout <<endl;
}

void assign_points(vector<centroid> &C, vector<fpSignalFrame>* f, int start, int end){
    for(int i=start; i< end; i++){
            float max = FLT_MIN;
            int shiftTo = (*f)[i].fpStartIndex;
            bool doFlip = (*f)[i].flip;
            for(int j=0; j< C.size(); j++){
                //NOTE: doesn't keep track of which centroid the poitn is assigned to (just the distance to the closest centroid)
                //flip
                bool testFlip = false;
                //OFFSET range
                int si = (*f)[i].fpStartIndex - MAX_OFFSET;
                if(si<0) 
                  si=0;
                int ei = (*f)[i].fpStartIndex + MAX_OFFSET;
                if(ei>= (FRAME_SIZE-WINDOW_SIZE))
                  ei = FRAME_SIZE-WINDOW_SIZE-1;


                for(int offset= si ; offset < ei; offset++){
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
      
                 int flipped_signal[FRAME_SIZE];       //holding the signal level for each basepair in the frame
                 // get flipped signal from (*f)[i]
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
            (*f)[i].fpStartIndex = shiftTo;
            (*f)[i].flip         = doFlip;
    }

}

//apply k-means clustering on all the fp points
void k_mean_correlation(vector<centroid> &C, vector<fpSignalFrame>* f)
{
    cout <<"[DEBUG]Start initializing centroids"<<endl;
    //initialize random variable
    srand (time(NULL) );
    //initializing k clusters using k-mean++ algorithm
    int domainSize = f->size();
    //step 1) pick a random X from the domain
    int s1Index= rand() % domainSize;  
    cout <<"[DEBUG]domainSize: "<<domainSize<<endl;
    cout <<"[DEBUG]push first random centroid: "<<s1Index <<endl;
    C.push_back(getCentroid(f, s1Index));
    
    
    vector<float> probPerPoint;
    while(C.size() < CLUSTER_NUM){
      int sum=0;
      //step 2a. calculate shortest distance between Xi and the nearest centroid
      for(int i=0; i< domainSize; i++){
            int shiftTo = (*f)[i].fpStartIndex;
            probPerPoint.push_back(0);
            int min = INT_MAX;
            for(int j=0; j< C.size(); j++){
                //NOTE: doesn't keep track of which centroid the poitn is assigned to (just the distance to the closest centroid)
                //OFFSET range
                int si = (*f)[i].fpStartIndex - MAX_OFFSET;
                if(si<0) 
                  si=0;
                int ei = (*f)[i].fpStartIndex + MAX_OFFSET;
                if(ei>= (FRAME_SIZE-WINDOW_SIZE))
                  ei = FRAME_SIZE-WINDOW_SIZE-1;


                for(int offset= si ; offset < ei; offset++){
                    int d = distance((*f)[i].signal, offset, C[j].signal);
                    //if distance is smaller than prev min
                    if(d < min){
                        min = d;
                        //change the offset
                        shiftTo = offset;
                        probPerPoint[i] = (float)min;
                    }
                }
            }
            (*f)[i].fpStartIndex = shiftTo;
            sum += (int)probPerPoint[i];
        }
        
        //step 2b. find the probability for each point
        for(int i=0; i<domainSize; i++){
            probPerPoint[i] = probPerPoint[i]/sum;
        }
        
        //step 2c. generate a random decimal and use it to assign a new centroid point
        float r = (float)rand()/(float)RAND_MAX;
        float total_r=0;
        for(int i=0; i<domainSize; i++){
            total_r+=probPerPoint[i];
            if(r<total_r){
                cout <<"push next centroid: "<<i <<endl;
                C.push_back(getCentroid(f, i));
                break;
            }
        }
        probPerPoint.clear();
    }
    //DEBUG
    /*
    for(int i=0; i< C.size(); i++){
      cout <<"\"[DEBUG]index " << i <<" centroid... \"\n";
                  
        for(int j=0; j<WINDOW_SIZE; j++){
          cout <<j <<" "<<C[i].cons_signal[j];
        }
        cout <<endl;
    
    }*/
    
    
    //clustering...K-mean
    cout <<"[DEBUG]start clustering... " <<endl;
    bool converge = false;
    vector<int> ptCount2;
    
    //while the centroids are not converging
    for(int x =0 ; x< MAX_ITERATION; x++){
       converge = true;
       cout <<"[DEBUG]cluster iteration : " << x<<endl;
       
       //step 1. assign each point in the domain to a closest centroids
       
       //try wih thread
       int size = domainSize/MAX_THREADS;
       thread t1(assign_points, C, f, 0, size);       
       thread t2(assign_points, C, f, size, size*3);       
       thread t3(assign_points, C, f, size*2, size*3);       
       thread t4(assign_points, C, f, size*3, domainSize);       
       
       t1.join();
       t2.join();
       t3.join();
       t4.join();
 
       //make a empty copy of temporary centroids
       vector<centroid> temp_C;
       vector<int> ptCount;
       for(int i=0; i<C.size(); i++){
     
          centroid s;
          for(int j=0 ; j<WINDOW_SIZE; j++){
             s.signal[j] = 0;
             s.cons_signal[j] = 0;
          }
          temp_C.push_back(s);
          ptCount.push_back(0);
       }
      
       //step 2. reposition each centroid to the average of all the points assigned to iterate
       for(int i=0; i<domainSize; i++){
          int clusterIndex = (*f)[i].clusterAssigned;
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
       //average the signals and reposition the centroid
       for(int i=0; i<C.size(); i++){
          //cout <<"index " << i <<" point count: " << ptCount[i]<<endl;
          for(int j=0; j<WINDOW_SIZE; j++){
            temp_C[i].signal[j] = (float)(temp_C[i].signal[j] / ptCount[i]);
            temp_C[i].cons_signal[j] = (float)(temp_C[i].cons_signal[j] / ptCount[i]);
             if(C[i].signal[j] != temp_C[i].signal[j]){
                C[i].signal[j]=temp_C[i].signal[j];
                converge = false;            
             }
          }

       }
       temp_C.clear();
       ptCount2 = ptCount; 
       
      /*cout <<"[DEBUG] SIGNALS output" <<endl;
      for(int i=0; i<C.size(); i++){
        cout <<"\"[DEBUG]index " << i <<" centroid... \"\n";
                  
        for(int j=0; j<WINDOW_SIZE; j++){
          cout <<j <<" "<<C[i].signal[j] << endl;
        }
      } 
      cout <<"[DEBUG] cons output" <<endl;
      for(int i=0; i<C.size(); i++){
        cout <<"\"[DEBUG]index " << i <<" centroid... \"\n";
                  
        for(int j=0; j<WINDOW_SIZE; j++){
          cout <<j <<" "<<C[i].cons_signal[j] << endl;
        }
      } 
      cout <<"[DEBUG]count"<<endl;
      for(int i=0; i<C.size(); i++){
        cout << i << " : " << ptCount[i] <<endl;
       
      }*/
      cout <<"[DEBUG]output result" <<endl;
      struct stat st = {0};
      //create directory if it doesn't exists
      if (stat("tmp/test_cluster2", &st) == -1) {
        mkdir("tmp/test_cluster2", 0700);
      }
      ofstream outputFile("tmp/test_cluster2/fpsig", ios_base::trunc );
      ofstream consFile("tmp/test_cluster2/consSig", ios_base::trunc );
      ofstream count("tmp/test_cluster2/count.txt", ios_base::trunc );
      for(int i=0; i<C.size(); i++){
        outputFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
        consFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
                  
        for(int j=0; j<WINDOW_SIZE; j++){
          outputFile <<j <<" "<<C[i].signal[j] << endl;
          consFile <<j <<" "<<C[i].cons_signal[j] << endl;
        }
        count << i << " : " << ptCount2[i] <<endl;
      } 
       
    }
    /*for(int i=0; i<C.size(); i++){
         cout <<"\"[DEBUG]index " << i <<" centroid... \"\n";
         print_frame(C[i].signal);
         cout <<"\n\n\n";
       } */
    //print the signals out to data
    

  cout <<"[DEBUG]output result" <<endl;
  struct stat st = {0};
  //create directory if it doesn't exists
  if (stat("tmp/test_cluster2", &st) == -1) {
    mkdir("tmp/test_cluster2", 0700);
  }
  ofstream outputFile("tmp/test_cluster2/fpsig", ios_base::trunc );
  ofstream consFile("tmp/test_cluster2/consSig", ios_base::trunc );
  ofstream count("tmp/test_cluster2/count.txt", ios_base::trunc );
  for(int i=0; i<C.size(); i++){
    outputFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
    consFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
              
    for(int j=0; j<WINDOW_SIZE; j++){
      outputFile <<j <<" "<<C[i].signal[j] << endl;
      consFile <<j <<" "<<C[i].cons_signal[j] << endl;
    }
    count << i << " : " << ptCount2[i] <<endl;
  } 
  
}

int main()
{
   vector<fpSignalFrame>* fpF = new vector<fpSignalFrame>;
   string append_name = "_clusterSize25";
   vector<centroid> C;
   string chromosomes[] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                           "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                           "chr21", "chr22","chrX"};
   cout <<"[DEBUG]start getting footprint data"<<endl;
   //23
   for(int i=0; i<23; i++){
      getFootPrint("K562", "K562Sig_filter", chromosomes[i], (*fpF), 0.5, true);

   }
   /*
   printf("start: %d\n", (*fpF)[0].fpStart);
   printf("frame start: %d\n",(*fpF)[0].startSeq);
   for(int i=0; i<FRAME_SIZE; i++){
    printf("%d (%d): %d , %f\n", i, i+(*fpF)[0].startSeq,(*fpF)[0].signal[i], (*fpF)[0].con_level[i]);
   
   }*/

   k_mean_correlation(C, fpF);

   delete fpF;
   return 0;

}

