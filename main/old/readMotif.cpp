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
#include "readfp.cpp"
using namespace std;
int total_motif=0;
int real_motif = 0;

bool allowFlip = true;
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
   
   //also contain fp information that it contains?  actualy, no it doesn't matter
};
struct centroid{
   float signal[WINDOW_SIZE];      //signal level for the cluster of fixed window size
   float cons_signal[WINDOW_SIZE];
   float cons_CI[WINDOW_SIZE];
   char* motifName;
   int count;
};
void getMotifSignal(string cell, string signalFile, string chromosome, vector<MotifSignalFrame> &f,double score, bool getCons);
void generateCentroids(vector<centroid> &C, int n);
void fpMatch(vector<centroid> &C, vector<MotifSignalFrame>* f, bool getCons);
void getMotifGraph(vector<MotifSignalFrame>* m);
const int MAX_WIDTH = 20;
const int MAX_HEIGHT = 3;
const int MAX_HEIGHT_FACTOR = 9;
bool sortByCount (centroid i,centroid j) { return (i.count>j.count); }

int main(){

    vector<MotifSignalFrame>* mt = new vector<MotifSignalFrame>;
  string chromosomes[] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                           "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                           "chr21", "chr22","chrX"};

  int max = 23;
  for(int i=0; i<max; i++){
     getMotifSignal("K562", "K562Sig_filter", chromosomes[i], (*mt), 0, true);

  }
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
  //method 1 : pass it to same template clsasfier
  /*vector<centroid> C;
  generateCentroids(C, 0);
  fpMatch(C, mt, true);
  */
  ///method 2: average the signal for the same motif instance
  getMotifGraph(mt);
  
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
centroid getCentroid(int width, int left_peak, int right_peak)
{
    centroid s;
    int center = WINDOW_SIZE/2;
    int left_half_width = width/2;
    int right_half_width = left_half_width;
    if(width%2 == 1){
      left_half_width++;
    }
    for(int i=0; i< WINDOW_SIZE; i++){
      if(i == center-left_half_width){
        s.signal[i] = left_peak;
      
      }else if(i== center+right_half_width){
        s.signal[i] = right_peak;
      
      }else{
        s.signal[i] = 0;
      }
    }
    return s;
}

void generateCentroids(vector<centroid> &C, int n){
  vector<centroid> S;
  S.push_back(getCentroid(0, 1*MAX_HEIGHT_FACTOR , 1));
  S.push_back(getCentroid(0, 2*MAX_HEIGHT_FACTOR , 1));
  S.push_back(getCentroid(0, 3*MAX_HEIGHT_FACTOR , 1));
  //create one-sided clusters
  for(int i=2; i<=MAX_WIDTH; i++){  
    if(!allowFlip){
      for(int j=1; j<= MAX_HEIGHT; j++){
       for(int k=1; k<=MAX_HEIGHT; k++){
          S.push_back(getCentroid(i, j*MAX_HEIGHT_FACTOR , k*MAX_HEIGHT_FACTOR));
       }
      }
    }else{
      //S.push_back(getCentroid(i, 1*MAX_HEIGHT_FACTOR , 1*MAX_HEIGHT_FACTOR));
      S.push_back(getCentroid(i, 2*MAX_HEIGHT_FACTOR , 2*MAX_HEIGHT_FACTOR));
      //S.push_back(getCentroid(i, 3*MAX_HEIGHT_FACTOR , 3*MAX_HEIGHT_FACTOR));
      S.push_back(getCentroid(i, 3*MAX_HEIGHT_FACTOR , 1*MAX_HEIGHT_FACTOR));
      S.push_back(getCentroid(i, 3*MAX_HEIGHT_FACTOR , 2*MAX_HEIGHT_FACTOR));
      S.push_back(getCentroid(i, 2*MAX_HEIGHT_FACTOR , 1*MAX_HEIGHT_FACTOR));
    }
  }
  C = S;
  return;
  
}


//read the base level conservational energy data for each cluster points
void getMotifConsLevel(string chr, vector<MotifSignalFrame> &temp_con_f, vector<MotifSignalFrame> &f)
{
   cout <<"[DEBUG]reading conservation data for " << chr <<endl;
   string cons_file = "cons_data/"+chr+".phyloP46way.placental.wigFix";
   ifstream Cinfile(cons_file.c_str());
   string a;
   string ss;
   int start;
   Cinfile >> a>> a >> ss >> a;
   
   start = atoi(ss.substr(6).c_str());
   start--;
   //cout << "cons start at " << start <<endl;
   string value;

   vector<MotifSignalFrame>::iterator it;
   while(Cinfile >> value){
      //cout <<"start: " << start << ", cons: " << value <<endl;
      if(value == "fixedStep"){
        //cout <<"hit new start: " << start << endl;
        
        Cinfile >> a >> ss >> a;
        int new_start = atoi(ss.substr(6).c_str());
        new_start--;
       
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
                  f.push_back(*it);
                  temp_con_f.erase(it);
                  it--;
                  if(temp_con_f.size() ==0 )
                     return;
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
                  f.push_back(*it);
                  temp_con_f.erase(it);
                  it--;
                  if(temp_con_f.size() ==0 )
                     return;
               }
         }
      }
      }
      start++;
   }
   Cinfile.close();
   return;
}
//read the footprint file for foorprint sequences
void getMotifSignal(string cell, string signalFile, string chromosome, vector<MotifSignalFrame> &f,double score, bool getCons)
{
   string motifFIle = "hg19_instances/"+chromosome+".hg19_instances";
   vector<MotifSignalFrame> temp_m;
   int s, e;   //start sequence, end sequence
   string chr, mname, strand; //chromosome number, motif name
   cout <<"[DEBUG] get motif instance for " << chromosome<<endl;

   ifstream MTinFile(motifFIle.c_str());
   //NOTE: motif instance interval is inclusive
   int z=0;
   while(MTinFile >> mname >> chr >> s >> e >> strand){
      //get a specific cell type     
      if(chr == chromosome){
           e=e+1;
           MotifSignalFrame mtt;
           mtt.startSeq = s - (FRAME_SIZE-(e-s))/2;
           mtt.length   = 0;
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
   string FPfile = "chr.footprints/"+chromosome+".footprints";
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
   string sFile = signalFile+"/"+chromosome+"."+signalFile;
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
                        //c/out <<remove<<endl;
                        if(getCons){
                           temp_con_f.push_back(*it);
                            real_m.erase(it);
                           //it--;
                           
                           it_start++;
                           if(remove >=real_m.size() ){
                              //read in conservational level for each point
                              getMotifConsLevel(chromosome, temp_con_f, f);
                              return;
                           }
                        }else{
                           f.push_back(*it);
                            real_m.erase(it);
                           it--;
                           if(remove >=real_m.size() )
                              return;
                  
                        }
                     }
                  }
                  
               }
                START++;            
            }
            START = s+1;
         }
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
                  if(getCons){
                     temp_con_f.push_back(*it);
                     //real_m.erase(it);
                     //it--;
                      remove++;
                     it_start++;
                     if(remove >=real_m.size() ){
                        //read in conservational level for each point
                        getMotifConsLevel(chromosome, temp_con_f, f);
                        return;
                     }
                  }else{
                     
                     f.push_back(*it);
                     real_m.erase(it);
                     it--;
                     if(remove >=real_m.size())
                        return;
                  
                  }
               }
            }
            
         }         

      }
      
   }
   
  singalinfile.close();
}
void getMotifGraph(vector<MotifSignalFrame>* mt){

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
  cout<<"[DEBUG]average teh aggregation" <<endl;
  //average the signals and reposition the centroid
  for(int i=0; i<C.size(); i++){
    for(int j=0; j<WINDOW_SIZE; j++){
      C[i].signal[j] = (float)(C[i].signal[j] / C[i].count);
      C[i].cons_signal[j] = (float)(C[i].cons_signal[j] / C[i].count);
    }
  }
  
  cout <<C.size() << " , " << (*mt).size()<<endl;
  cout <<"[DEBUG]sort"<<endl;
  sort (C.begin(), C.end(), sortByCount); 
  cout <<C.size() << " , " << (*mt).size()<<endl;
  
  //print the signals out to data
  cout <<"[DEBUG]output result" <<endl;
  struct stat st = {0};
  //create directory if it doesn't exists
  if (stat("tmp/test_mp_result", &st) == -1) {
    mkdir("tmp/test_mp_result", 0700);
  }
  ofstream outputFile("tmp/test_mp_result/fpsig");
  ofstream consFile("tmp/test_mp_result/consSig");
  ofstream motifFile("tmp/test_mp_result/motifOrder.txt");
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

 
  
}


void fpMatch(vector<centroid> &C, vector<MotifSignalFrame>* f, bool getCons ){
   cout <<"[DEBUG]Assignment mps " <<endl;
  int fpSize = f->size();
  //match each footprint to a template structure
  for(int i=0; i<fpSize; i++){
    int shiftTo = (*f)[i].mtStartIndex;
    bool doFlip = (*f)[i].flip;
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
      s.cons_CI[j] =0;
    }
    temp_C.push_back(s);
    ptCount.push_back(0);
  }
      
      
  //step 2. reposition each centroid to the average of all the points assigned to iterate
  for(int i=0; i<fpSize; i++){
    int clusterIndex = (*f)[i].clusterAssigned;
    int fpIndex      = (*f)[i].mtStartIndex;
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
  for(int i=0; i<fpSize; i++){
    int clusterIndex = (*f)[i].clusterAssigned;
    int fpIndex      = (*f)[i].mtStartIndex;
    float flip_consSignal[FRAME_SIZE];
    
    for(int j=0; j<WINDOW_SIZE; j++){
        temp_C[clusterIndex].cons_CI[j] += (float)(flip_consSignal[j+fpIndex]-(float)temp_C[clusterIndex].cons_signal[j])*
                                           (float)(flip_consSignal[j+fpIndex]-(float)temp_C[clusterIndex].cons_signal[j]);

    }
    
  }
  //finish the CI: caluclaute division by (N-1) and square root
  for(int i=0; i<temp_C.size(); i++){
    //cout <<"index " << i <<" point count: " << ptCount[i]<<endl;
    for(int j=0; j<WINDOW_SIZE; j++){
      temp_C[i].cons_CI[j] =  sqrt(temp_C[i].cons_CI[j]/(ptCount[i]-1));
    }
  }
    

  //print the signals out to data
  cout <<"[DEBUG]output result" <<endl;
  struct stat st = {0};
  //create directory if it doesn't exists
  if (stat("tmp/test_mp", &st) == -1) {
    mkdir("tmp/test_mp", 0700);
  }
  ofstream outputFile("tmp/test_mp/testing_flip_fpsig");
  ofstream consFile("tmp/test_mp/testing_flip_consSig");
  ofstream tempFile("tmp/test_mp/testing_filp_temp");
  ofstream anticFile("tmp/test_mp/anticorrelation_lv.txt");
  ofstream countFile("tmp/test_mp/count.txt");
  ofstream CIFile("tmp/test_mp/cons_CI.txt");
  for(int i=0; i<temp_C.size(); i++){
    outputFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
    tempFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
    CIFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";

    if(getCons){
      consFile <<"\"[DEBUG]index " << i <<" centroid... \"\n";
              
    }
    for(int j=0; j<WINDOW_SIZE; j++){
      outputFile <<j <<" "<<temp_C[i].signal[j] << endl;
      
      //CI range
      CIFile <<j <<" "<<temp_C[i].cons_signal[j]-1.96*temp_C[i].cons_CI[j] 
             << " ~ " <<temp_C[i].cons_signal[j]+1.96*temp_C[i].cons_CI[j] << endl;
      
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
       
       
      C = temp_C;
}
