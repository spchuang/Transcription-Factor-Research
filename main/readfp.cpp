#ifndef READFP_CPP           
#define READFP_CPP 
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
#include "struct.h"
using namespace std;

void getFootPrint(string cell, string chromosome, vector<fpSignalFrame> &f, double score, bool getCons, bool consOffset);
void getConsLevel(string chr, vector<fpSignalFrame> &temp_con_f, vector<fpSignalFrame> &f, bool offset);
void printFPsignals(fpSignalFrame *fp);


//read the base level conservational energy data for each cluster points
void getConsLevel(string chr, vector<fpSignalFrame> &temp_con_f, vector<fpSignalFrame> &f, bool offset)
{
   cout <<"[DEBUG]reading conservation data for " << chr <<endl;
   string cons_file = "../data/cons_data/"+chr+".phyloP46way.placental.wigFix";
   ifstream Cinfile(cons_file.c_str());
   string a;
   string ss;
   int start;
   Cinfile >> a>> a >> ss >> a;
   //cout <<"test1"<<endl;
   cout << ss << endl;
   start = atoi(ss.substr(6).c_str());
   if(offset)
      start--;
   //cout << "cons start at " << start <<endl;
   string value;

   int remove = 0;
   vector<fpSignalFrame>::iterator it;
   vector<fpSignalFrame>::iterator it_start;
   it_start = temp_con_f.begin();
   while(Cinfile >> value){
      //cout <<"start: " << start << ", cons: " << value <<endl;
      if(value == "fixedStep"){
        //cout <<"hit new start: " << start << endl;
        
        Cinfile >> a >> ss >> a;
        //cout <<"test2"<<endl;
        int new_start = atoi(ss.substr(6).c_str());
        if(offset)
           new_start--;
        // <<"new start is ... " << new_start <<endl;
        while(start<new_start){

          for(it=it_start; it<temp_con_f.end(); it++){
            if(it->chr == chr && it->startSeq+it->length > start)
              break;
            if(it->chr == chr && it->length+it->startSeq == start){
              //add 0 for non specified conservation level
              it->con_level[it->length] = 0;
              //cout << "start seq: " << it->startSeq << " , length: " << it->length << ", cons: " << it->con_level[it->length]<<endl;
              it->length++;
              if(it->length >= FRAME_SIZE){
                  f.push_back(*it);
                  //temp_con_f.erase(it);
                  //it--;
                  it_start++;
                  remove++;
                  //if(temp_con_f.size() ==0 )
                  if(remove >=temp_con_f.size() )
                     return;
               }
            }
           
          }
          start++;
       } 
        //cout <<"new start: " << start <<endl;
        continue;
      }
      for(it=it_start; it<temp_con_f.end(); it++){
         if(it->chr == chr && it->startSeq+it->length > start)
            break;
         if(it->chr == chr && it->length+it->startSeq == start){
             
            it->con_level[it->length] = atof(value.c_str());
            //cout << "start seq: " << it->startSeq << " , length: " << it->length << ", cons: " << it->con_level[it->length]<<endl;
            it->length++;
            if(it->length >= FRAME_SIZE){
                  f.push_back(*it);
                  //temp_con_f.erase(it);
                  //it--;
                  it_start++;
                  remove++;
                  //if(temp_con_f.size() ==0 )
                  if(remove >=temp_con_f.size() )
                     return;
               }
         }
      }
      start++;
   }
   Cinfile.close();
   return;
}


//read the footprint file for foorprint sequences
void getFootPrint(string cell, string chromosome, vector<fpSignalFrame> &f,double score, bool getCons, bool consOffset)
{
	//current directory is "main"
   string FPfile = "../data/chr.footprints/"+chromosome+".footprints";
   //temporarily holding the fp data
   vector<fpSignalFrame> temp_f;
   vector<fpSignalFrame> temp_con_f;
   //variables to read in
   int s, e;   //start sequence, end sequence
   double fos; //footprint occupancy score
   string chr, ctype; //chromosome number, cell type
   cout <<"[DEBUG]reading in " << chromosome <<endl;

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
   //READ FP signals
   int START,signal;
   bool first = false;
   //string sFile = signalFile+"/"+chromosome+"."+signalFile;
   string sFile = "../data/signals/"+cell+"_split/"+chromosome+"."+cell;
   ifstream singalinfile(sFile.c_str());
   int remove =0;
   int size = temp_f.size();
   vector<fpSignalFrame>::iterator it_start;
   vector<fpSignalFrame>::iterator it;
   it_start = temp_f.begin();
   //read through the signal data, if the window contains footprint sequence
   //first iteration, find the start of sequence
   while(singalinfile >> chr >> s >> e >> signal){
      if(chr == chromosome){

         if(!first){
            first = true;
            START = s+1;
         }else{
            //skip useless zero signals...
            //int firstFPSeq = ((temp_f.begin())->startSeq)+((temp_f.begin())->length);
            int firstFPSeq = (it_start->startSeq)+(it_start->length);
            if(START < firstFPSeq && firstFPSeq < s)
               START = firstFPSeq;
            else if(START<firstFPSeq && s < firstFPSeq)
               START = s+1;
            //insert zero signal
            while(START < s){
               //iterate through the temp footprint vector and insert the signalFile  
               //it = temp_f.begin()
               for(it=it_start; it<temp_f.end(); it++){
                  //break if the seq is bigger
                  if( (it->startSeq+it->length) > START)
                     break; 

                  //if it's the sequence to enter
                  if( (it->startSeq+ it->length) == START){
                     it->signal[it->length] = 0; //ADD THE EMPTY SIGNAL
                     it->length++;
                     //if the point is full, push it to vector f 
                     if(it->length >= FRAME_SIZE){
                        it->length = 0;
                        remove++;
                        if(getCons){
                           temp_con_f.push_back(*it);
                           //temp_f.erase(it);
                           //it--;
                           it_start++;
                           //if(temp_f.size() ==0 ){
                           if(remove >=size ){
                              //read in conservational level for each point
                         
                              getConsLevel(chromosome, temp_con_f, f, consOffset);
                              return;
                           }
                        }else{
                           f.push_back(*it);
                           temp_f.erase(it);
                           it--;
                           if(temp_f.size() ==0 )
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
         //for(it=temp_f.begin(); it<temp_f.end(); it++){
         for(it=it_start; it<temp_f.end(); it++){
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
                     //temp_f.erase(it);
                     //it--;
                      remove++;
                     it_start++;
                     //if(temp_f.size() ==0 ){
                     if(remove >=size){
                        //read in conservational level for each point
                        getConsLevel(chromosome, temp_con_f, f, consOffset);
                        return;
                     }
                  }else{
                     
                     f.push_back(*it);
                     temp_f.erase(it);
                     it--;
                     if(temp_f.size() ==0 )
                        return;
                  
                  }
               }
            }
            
         }         

      }
      
   }
   
  singalinfile.close();
}

void printFPsignals(fpSignalFrame *fp){
  if(DEBUG){
      cout <<"---------------------" <<endl;       
        
      cout << fp->startSeq <<endl;
      for(int j=0; j<FRAME_SIZE; j++){
        cout <<fp->startSeq+j <<" : " << fp->signal[j] <<endl;
        
      }
      cout <<"-----------------------" <<endl;
   }

}
float correlation(float pSig[], int pIndex, float cSig[])
{
   float pSigSum=0, cSigSum=0;
   for(int i=0; i<WINDOW_SIZE; i++){
      pSigSum+= pSig[i+pIndex];
      cSigSum+= cSig[i];
   }
   pSigSum = pSigSum/(float)WINDOW_SIZE;
   cSigSum = cSigSum/(float)WINDOW_SIZE;
   //cout <<"psum: " <<pSigSum << " , cSum: " << cSigSum <<endl;
   float top=0;
   float stdA=0, stdB=0;
   for(int i=0; i<WINDOW_SIZE; i++){
      top+= ((float)pSig[i+pIndex]-pSigSum)*(cSig[i]-cSigSum);
      stdA+= ((float)pSig[i+pIndex]-pSigSum)*((float)pSig[i+pIndex]-pSigSum);
      stdB+= (cSig[i]-cSigSum)*(cSig[i]-cSigSum);
   }
   return (float)top/sqrt(stdA*stdB);
}
float correlation(int pSig[], int pIndex, float cSig[])
{
   float pSigSum=0, cSigSum=0;
   for(int i=0; i<WINDOW_SIZE; i++){
      pSigSum+= pSig[i+pIndex];
      cSigSum+= cSig[i];
   }
   pSigSum = pSigSum/(float)WINDOW_SIZE;
   cSigSum = cSigSum/(float)WINDOW_SIZE;
   //cout <<"psum: " <<pSigSum << " , cSum: " << cSigSum <<endl;
   float top=0;
   float stdA=0, stdB=0;
   for(int i=0; i<WINDOW_SIZE; i++){
      top+= ((float)pSig[i+pIndex]-pSigSum)*(cSig[i]-cSigSum);
      stdA+= ((float)pSig[i+pIndex]-pSigSum)*((float)pSig[i+pIndex]-pSigSum);
      stdB+= (cSig[i]-cSigSum)*(cSig[i]-cSigSum);
   }
   return (float)top/sqrt(stdA*stdB);
}


//return the eclidean distance ..between a point and a centroid
int distance(int pSig[], int pIndex, float cSig[])
{
    float d = 0;
    for(int i=0; i<WINDOW_SIZE; i++){
        d+= pow((pSig[i+pIndex]- cSig[i]), 2);
    }
    return (int)sqrt(d);
}

#endif 

