#ifndef TEMPLATE_UTILITY_CPP           
#define TEMPLATE_UTILITY_CPP

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
#include <sys/stat.h>
#include "config.h"
#include "helper_functions.cpp"

using namespace std;

struct template_segment{
   int signal[FRAME_SIZE];
   float assigned_signal[FRAME_SIZE];
   float assigned_correlation;
   int assigned_count;
};

class TemplateGenerator
{
  public:                    
    TemplateGenerator();     
    ~TemplateGenerator();                 
	void generateTemplates(int n);
	void readTemplatesFromFile(string src);
	void printTemplates();
	void assignToTemplate(segment s);
	void averageAssignedSignals();
	void initializeAssigned();
	
	vector<template_segment> templates;
	int total_assigned;
    
 	
 	

};

TemplateGenerator::TemplateGenerator(){

}

TemplateGenerator::~TemplateGenerator(){
   
}

void TemplateGenerator::generateTemplates(int n){
   //TODO
}

void TemplateGenerator::initializeAssigned(){
   for(int i=0; i<templates.size(); i++){
      if(templates[i].assigned_count >0){
         
         for(int j=0; j<FRAME_SIZE; j++){
   			templates[i].assigned_signal[j] = 0;
   		}
   		templates[i].assigned_correlation = 0;
   		templates[i].assigned_count = 0;
      }
		
   }
   total_assigned=0;
}

void TemplateGenerator::readTemplatesFromFile(string src){
   ifstream templateFile(src.c_str());
   int index, value;
   template_segment t;
   
   while(templateFile >> index >> value){
      t.signal[index] = value;
      t.assigned_signal[index] = 0;
      
      //push the new segment to vector
      if(index == FRAME_SIZE-1){
         t.assigned_count = 0;
         t.assigned_correlation = 0;
         templates.push_back(t);   
         //put in empty templates 
      }
   }
   cout << "[DEBUG]Total template: " << templates.size() << endl;
}

void TemplateGenerator::printTemplates(){
   for(int i=0; i<templates.size(); i++){
      //show percentage
      cout <<"Template " << i << " ( " <<(templates[i].assigned_count/(float)total_assigned)*100<< " , " << templates[i].assigned_correlation << " ) " <<endl;
      cout <<"Profile:  ";
      for(int w=0; w<FRAME_SIZE; w++){
            cout << templates[i].signal[w] << " ";
      }
      cout <<endl;
      cout <<"Assigned: ";
      for(int w=0; w<FRAME_SIZE; w++){
            cout << templates[i].assigned_signal[w] << " ";
      }
      cout <<endl;
   }
   
}

void TemplateGenerator::assignToTemplate(segment s){

   float max_d = FLT_MIN;
   int max_index = -1;
   
   for(int i=0; i<templates.size(); i++){
   
      float d = correlation(templates[i].signal, 0, s.signal, FRAME_SIZE);
      /*
      print_signals(templates[i].signal,FRAME_SIZE);
      print_signals(s.signal, FRAME_SIZE);
      if(isnan(d)){
         cout <<"NOT VALID"<<endl;
      }else{
         cout << "correlation: " << d<<endl;
      }
      */
      //if distance is smaller than prev min
      
      //TODO: there could be a minimum correlation level
      if(d > max_d){
         max_d = d;
         max_index = i;
      }
      
   }
   //if there is a valid correlation
   if(max_index != -1){
      for(int i=0; i<FRAME_SIZE; i++){
         templates[max_index].assigned_signal[i] += s.signal[i];
      }
      templates[max_index].assigned_correlation+=max_d;
      templates[max_index].assigned_count++;
      total_assigned++;
   }
}

void TemplateGenerator::averageAssignedSignals(){
   int total_count =0;
   for(int i=0; i<templates.size(); i++){
      //only calculate the average if there are actually assigned signals
      if(templates[i].assigned_count >0){
         
         for(int j=0; j<FRAME_SIZE; j++){
   			templates[i].assigned_signal[j] = (float)(templates[i].assigned_signal[j] / templates[i].assigned_count);
   		}
   		templates[i].assigned_correlation = (float)(templates[i].assigned_correlation / templates[i].assigned_count);   		
      }
		
   }
   cout <<"[DEBUG]total averaged motifs: " <<total_assigned<<endl;
   
}

/*
int main(){
   TemplateGenerator a;
   a.readTemplatesFromFile("template_profile_100");
   a.printTempaltes();
   
}*/



#endif