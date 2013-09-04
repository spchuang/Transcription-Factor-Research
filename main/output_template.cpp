#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <list>

#include <sys/types.h>
#include <sys/stat.h>
#include <sstream>
#include "readfp.cpp"

using namespace std;

//reduce the vector S into a most represented subset of size n
void reduceVector(vector<centroid> &S, int n){

  //pick one points
  vector<centroid> R;
  R.push_back(S[0]);
  S[0] = S[S.size()-1];
  S.pop_back();
  int count = 2;
  //iterate until fine one that maximizes the minimum distace
  while(count <= n){
    cout << "[DEBUG]creating: "<<count <<endl;
    vector<float> tracker;
   
    for(int i=0; i<S.size(); i++){
      tracker.push_back(FLT_MAX);
      for(int j=0; j<R.size(); j++){
          float d = 1-correlation(S[i].signal, 0, R[j].signal);
          if(d < tracker[i]){
            tracker[i] = d;
          }
      }
    }
    float max = FLT_MIN;
      int add_i = -1;
      //find maximum of minimum
      for(int j=0; j<tracker.size();j++){
        if(tracker[j] > max){
          max = tracker[j];
          add_i = j;
        }
      }
    R.push_back(S[add_i]);
    S.erase(S.begin()+add_i);
    count++;
  }
  S = R;
  
}

centroid readTemplate(long long int index, string dirname){
	ostringstream convert;   // stream used for the conversion

  convert << index;      // insert the textual representation of 'Number' in the characters in the stream
  string filename = dirname+"/"+convert.str();
  centroid s;
  int value;
  int i=0;
  ifstream centroidFile(filename.c_str());
  while(centroidFile >> value && i<WINDOW_SIZE){
  	
    s.signal[i]=(int)value;
    //cout << value <<endl;
    i++;
  }
  return s;
}

void saveTemplate(centroid s, long long int index, string dirname){
  ostringstream convert;   // stream used for the conversion

  convert << index;      // insert the textual representation of 'Number' in the characters in the stream
  string filename = dirname+"/"+convert.str();
  ofstream centroidFile(filename.c_str(), ios_base::trunc );
  //cout <<"saving.."<<endl;
  for(int i=0; i< WINDOW_SIZE; i++)
  {
  	//cout <<s.signal[i]<<endl;
  	centroidFile << s.signal[i] << " ";
  }

}

void generateTemplate( int n, int start, int length){
 
  //n degree of movement
  //y>=0
  //start = end = 0
  string dirname = "../tmp/output_template/tmp";
  
  centroid s;
  for(int i=0; i<WINDOW_SIZE; i++){
    s.signal[i]=0;
  }
  s.signal[start] = 2*MAX_TEMPLATE_HEIGHT;
  long long int size = 1;
  saveTemplate(s, 0, dirname);
  
  //S.push_back(s);
  for(int i=1; i<=length; i++){
 	cout << "size: " << size <<endl;
    long long int new_size=0;
    for(int j= 0; j<size; j++){
      //different levels
      for(int z=n; z>=(n*-1);z--){
      	//read in centroid
      	
      	
      	centroid ss = readTemplate(j, dirname);
      /*	
      	cout <<endl<<"signals read.."<<endl;
      	for(int kk=0; kk< WINDOW_SIZE; kk++)
  {
  	cout <<ss.signal[kk]<<endl;
  	;
  }
  cout <<endl<<endl<<".. .."<<endl;*/
        //centroid ss = S[j];
        int next_val=ss.signal[i-1+start] + MAX_TEMPLATE_HEIGHT*z;
        
        if(i==length && next_val!=0) continue;
        ss.signal[i+start]= next_val;
       
        
        //S.push_back(ss);
        saveTemplate(ss, size+new_size, dirname);
        
        new_size++;
      }
    }
   
    //swap new content to front
    for(int t=0; t<new_size; t++){
      //S[t] = S[size+t];
      centroid ss = readTemplate(size+t, dirname);
      saveTemplate(ss, t, dirname);
      //rename("oldname.txt", "newname.txt");
    }
    
    
    while(size){
      //S.pop_back();
      ostringstream convert;   // stream used for the conversion
 	  convert << (new_size+size-1);      // insert the textual representation of 'Number' in the characters in the stream
  
      string filename = dirname+"/"+convert.str();
      //cout <<"removing: " <<filename<<endl;
      if( remove(filename.c_str() ) != 0 ){
	    cout << "Error deleting file: " <<  filename << endl;
	  }
	  size--;
    }
    
    size = new_size;
  }
  
  cout <<"[DEBUG}Total template size: " << size<<endl;

}


int main()
{

struct stat st = {0};
   	//create directory if it doesn't exists
	if (stat("../tmp/output_template", & st) == -1) {
		if(mkdir("../tmp/output_template", 0700) == -1){
			cout <<"[DEBUG]Error creating output template folder" <<endl;
		}
		
	}
	if (stat("../tmp/output_template/tmp", & st) == -1) {
		if(mkdir("../tmp/output_template/tmp", 0700) == -1){
			cout <<"[DEBUG]Error creating output template folder" <<endl;
		}
		
	}

  vector<centroid> C;
  vector<fpSignalFrame>* fpF = new vector<fpSignalFrame>;
  //generateTemplate(C, 1, 7, 19);
  generateTemplate(1, 7, 19);
  //generateTemplate(C, 1, 9 , 13);
  //reduceVector(C,100);
	for(int i=0; i< C.size(); i++){
	
	}

}

