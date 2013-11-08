
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


#define min(a,b)  ((a < b) ? a : b)
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


static unsigned int split(const std::string &txt, std::vector<std::string> &strs, char ch)
    {
        unsigned int pos = txt.find( ch );
        unsigned int initialPos = 0;
        strs.clear();

        // Decompose statement
        while( pos != std::string::npos ) {
			cout<< txt.substr( initialPos, pos - initialPos + 1 )<<endl;
            strs.push_back( txt.substr( initialPos, pos - initialPos + 1 ) );
            initialPos = pos + 1;

            pos = txt.find( ch, initialPos );

        }

        // Add the last one
        strs.push_back( txt.substr( initialPos,min( pos, txt.size() ) - initialPos + 1 ) );

        return strs.size();
    }
    
centroid readTemplate(ifstream& centroidFile){
	string line;
	
	if(!centroidFile.eof()){
		getline(centroidFile,line);
	}


	stringstream stream(line);
		
	//cout<<'X'<<line<<'X'<<endl;
	vector<string> vals;
	/*if(split(line, vals, ' ') !=WINDOW_SIZE){
		cerr <<"the read is wrong\n";
		exit(1);
	}
	cout<<vals.size()<<endl;*/
	centroid s;
	int value;
	int i=0;
	while(stream >> value && i<WINDOW_SIZE){
		
		
		s.signal[i]=(int)value;
		i++;
		stream.get();
	}
	
	/*
	for(int i=0; i<vals.size(); i++){
		s.signal[i] = atoi(vals[i].c_str());
	}*/
	
  return s;
}

void saveFirstTemplate(centroid s, string name, string dirname){
  string filename = dirname+"/"+ name;
  ofstream centroidFile(filename.c_str(), ios_base::trunc );
  //cout <<"saving.."<<endl;
  for(int i=0; i< WINDOW_SIZE; i++)
  {
  	//cout <<s.signal[i]<<endl;
  	centroidFile << s.signal[i] << " ";
  }
}

void pushTemplateToFile(centroid s, ofstream& centroidFile){
  for(int i=0; i< WINDOW_SIZE; i++)
  {
  	//cout <<s.signal[i]<<endl;
  	centroidFile << s.signal[i] << " ";
  }
  centroidFile << "\n";
}

void generateTemplate( int n, int start, int length, string dirname){
 
  //n degree of movement
  //y>=0
  //start = end = 0
  
  //start with 
  centroid s;
  for(int i=0; i<WINDOW_SIZE; i++){
    s.signal[i]=0;
  }
  s.signal[start] = 0;
  long long int size = 1;
  
  string old_file = "a";
  string new_file = "b";
  saveFirstTemplate(s, old_file, dirname);

  
  //S.push_back(s);
  for(int i=0; i<=length; i++){
  	size=0;
  	cout<<"Loop: "<<i<<endl;
  	//open a new stream for next round of loop
  	string new_file_name = dirname +"/"+new_file;
  	string old_file_name = dirname +"/"+old_file;
  	cout <<"newFile: "<<new_file_name<<endl;
  	cout <<"oldFile: "<<old_file_name<<endl;
  	ofstream writeCentroidFile(new_file_name.c_str(), ios_base::trunc );
  	ifstream readCentroidFile(old_file_name.c_str());
  	
  	//read from each centroid from old file then write to new file
  	while(!readCentroidFile.eof()){
  		//read in centroid
  		centroid s = readTemplate(readCentroidFile);
  		
  		//different levels
  		for(int z=n; z>=(n*-1);z--){
      		centroid ss = s;

  			int next_val=ss.signal[i-1+start] + MAX_TEMPLATE_HEIGHT*z;
  			if(i==length && next_val!=0) continue;
  			ss.signal[i+start]= next_val;
        	/*
        	for(int g=0; g<WINDOW_SIZE;g++)
  				cout <<ss.signal[g] <<" ";
  		cout <<endl;
  		*/
  			pushTemplateToFile(ss, writeCentroidFile);
        
  			size++;
  		}		
  		
  		
  		
  	}
  	
  	//close the file streams
  	writeCentroidFile.close();
  	readCentroidFile.close();
  	
  	//swap old and new file
  	if( remove(old_file_name.c_str() ) != 0 ){
	    cout << "Error deleting file: " <<  old_file_name << endl;
	 }
	string tmp;
	tmp = old_file;
	old_file = new_file;
	new_file = tmp;
  	
  	/*
 	cout << "size: " << size <<endl;
    long long int new_size=0;
    
    for(int j= 0; j<size; j++){
      //different levels
      for(int z=n; z>=(n*-1);z--){
      	//read in centroid
      	
      	
      	centroid ss = readTemplate(j, dirname);
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
    */
  }
  
  cout <<"[DEBUG}Total template size: " << size<<endl;

}

void reduceVectorFromFile(vector<centroid>& corr_result, ifstream& centroidFile, int n){
	//reset
	centroidFile.clear();                 // clear fail and eof bits
	centroidFile.seekg(0, std::ios::beg); // back to the start!

  //pick one points
  corr_result.clear();
  corr_result.push_back(readTemplate(centroidFile));
  vector<int> centroid_index;
  centroid_index.push_back(0);
  cout << "[DEBUG]creating: 1" <<endl;


  int count = 2;
  //iterate until fine one that maximizes the minimum distace
  while(count <= n){
    cout << "[DEBUG]creating: "<<count <<endl;
    //reset
	centroidFile.clear();                 // clear fail and eof bits
	centroidFile.seekg(0, std::ios::beg); // back to the start!
	int i=0;
	vector<float> tracker;
	while(!centroidFile.eof()){
		centroid s = readTemplate(centroidFile);
		tracker.push_back(FLT_MAX);
		//check if the centroid is already in the list
		bool centroidAdded = false;
		for(int z=0; z<centroid_index.size(); z++){
			if(centroid_index[z] == i){
				centroidAdded = true;
				tracker[i] = -1;
			}
			
		}
		if(centroidAdded){
			continue;
		}
		
		
		for(int j=0; j<corr_result.size(); j++){
	          float d = 1-correlation(s.signal, 0, corr_result[j].signal);
	          if(d < tracker[i]){
	            tracker[i] = d;
	          }
	    }
		i++;
	}
    
    
    float max = FLT_MIN;
	int add_i = -1;
	//find maximum of minimum
	for(int j=0; j<tracker.size();j++){
		if(tracker[j]!=-1 && tracker[j] > max){
		  max = tracker[j];
		  add_i = j;
		}
	}
	//reset
	centroidFile.clear();                 // clear fail and eof bits
	centroidFile.seekg(0, std::ios::beg); // back to the start!
	
	
	//add centroid index to result list
	int c_index=0;
	while(!centroidFile.eof()){
		centroid s = readTemplate(centroidFile);
		if(c_index == add_i){
			corr_result.push_back(s);
			centroid_index.push_back(c_index);
			break;
		}

		c_index++;
	}
    
    count++;
  }

  
}

//use recurssion
void generateTemplateWIthRescurrsion(centroid s,int n, int start, int end, ofstream& centroidFile){

	//end condition
	if(start ==  end){
		//converges back then we save it as a potential template
		if(s.signal[end] ==0){
			/*cout<<"one-----------------"<<endl;
			for(int i=0; i<WINDOW_SIZE; i++)
				cout <<s.signal[i] << " ";
			cout <<endl;*/
			pushTemplateToFile(s, centroidFile);
			return;
			
		
		}else{
			return;
		}
	}

	for(int z=1; z>=(n*-1);z--){
		s.signal[start+1] = s.signal[start] + MAX_TEMPLATE_HEIGHT*z;
		generateTemplateWIthRescurrsion(s, n, start+1, end, centroidFile);
	
	}



}

int main()
{
	string dirname = "../tmp/output_template_9_4";
	struct stat st = {0};
   	//create directory if it doesn't exists
	if (stat(dirname.c_str(), & st) == -1) {
		if(mkdir(dirname.c_str(), 0700) == -1){
			cout <<"[DEBUG]Error creating output template folder" <<endl;
		}
		
	}
	string tmpDir = dirname +"/tmp";
	if (stat(tmpDir.c_str(), & st) == -1) {
		if(mkdir(tmpDir.c_str(), 0700) == -1){
			cout <<"[DEBUG]Error creating output template folder" <<endl;
		}
		
	}

  vector<centroid> C;
  vector<fpSignalFrame>* fpF = new vector<fpSignalFrame>;
  
  centroid s;
  for(int i=0;i<WINDOW_SIZE;i++){
  	s.signal[i] = 0;
  }
  string new_file_name = tmpDir+"/test";
  ofstream writeCentroidFile(new_file_name.c_str(), ios_base::trunc );
  generateTemplateWIthRescurrsion(s,1,1,28,writeCentroidFile);
  //generateTemplate(C, 1, 7, 19);
  //generateTemplate(1, 7, 19, tmpDir);
  //generateTemplate(1, 1, 28, tmpDir);

  //generateTemplate(C, 1, 9 , 13);
  //reduceVector(C,100);
	for(int i=0; i< C.size(); i++){
	
	}

}
