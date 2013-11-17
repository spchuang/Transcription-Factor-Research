#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <list>

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
   /* cout <<endl<<endl<<endl;
    cout << count <<endl;
    cout << S.size() <<" , " <<R.size()<<endl;*/
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
  /*cout <<"subset size: " << R.size()<<endl;
  for(int i=0; i<R.size(); i++){
    for(int j=0; j<WINDOW_SIZE; j++){
      cout << R[i].signal[j] << " ";
    
    }
    cout <<endl;
  
  }
  exit(1);*/
  S = R;
  
}

void readTemplateFromFile(vector<centroid> &S ,ifstream& centroidFile)
{
	S.clear();
	int i=0;
	centroid s;
	while(!centroidFile.eof()){
		string a;
		getline(centroidFile,a);
		
		if(a[0]!='"'){
			stringstream stream(a);
			
			int junk, val;
			stream >> junk >> val;
		

			s.signal[i]=(int)val;
			stream.get();
			
			i++;
			
			if(i==30){
				i=0;

				S.push_back(s);	
			} 
		}
			
		
		
	}
	
	

}

void generateTemplate(vector<centroid> &S, int n, int start, int length){
  S.clear();
  //n degree of movement
  //y>=0
  //start = end = 0
  
  centroid s;
  for(int i=0; i<WINDOW_SIZE; i++){
    s.signal[i]=0;
  }
  s.signal[start] = 2*MAX_TEMPLATE_HEIGHT;
  S.push_back(s);
  for(int i=1; i<=length; i++){
 
    int size = S.size();
    int number=0;
    for(int j= 0; j<size; j++){
      //different levels
      for(int z=n; z>=(n*-1);z--){
        centroid ss = S[j];
        int next_val=ss.signal[i-1+start] + MAX_TEMPLATE_HEIGHT*z;
        
        //if the value is less than 0
        //if(next_val<0) continue;
        if(i==length && next_val!=0) continue;
        ss.signal[i+start]= next_val;
        number++;
       // cout <<"pushing: " <<S.capacity()<< ", " <<S.size() <<endl;
        S.push_back(ss);
        
      }
    }
    
    //swap new content to front
    for(int t=0; t<number; t++){
      S[t] = S[size+t];
    }
    //cout <<"before: " <<S.capacity()<< ", " <<S.size() <<endl;
    while(S.size() > number){
      S.pop_back();
    }
    
   // cout <<"after: " <<S.capacity()<< ", " << S.size() <<endl;
  }
  
  cout <<"[DEBUG}Total template size: " << S.size()<<endl;
  /*for(int i=0; 3; i++){
    bool test = true;
    
      for(int j=0; j<WINDOW_SIZE;j++)
         cout << S[i].signal[j] << " ";
 
    cout <<endl;
  }*/
  
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

centroid getCentroid2(int width, int *left_struc, int left_l, int *right_struc, int right_l)
{
    centroid s;
    int center = WINDOW_SIZE/2;
    int left_half_width = width/2;
    int right_half_width = left_half_width;
    if(width%2 == 1){
      left_half_width++;
    }
    int l=0;
    int r=0;
    for(int i=0; i< WINDOW_SIZE; i++){
      if(left_l>0 && i == center-left_half_width-left_l){
        s.signal[i] = left_struc[l];
        l++;
        left_l--;
      
      }else if(r<right_l && i== center+right_half_width+r){
        s.signal[i] = right_struc[r];
        r++;
      }else{
        s.signal[i] = 0;
      }
    }
    return s;


}


void generateCentroids(vector<centroid> &S, int n){
  S.clear();
  S.push_back(getCentroid(0, 1*MAX_HEIGHT_FACTOR , 1));
  S.push_back(getCentroid(0, 3*MAX_HEIGHT_FACTOR , 1));
  
  
  int s[3][1];
  for(int i=1; i<=3; i++){
    s[i-1][0] = i*MAX_HEIGHT_FACTOR;
    
  }
  //structure for 2 consecutive signals
  int d[9][2] ;
  for(int i=1; i<=3; i++){
    for(int k=1; k<=3; k++){
      d[(i-1)*3+k-1][0] = i*MAX_HEIGHT_FACTOR;
      d[(i-1)*3+k-1][1] = k*MAX_HEIGHT_FACTOR;
    }
  }  
  
  
  int t[27][3];
  for(int i=1; i<=3; i++){
    for(int j=1; j<=3; j++){
      for(int k=1; k<=3; k++){
        t[(i-1)*9+(j-1)*3+k-1][0] = i*MAX_HEIGHT_FACTOR;
        t[(i-1)*9+(j-1)*3+k-1][1] = j*MAX_HEIGHT_FACTOR;
        t[(i-1)*9+(j-1)*3+k-1][2] = k*MAX_HEIGHT_FACTOR;
      }
    }
  }
  
  int f[81][4];
  for(int i=1; i<=3; i++){
    for(int j=1; j<=3; j++){
      for(int k=1; k<=3; k++){
        for(int z=1; z<=3; z++){ 
          f[(i-1)*27+(j-1)*9+(k-1)*3+(z-1)][0] = i*MAX_HEIGHT_FACTOR;
          f[(i-1)*27+(j-1)*9+(k-1)*3+(z-1)][1] = j*MAX_HEIGHT_FACTOR;
          f[(i-1)*27+(j-1)*9+(k-1)*3+(z-1)][2] = k*MAX_HEIGHT_FACTOR;
          f[(i-1)*27+(j-1)*9+(k-1)*3+(z-1)][3] = z*MAX_HEIGHT_FACTOR;
        }
      }
    }
  }
  /*for(int i=0; i< 81; i++){
    cout <<f[i][0]<<" "<<f[i][1]<<" "<<f[i][2]<<" "<<f[i][3]<<endl;
  
  }
  
  exit(1);*/

  //structure for 3 consecutie signals;
  
  //create one-sided clusters
  for(int w=2; w<=MAX_WIDTH; w++){  
    if(!allowFlip){
      for(int j=1; j<= MAX_HEIGHT; j++){
       for(int k=1; k<=MAX_HEIGHT; k++){
          S.push_back(getCentroid(w, j*MAX_HEIGHT_FACTOR , k*MAX_HEIGHT_FACTOR));
       }
      }
    }else{
      for(int i=0; i< 120 ; i++){
        //for single struct
        
        for(int j=0; j<120; j++){
          if(i<3 && j<3){
            S.push_back(getCentroid2(w, s[i], 1 , s[j], 1));
          
          }else if(i<3&& j<12){
            S.push_back(getCentroid2(w, s[i], 1 , d[j-3], 2));

          }else if(i<3 && j<39){
            S.push_back(getCentroid2(w, s[i], 1 , t[j-12], 3));
          
          }else if(i<3 && j<120){
            S.push_back(getCentroid2(w, s[i], 1 , f[j-39], 4));
          
          }else if(i<12&& j<3){
            S.push_back(getCentroid2(w, d[i-3], 2 , s[j], 1));
            
          }else if(i<12 && j<12){
            S.push_back(getCentroid2(w, d[i-3], 2 , d[j-3], 2)); 
          }else if(i<12&& j<39){
            S.push_back(getCentroid2(w, d[i-3], 2 , t[j-12], 3));
          
          }else if(i<12&& j<120){
            S.push_back(getCentroid2(w, d[i-3], 2 , f[j-39], 4));
          
          }else if(i<39 && j<3){
            S.push_back(getCentroid2(w, t[i-12], 3 , s[j], 1));
          
          }else if(i<39&& j<12){
            S.push_back(getCentroid2(w, t[i-12], 3 , d[j-3], 2));
          
          }else if(i<39 && j<39){
            S.push_back(getCentroid2(w, t[i-12], 3 , t[j-12], 3));
          
          }else if(i<39 && j<120){
            S.push_back(getCentroid2(w, t[i-12], 3 , f[j-39], 4));
          
          }else if(i<120 && j<3){
              S.push_back(getCentroid2(w, f[i-39], 4 , s[j], 1));
          
          
          }else if(i<120 && j<12){
              S.push_back(getCentroid2(w, f[i-39], 4 , d[j-3], 2));
          
          
          }else if(i<120 && j<39){
              S.push_back(getCentroid2(w, f[i-39], 4 , t[j-12], 3));
          
          
          }else if(i<120 && j<120){
              S.push_back(getCentroid2(w, f[i-39], 4 , f[j-39], 4));
          }
             
        }      
      
      }
      //S.push_back(getCentroid(i, 1*MAX_HEIGHT_FACTOR , 1*MAX_HEIGHT_FACTOR));
      /*S.push_back(getCentroid(i, 2*MAX_HEIGHT_FACTOR , 2*MAX_HEIGHT_FACTOR));
      //S.push_back(getCentroid(i, 3*MAX_HEIGHT_FACTOR , 3*MAX_HEIGHT_FACTOR));
      S.push_back(getCentroid(i, 3*MAX_HEIGHT_FACTOR , 1*MAX_HEIGHT_FACTOR));
      S.push_back(getCentroid(i, 3*MAX_HEIGHT_FACTOR , 2*MAX_HEIGHT_FACTOR));
      S.push_back(getCentroid(i, 2*MAX_HEIGHT_FACTOR , 1*MAX_HEIGHT_FACTOR));*/

    }
  }
    cout <<"size: " << S.size()<<endl;
     return;
  
}
