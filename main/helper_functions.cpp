#ifndef HELPER_FUNCTION_CPP           
#define HELPER_FUNCTION_CPP 

#include <cmath>
/*
   return the correlation of two list
   [in]pSig: list 1
   [in]cSig: list 2
   [in]size: size of two lists (usually WINDOW_SIZE)
   [int]pIndex: offset for pSig (since sometimes psig could be a longer array then cSig)
*/
float correlation(float pSig[], int pIndex, float cSig[], int size)
{
   float pSigSum=0, cSigSum=0;
   for(int i=0; i<size; i++){
      pSigSum+= pSig[i+pIndex];
      cSigSum+= cSig[i];
   }
   pSigSum = pSigSum/(float)size;
   cSigSum = cSigSum/(float)size;
   //cout <<"psum: " <<pSigSum << " , cSum: " << cSigSum <<endl;
   float top=0;
   float stdA=0, stdB=0;
   for(int i=0; i<size; i++){
      top+= ((float)pSig[i+pIndex]-pSigSum)*(cSig[i]-cSigSum);
      stdA+= ((float)pSig[i+pIndex]-pSigSum)*((float)pSig[i+pIndex]-pSigSum);
      stdB+= (cSig[i]-cSigSum)*(cSig[i]-cSigSum);
   }
   return (float)top/sqrt(stdA*stdB);
}
//integer version
float correlation(int pSig[], int pIndex, int cSig[], int size)
{
   float pSigSum=0, cSigSum=0;
   for(int i=0; i<size; i++){
      pSigSum+= pSig[i+pIndex];
      cSigSum+= cSig[i];
   }
   pSigSum = pSigSum/(float)size;
   cSigSum = cSigSum/(float)size;
   //cout <<"psum: " <<pSigSum << " , cSum: " << cSigSum <<endl;
   float top=0;
   float stdA=0, stdB=0;
   for(int i=0; i<size; i++){
      top+= ((float)pSig[i+pIndex]-pSigSum)*(cSig[i]-cSigSum);
      stdA+= ((float)pSig[i+pIndex]-pSigSum)*((float)pSig[i+pIndex]-pSigSum);
      stdB+= (cSig[i]-cSigSum)*(cSig[i]-cSigSum);
   }
   return (float)top/sqrt(stdA*stdB);
}

void print_signals(int pSig[], int size){
   cout <<"Print signal" <<endl;
   for(int i=0; i<size; i++){
      cout << pSig[i]<< " ";
   }
   cout <<endl;
   
}

#endif