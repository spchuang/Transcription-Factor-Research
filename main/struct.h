#ifndef STRUCT_H            
#define STRUCT_H          

const int WINDOW_SIZE = 30;
const int FRAME_SIZE  = 60;
const int MAX_OFFSET  = 30;
const int MAX_ITERATION = 1000;

bool DEBUG = false;

//footprint data and the signals
//This represents a single point in the domain that will be used for clustering
struct fpSignalFrame{
   float con_level[FRAME_SIZE];
   int signal[FRAME_SIZE];       //holding the signal level for each basepair in the frame
   int fpStart;
   int startSeq;                 //the start DNA location for the frame (including the actual footprint sequence)
   int length;                   //used for grabbing signal data
   char* chr;                    //the chromosome
   int fpLength;                 //length of the footprint (from the fp data)
   int fpStartIndex;             //the offset for the fixed window in the frame
   int clusterAssigned;          //the cluster this point is assigned to
   bool flip;
};

struct centroid{
   float signal[WINDOW_SIZE];      //signal level for the cluster of fixed window size
   float cons_signal[WINDOW_SIZE];
  // float cons_CI[WINDOW_SIZE];
   
   //for motif clustering
   //char* motifName;
   //int count;
};

const int MAX_WIDTH = 20;
const int MAX_HEIGHT = 6;
const int MAX_HEIGHT_FACTOR = 6;
const int MAX_TEMPLATE_HEIGHT = 3;
bool allowFlip = true;
#endif 
