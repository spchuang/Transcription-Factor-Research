#ifndef CONFIG_H    
#define CONFIG_H 

const int WINDOW_SIZE = 30;
const int FRAME_SIZE  = 30;
      
struct segment{
   float con_level[FRAME_SIZE];
   int 	 signal[FRAME_SIZE];       //holding the signal level for each basepair in the frame
   int   segStart;
   int   frameStart;                 //the start DNA location for the frame (including the actual footprint sequence)
   int   length;                   //used for grabbing signal data
   char* chr;                    //the chromosome
   char* cell;
   char* name;
   int   segLength;                 //length of the footprint (from the fp data)
   int   startIndex;             //the offset for the fixed window in the frame
   bool flip;
   int segIndex;
};




#endif