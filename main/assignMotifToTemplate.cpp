#include "readMotifMaster.cpp"
#include "templatesUtility.cpp"


void saveMotifAssignment(string dirname, vector<template_segment> &templates, int total_assigned){
   
	string filename = dirname+"/fpsig";
	cout << filename <<endl;
	ofstream outputFile(filename.c_str(), ios_base::trunc );
	
	/*
	filename = dirname+"/consSig";
	ofstream consFile(filename.c_str(), ios_base::trunc );
	filename = dirname+"/motifOrder";
	ofstream motifFile(filename.c_str(), ios_base::trunc );*/
	
	
	for(int i=0; i<templates.size(); i++){
	   float ratio = (templates[i].assigned_count/(float)total_assigned)*100;
	   float corr  = templates[i].assigned_correlation;
	   
	
      outputFile << i << " " << ratio << " " <<corr<<endl;      	
      for(int j=0; j<FRAME_SIZE; j++){
         outputFile << templates[i].assigned_signal[j] << endl;	
		}
   	
	}
	
   outputFile.close();
   
}

void createMotifAssignmentFolders(string rootDir, string motif_name, string cell){
	struct stat st = {0};
   //TF/Cell/fpsig -> p1...p100
   
   //create the motif/ directory
   string dirname = rootDir+"/motif_assignment/"+motif_name;
	
	cout <<"[DEBUG]create direcotry:" + dirname <<endl;
	if (stat(dirname.c_str(), &st) == -1) {
		if (mkdir(dirname.c_str(), 0700) != 0){
			cout <<"PROBLEM FUCK.."<<endl;
		}
	}
	
	//create the motif/cell directory
	dirname = rootDir+"/motif_assignment/"+motif_name+"/"+cell;
	
	cout <<"[DEBUG]create direcotry:" + dirname <<endl;
	if (stat(dirname.c_str(), &st) == -1) {
		if (mkdir(dirname.c_str(), 0700) != 0){
			cout <<"PROBLEM FUCK.."<<endl;
		}
	}
}

int main(){
   /*
      Steps:
      1)generate the 100 profiles (done: templatesUtility.cpp)
      2)read in all the motifs for the cell type
      3)find out all the unique motifs
      4)for each instances of unique motifs, assign it to a profile (add the signals and increment the count)
      5)average the signals for that motif for each profile
   
   NOTE: naive solution with slow efficiency
   */


	string cellTypes[] = {"AG10803", "AoAF", "CD20+", "GM06990", "GM12865","H7-hESC","HAEpiC","HA-h","HCF","HCM","HCPEpiC","HEEpiC","HepG2","HFF","HIPEpiC","HMF","HMVEC-dBl-Ad","HPAF","HPdLF","HPF","HRCEpiC","HSMM","HVMF","K562","NB4","NH-A","NHDF-Ad","NHDF-neo","NHLF","SAEC","SKMC","Th1"};
	string chromosomes[] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                           "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                           "chr21", "chr22","chrX"};
	int totalChr = 23;//23;     
   int totalCell = 31;//31;       
   string rootDir = "all_chr"; 
   
   //generate the 100 profiles
   TemplateGenerator TG;
   TG.readTemplatesFromFile("template_profile_100");
   
   ReadMotifMaster ReadMaster(rootDir);
    
   for(int cell=0; cell<totalCell; cell++){
      //read in all the motifs for the cell type
      ReadMaster.readMotif(cellTypes[cell], chromosomes, totalChr);
      
      //find out all the unque motifs
      ReadMaster.computeUniqueMotifs();
   
      
      //for each unique motifs
      for(int i=0; i<ReadMaster.unique_motifs.size(); i++){
         //reset the template assigned values
         TG.initializeAssigned();
        
         //for each instances of unique motifs, assign it to a profile (add the signals and increment the count)
         int count =0;
         for(int j=0; j<ReadMaster.segs.size(); j++){
        
            if(strcmp(ReadMaster.segs[j].name, ReadMaster.unique_motifs[i].name) == 0){
               count++;
               TG.assignToTemplate(ReadMaster.segs[j]);
            }
         }
         cout <<"Total: " <<count<<endl;
         
         //average all the assigned signals
         TG.averageAssignedSignals();
         
         //OUTPUT
         string motif_name = string(ReadMaster.unique_motifs[i].name);
         createMotifAssignmentFolders(rootDir, motif_name, cellTypes[cell]);
         
         string saveDir = rootDir+"/motif_assignment/"+motif_name+"/"+cellTypes[cell];
         saveMotifAssignment(saveDir, TG.templates, TG.total_assigned);  
  
      }
      

   }
    			
}