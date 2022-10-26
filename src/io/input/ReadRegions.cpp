#include "FemTech.h"
#include "utilities.h"
#include <assert.h>

int *region_array;
int mel, sel, cel, fel, pel, oel, tel;
//-------------------------------------------------------------------------------------------
void ReadRegions(const char *FileName) {
  // Checking MPI variables validity
  if (world_size < 1 || world_rank < 0 || world_rank >= world_size) {
    FILE_LOG_SINGLE(ERROR, "'World_size' and/or 'world_rank' variable is not valid");
    TerminateFemTech(3);
  }
  FILE *fp;
  region_array = (int *)malloc(nelements * sizeof(int));
  char region[20]; 
  int i = 0;
  mel = 0, sel = 0, cel = 0, fel = 0, pel = 0, oel = 0, tel = 0;
  fp = fopen("coarse_cellcentresandvol.txt", "r");
  while (fscanf(fp, "%*d\t%*f\t%*f\t%*f\t%s\t%*f", region)!=EOF){
        for(int j = 0; j<nelements; j++)
	    if(global_eid[j]==i){
		    if(strcmp(region, "skull")==0){
			region_array[j]=1;}
		    else if(strcmp(region, "csf")==0){
			region_array[j]=2;}
		    else if(strcmp(region, "msc")==0){
			region_array[j]=3;
			mel++;}
		    else if(strcmp(region, "stem")==0){
			region_array[j]=4;
			sel++;}
		    else if(strcmp(region, "cerebellum")==0){
			region_array[j]=5;
			cel++;}
		    else if(strcmp(region, "frontal")==0){
			region_array[j]=6;
			fel++;}
		    else if(strcmp(region, "parietal")==0){
			region_array[j]=7;
			pel++;}
		    else if(strcmp(region, "occipital")==0){
			region_array[j]=8;
			oel++;}
		    else if(strcmp(region, "temporal")==0){
			region_array[j]=9;
			tel++;}
		    break;
	    }	
	i++;
  }
  fclose(fp);
}
