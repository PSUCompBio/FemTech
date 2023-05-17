#include "FemTech.h"

#include <assert.h>

void ReadMapping(const char *FileName) {
  // Checking MPI variables validity
  if (world_size < 1 || world_rank < 0 || world_rank >= world_size) {
    FILE_LOG_SINGLE(ERROR, "'World_size' and/or 'world_rank' variable is not valid");
    TerminateFemTech(3);
    }
  FILE *File;
  int id, hostid, embedelid, j=0;
  if ((File = fopen(FileName, "rb")) == NULL) {
    FILE_LOG_SINGLE(ERROR, "Cannot open input mesh file");
    TerminateFemTech(3);
    }
  // Checking file name validity
  if (FileName == NULL || strlen(FileName) == 0) {
    FILE_LOG_SINGLE(ERROR, "Mapping file name is empty");
    TerminateFemTech(3);
    }
  while (fscanf(File, "%d\t%d", &id, &hostid)!=EOF){
	for (int i=0; i<nelements; i++){
		if(id==global_eid[i])	{
			for(int k = 0; k<nelements; k++){
				if(hostid == global_eid[k]) {
					embedinfo[j] = k;
					nodeconstrain[connectivity[eptr[i]]] = k;
					nodeconstrain[connectivity[eptr[i]+1]] = k;
					j=j+1;
					break;
				}
				else { 
					nodeconstrain[connectivity[eptr[i]]] = -1;
					nodeconstrain[connectivity[eptr[i]+1]] = -1;
				}
			}
		}
	}
  }
  return;
}
