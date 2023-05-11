#include "FemTech.h"

#include <assert.h>

void ReadMapping(const char *FileName) {
  // Checking MPI variables validity
  if (world_size < 1 || world_rank < 0 || world_rank >= world_size) {
    FILE_LOG_SINGLE(ERROR, "'World_size' and/or 'world_rank' variable is not valid");
    TerminateFemTech(3);
    }
  FILE *File;
  int id=0, hostid;
  if ((File = fopen(FileName, "rb")) == NULL) {
    FILE_LOG_SINGLE(ERROR, "Cannot open input mesh file");
    TerminateFemTech(3);
    }
  // Checking file name validity
  if (FileName == NULL || strlen(FileName) == 0) {
    FILE_LOG_SINGLE(ERROR, "Mapping file name is empty");
    TerminateFemTech(3);
    }
  while (fscanf(File, "%*d\t%d", &hostid)!=EOF){
	embedinfo[id] = hostid-1;
	nodeconstrain[connectivity[eptr[id+nelements-nembedel]]] = hostid-1;
	nodeconstrain[connectivity[eptr[id+nelements-nembedel]+1]] = hostid-1;
	id++;
	}
  }
