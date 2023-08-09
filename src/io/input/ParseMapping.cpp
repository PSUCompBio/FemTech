#include "FemTech.h"

#include <assert.h>

int nembedel = 0;
int *embedinfo, *embedelID, *embedproc;

void ParseMapping(const char *FileName) {
  // Checking MPI variables validity
  if (world_size < 1 || world_rank < 0 || world_rank >= world_size) {
    FILE_LOG_SINGLE(ERROR, "'World_size' and/or 'world_rank' variable is not valid");
    TerminateFemTech(3);
    }
  FILE *File;
  int id, hostid, count=0;
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
	nembedel = nembedel+1;
  }
  fclose(File);
  File = fopen(FileName, "rb");
  embedinfo = (int*)calloc(nembedel, sizeof(int));
  embedelID = (int*)calloc(nembedel, sizeof(int));
  embedproc = (int*)calloc(nembedel, sizeof(int));
  while (fscanf(File, "%d\t%d", &id, &hostid)!=EOF){
	embedinfo[count] = hostid-1;
	embedelID[count] = id-1;
	embedproc[count] = -1;
	count++;
  }
  fclose(File);
  return;
}
