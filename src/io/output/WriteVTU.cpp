#include "digitalbrain.h"

static void strip_ext(char *);

void WriteVTU(const char* FileName){
    static const int ARR_SIZE = 1000;
    
	FILE *fp;
	int i,j;
	char s[ARR_SIZE], outfile[ARR_SIZE] = {0};
	if (strlen(FileName) < ARR_SIZE)
	{
	    strcpy(outfile, FileName);
	}
	strip_ext(outfile);

#if PARALLEL
	printf("\nwrite_VTU partition: %d\n", world_rank);
	sprintf(s, ".vtu.%04d",world_rank);
	strcat(outfile, s);
#else
	strcat(outfile, ".vtu");
#endif

	printf("\nnew name: %s\n",outfile);
	fp=fopen(outfile,"w");

	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(fp,"\t<UnstructuredGrid>\n");
	fprintf(fp,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",nnodes,nelements);
	// write coordinates
	fprintf(fp,"\t\t\t<Points>\n");
	fprintf(fp,"\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"%d\" format=\"ascii\">\n",ndim);
	for(i=0;i<nnodes;i++){
		for(j=0;j<ndim;j++){
			fprintf(fp,"\t\t\t\t\t%.6e",coordinates[ndim*i+j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\t\t\t\t</DataArray>\n");
	fprintf(fp,"\t\t\t</Points>\n");
	//element connectivity
	fprintf(fp,"\t\t\t<Cells>\n");
	fprintf(fp,"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
	for (i = 0; i < nelements; i++) {
		for (j = eptr[i]; j < eptr[i + 1]; j++) {
			fprintf(fp, "\t\t\t\t\t%d", connectivity[j]);
		}
		fprintf(fp, "\n");
	}

	// write offsets
	fprintf(fp,"\t\t\t\t</DataArray>\n");
	fprintf(fp,"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
	for(i=0;i<nelements;i++){
		fprintf(fp, "\t\t\t\t\t%d\n",eptr[i+1]);
	}
	fprintf(fp,"\t\t\t\t</DataArray>\n");

	// write cell types
	fprintf(fp,"\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n");
	for(i=0;i<nelements;i++){
		if( strcmp(ElementType[i], "C3D8") == 0){
			fprintf(fp,"\t\t\t\t\t%d\n",12);
		}
		else if (strcmp(ElementType[i], "C3D4") == 0){
			fprintf(fp, "\t\t\t\t\t%d\n", 10);
		}
	}

	fprintf(fp, "\t\t\t\t</DataArray>\n");
	fprintf(fp, "\t\t\t</Cells>\n");

	// Cell Data
	fprintf(fp, "\t\t\t<CellData>\n");
	// write part ID
	fprintf(fp,"\t\t\t\t<DataArray type=\"Int32\" Name=\"PartID\" format=\"ascii\">\n");
	for(i=0;i<nelements;i++){
		 fprintf(fp,"\t\t\t\t\t%d\n",pid[i]);
	}
	fprintf(fp,"\t\t\t\t</DataArray>\n");
	

	fprintf(fp, "\t\t\t</CellData>\n");
	fprintf(fp,"\t\t</Piece>\n");
	fprintf(fp,"\t</UnstructuredGrid>\n");
	fprintf(fp,"</VTKFile>\n");
	fclose(fp);
	return;
}
//-------------------------------------------------------------------------------------------
void strip_ext(char *fname){
    char *end = fname + strlen(fname);
    while (end > fname && *end != '.' && *end != '\\' && *end != '/') {
        --end;
    }
    if (end > fname && *end == '.') {
        *end = '\0';
    }
}
//-------------------------------------------------------------------------------------------
int compare (const void * a, const void * b) {
  return ( *(int*)a - *(int*)b );
}
//-------------------------------------------------------------------------------------------
int unique(int arr[], int n) {
	//int temp[n];
	int *temp;
	int j = 0;
	temp = (int *)calloc(n, sizeof(int));
    for (int i=0; i < n-1; i++) {
      if (arr[i] != arr[i+1]) {
        temp[j++] = arr[i];
      }
    }
    temp[j++] = arr[n-1];
    for (int i = 0; i < j; i++) {
      arr[i] = temp[i];
    }
	free(temp);
    return j;
}
//-------------------------------------------------------------------------------------------
void updateConnectivityGlobalToLocal(void) {
  int totalSize = eptr[nelements]-eptr[0];
  int *newConnectivity = (int*)malloc(totalSize*sizeof(int));   
  int *sorted = (int*)malloc(totalSize*sizeof(int));   
  memcpy(sorted, connectivity, totalSize*sizeof(int));
  qsort(sorted, totalSize, sizeof(int), compare);
  nCoordinates = unique(sorted, totalSize);
  for(int i = 0; i < totalSize; ++i) {
    int j;
    for (j = 0; j < nCoordinates; ++j) {
      if (sorted[j] == connectivity[i]) {
        break;
      }
    }
    newConnectivity[i] = j;
  }
  // Reoder co-ordinates 
  double *newCoordinates = (double*)malloc(ndim*nCoordinates*sizeof(double));
  for (int j = 0; j < nCoordinates; ++j) {
    int i;
    for(i = 0; i < totalSize; ++i) {
      if (sorted[j] == connectivity[i]) {
        break;
      }
    }
    memcpy(newCoordinates+ndim*j, coordinates+ndim*i, ndim*sizeof(double));
  }
  memcpy(connectivity, newConnectivity, totalSize*sizeof(int));
  memcpy(coordinates, newCoordinates, ndim*nCoordinates*sizeof(double));
  free(newCoordinates);
  free(newConnectivity);
  free(sorted);
}
