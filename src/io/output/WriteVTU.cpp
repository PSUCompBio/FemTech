#include "FemTech.h"

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
	//printf("\nwrite_VTU partition: %d\n", world_rank);
	sprintf(s, ".vtu.%04d",world_rank);
  char outfileP[ARR_SIZE] = {0};
  char outfileP2[ARR_SIZE] = {0};
  strcpy(outfileP, outfile);
  strcpy(outfileP2, outfile);
	strcat(outfile, s);
#else
	strcat(outfile, ".vtu");
#endif

	//printf("\nnew name: %s\n",outfile);
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

	// write proc ID
	fprintf(fp, "\t\t\t\t<DataArray type=\"Int32\" Name=\"ProcID\" format=\"ascii\">\n");
	for (i = 0; i < nelements; i++) {
		fprintf(fp, "\t\t\t\t\t%d\n", world_rank);
	}
	fprintf(fp, "\t\t\t\t</DataArray>\n");

	fprintf(fp, "\t\t\t</CellData>\n");
	fprintf(fp,"\t\t</Piece>\n");
	fprintf(fp,"\t</UnstructuredGrid>\n");
	fprintf(fp,"</VTKFile>\n");
	fclose(fp);

  // Write the pvtu file if you are rank zero and code in parallel
#if PARALLEL
  if (world_rank == 0) {
	  //printf("\nRank 0 Writing PVTU file\n");
	sprintf(s, ".pvtu");
	strcat(outfileP, s);
	fp=fopen(outfileP, "w");
	fprintf(fp,"<?xml version=\"1.0\"?>\n");
    fprintf(fp,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp,"\t<PUnstructuredGrid GhostLevel=\"0\">\n");
    fprintf(fp,"\t\t<PPoints>\n");
	fprintf(fp,"\t\t\t<PDataArray type=\"Float64\" NumberOfComponents=\"%d\" format=\"ascii\"/>\n",ndim);
    fprintf(fp,"\t\t</PPoints>\n");
    fprintf(fp,"\t\t<PCells>\n");
	fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>\n",ndim);
	fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\"/>\n",ndim);
	fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"types\" NumberOfComponents=\"1\"/>\n",ndim);
    fprintf(fp,"\t\t</PCells>\n");
    fprintf(fp,"\t\t<PCellData Scalars=\"PartID\">\n");
	fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"PartID\"/>\n",ndim);
	fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"ProcID\"/>\n",ndim);
    fprintf(fp,"\t\t</PCellData>\n");
    for (int i = 0; i < world_size; ++i) {
      fprintf(fp,"\t\t<Piece Source=\"%s.vtu.%.4d\"/>\n", outfileP2, i);
    }
    fprintf(fp,"\t</PUnstructuredGrid>\n");
    fprintf(fp,"</VTKFile>\n");
  }
#endif
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
