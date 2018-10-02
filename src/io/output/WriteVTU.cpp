#include "digitalbrain.h"


void strip_ext(char *);


void WriteVTU(char* outfile){
	FILE *fp;
	int i,j;

	strip_ext(outfile);
	strcat(outfile,".vtu");
	printf("new name: %s",outfile);
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
	for(i=0;i<nelements;i++){
		for(j=0;j<max_nnodes_per_element;j++){
			fprintf(fp,"\t\t\t\t\t%d",connectivity[max_nnodes_per_element*i+j]);
    }
    fprintf(fp,"\n");
	}
	// write offsets
	fprintf(fp,"\t\t\t\t</DataArray>\n");
	fprintf(fp,"\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
	for(i=0;i<nelements;i++){
		fprintf(fp,"\t\t\t\t\t%d\n",i*max_nnodes_per_element+8);
	}
	fprintf(fp,"\t\t\t\t</DataArray>\n");

	// write cell types
	fprintf(fp,"\t\t\t\t<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n");
	for(i=0;i<nelements;i++){
		if( strncmp(ElementType[i], "C3D8", strlen("C3D8")) == 0){
			fprintf(fp,"\t\t\t\t\t%d\n",12);
		}
		if (strncmp(ElementType[i], "C3D4", strlen("C3D4")) == 0) {
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


void strip_ext(char *fname){
    char *end = fname + strlen(fname);
    while (end > fname && *end != '.' && *end != '\\' && *end != '/') {
        --end;
    }
    if (end > fname && *end == '.') {
        *end = '\0';
    }
}
