#include "FemTech.h"

static void strip_ext(char *);

// Write the pvtu file if you are rank zero and code in parallel
void WritePVD(const char* FileName, int step,double time){
	static const int ARR_SIZE = 1000;
	//printf("\nRank 0 Writing PVTU file\n");
	FILE *fp;
	int i,j;
	char s[ARR_SIZE];
	char outfile[ARR_SIZE] = {0};
  char outfileP[ARR_SIZE] = {0};
 	char outfileP2[ARR_SIZE] = {0};
	if (strlen(FileName) < ARR_SIZE)
	{
			strcpy(outfile, FileName);
	}
	strip_ext(outfile);

	strcpy(outfileP, outfile);
  strcpy(outfileP2, outfile);
	sprintf(s, ".pvd");
	strcat(outfileP, s);
	fp=fopen(outfileP, "w");


	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\" >\n");
	fprintf(fp,"\t<Collection>\n");
	for(j = 0; j < step+1; j++){
		//fprintf(fp,"\t\t<Piece Source=\"%s.vtu.%.4d.%.4d\"/>\n", outfileP2,j,i);
		fprintf(fp,"\t\t<DataSet timestep=\"%d\" file=\"results/%s.%.04d.pvtu\" />\n",j,outfileP2,j);
	}
	fprintf(fp,"\t</Collection>\n");
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
