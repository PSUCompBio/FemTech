#include "digitalbrain.h"

/* Global Variables */
int nparts=0;
int nelements=0;
int nnodes=0;
int ndim=0;
double *coordinates;
int *connectivity;
int *eptr;
int *pid;
int *mid;
char **ElementType;
int MAX_ELEMENT_TYPE_SIZE = 10;

void ReadInputFile(const char* inputfile){
	char str[80];
	int i,j;
	// printf() displays the string inside quotation
	printf("Hello Digital Brain!!\n");
	FILE *fp;
	fp	= fopen(inputfile, "r");
	if (fp == NULL)
		exit(EXIT_FAILURE);

	printf("Will simulate using %s \n", inputfile);

	int flag;
	int dummyint;
	double dummyfloat;
	int counter;
	const int LENGTH_LINE = 100;
	char line[LENGTH_LINE];
	char * pch;
	int len=0;
	int nnodes_in_connectivity = 0;

/*	Parse Input file and get global vaiables */
	while (fscanf(fp,"%s",str) != EOF) {
		if(strcmp(str, "*PART")==0){
			//printf("%s\n",str);
			flag=0;
			while (flag!=1){
				fscanf(fp,"%s",str);
				//printf("i:%s\n",str);
					if (strcmp(str,"tmid")==0) {
						flag=1;
						nparts++;
					}
			}
		}
		if(strcmp(str, "*ELEMENT_SOLID")==0){
			//printf("%s\n",str);
			for(i=0;i<11;i++){
				fscanf(fp,"%s",str);
				//printf("%s\n",str);
			}
			flag=0;
			
			while ( flag != 1 ){
				fscanf(fp,"%s", str);
				//printf("%s\n",str);
				if (strcmp(str,"*NODE")==0) {
						flag=1;
				}
				else {
					nelements++;
					fgets(line, LENGTH_LINE, fp);
					printf("line: %s\n", line);
					pch = strtok(line, " ,.-");
					int j = 0;
					int k = 0;
					int dummyint_last = -1;
					while (pch != NULL){
						//printf("%s\n", pch);
						dummyint = atoi(pch);
						//printf("dummyint = %d\n", dummyint);
						if (k>0 && dummyint != dummyint_last) {
							// different vertex
							j++;
						}
						else if (k>0 && dummyint == dummyint_last) {
							//same vertex
							//printf("same vertex; vert not counted...\n");
						}
						if (k > 0) {
							dummyint_last = dummyint;
						}
						pch = strtok(NULL, " ,.-");
						k++;
					}
					printf("this element has %d different nodes\n", j);
					nnodes_in_connectivity = nnodes_in_connectivity + j;
				}
				printf("\n");
			}
		 }

		if(strcmp(str, "*NODE")==0){
			//printf("%s\n",str);
			for(i=0;i<7;i++){
				fscanf(fp,"%s",str);
				//printf("sq:%s\n",str);
			}
			ndim=3;
			flag=0;
			while( flag != 1 ){
				fscanf(fp,"%s", str);
				//printf("%s\n",str);
				if (strcmp(str,"*END")==0) {
					flag=1;
					//printf("%s\n",str);
				}
				else {
					nnodes++;
					for(i=0;i<ndim;i++){
							dummyfloat=0.0;
							fscanf(fp,"%lf",&dummyfloat);
							//printf("%f ",dummyfloat);
					}
					for(i=0;i<2;i++){
							fscanf(fp,"%d",&dummyint);
							//printf("%d ",dummyint);
					}
					//printf("\n");
				}
			}
		}
	}

	printf("ndim=%d\n",ndim);
	printf("nnodes=%d\n",nnodes);
	printf("nelements=%d\n",nelements);
	printf("nnodes_in_connect = %d\n", nnodes_in_connectivity);
	printf("nparts=%d\n",nparts);
	
	/* initalize arrays */
	coordinates = (double*)malloc((ndim*nnodes)* sizeof(double));
	pid = (int*)malloc((nelements)* sizeof(int));
	mid = (int*)malloc((nelements)* sizeof(int));
	eptr = (int*)malloc((nelements+1) * sizeof(int));
	connectivity = (int*)malloc(nnodes_in_connectivity * sizeof(int));
	ElementType = (char**)malloc((MAX_ELEMENT_TYPE_SIZE*nelements)* sizeof(char));

	//initalize coordinates
	for (i=0;i<nnodes;i++){
		for (j=0;j<ndim;j++){
			coordinates[ndim*i+j]=0.0;
		}
	}
	//initalize material, part and pointer arrays
	for (i=0;i<nelements;i++){
		mid[i]=0;
		pid[i]=0;
		eptr[i] = 0;
	}

	//initalize connectivity
	for (j = 0; j < nnodes_in_connectivity; j++) {
		connectivity[j] = 0;
	}

	/*rewind input file */
	rewind(fp);

	/* Read input and place values into arrays */
	while (fscanf(fp,"%s",str) != EOF) {
		if(strcmp(str, "*PART")==0){
			//printf("%s\n",str);
			flag=0;
			while (flag!=1){
				fscanf(fp,"%s",str);
				//printf("i:%s\n",str);
				if (strcmp(str,"tmid")==0) {
					flag=1;
				}
			}
		}
		if(strcmp(str, "*ELEMENT_SOLID")==0){
			//printf("%s\n",str);
			for(i=0;i<11;i++){
				fscanf(fp,"%s",str);
				//printf("%s\n",str);
			}
			flag=0;
			counter=0;
			int j = 0;
			int start_counter = 0;
			int end_counter = 0;
			while ( flag != 1 ){
				fscanf(fp,"%s", str);
				//printf("%s\n",str);
				if (strcmp(str,"*NODE")==0) {
						flag=1;
				}
				else {
					fgets(line, LENGTH_LINE, fp);
					printf("line: %s\n", line);
					pch = strtok(line, " ,.-");
					nnodes_in_connectivity = 0;
					int k = 0;
					int dummyint_last = -1;
					while (pch != NULL) {
						//printf("%s\n", pch);
						dummyint = atoi(pch);
						//printf("dummyint = %d\n", dummyint);
						if (k == 0) {
							pid[counter] = dummyint;
						}
						if (k > 0 && dummyint != dummyint_last) {
							// different vertex
							connectivity[j] = dummyint;
							j++;
							nnodes_in_connectivity++;
						}
						else if (k > 0 && dummyint == dummyint_last) {
							//same vertex
							//printf("same vertex; vert not counted...\n");
						}
						if (k > 0) {
							dummyint_last = dummyint;
						}
						pch = strtok(NULL, " ,.-");
						k++;
					}
					end_counter = start_counter + nnodes_in_connectivity;
					eptr[counter] = start_counter;
					eptr[counter+1] = end_counter;
					printf("this element has %d different nodes, start %d -> end %d\n", nnodes_in_connectivity, start_counter,end_counter);
					start_counter = end_counter;
				}
				//printf("\n");
				counter++;
			}
			}
		if(strcmp(str, "*NODE")==0){
			//printf("%s\n",str);
			for(i=0;i<7;i++){
				fscanf(fp,"%s",str);
				//printf("sq:%s\n",str);
			}
			flag=0;
			counter=0;
			while( flag != 1 ){
				fscanf(fp,"%s", str);
				//printf("%s\n",str);
				if (strcmp(str,"*END")==0) {
					flag=1;
					//printf("%s\n",str);
				}
				else {
					//nnodes++;
					dummyint = atoi(str);
					//printf("counter:%d, nid:%d ",counter,dummyint);
					for(i=0;i<ndim;i++){
							fscanf(fp,"%lf",&coordinates[ndim*counter+i]);
							//printf("%5.5f  ",coordinates[ndim*counter+i]);
					}
					//read tc and rc from k file
					for(i=0;i<2;i++){
							fscanf(fp,"%d",&dummyint);
							//printf("%d ",dummyint);
					}
					//printf("\n");
					counter++;
				}
			}
		}
	}

	// remove 1 from connectivity so node referencing
	// is correct
	for(i=0;i<nelements;i++){
		for(j=eptr[i];j<eptr[i+1];j++){
		 connectivity[j]--;
		}
	}

	//assign element type
	for (i = 0; i < nelements; i++) {
		if ((eptr[i + 1] - eptr[i]) == 8) {
			ElementType[i] = (char *)"C3D8";
		}
		if ((eptr[i + 1] - eptr[i]) == 4) {
			ElementType[i] = (char *)"C3D4";
		}
	}
	

	/*close input file */
	fclose(fp);
	return;
}
