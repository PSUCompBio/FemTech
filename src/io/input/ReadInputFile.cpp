#include "digitalbrain.h"

/* Globla Variables */
int nparts=0;
int nelements=0;
int nnodes=0;
int ndim=0;
int max_nnodes_per_element=0;
double *coordinates;
int *connectivity;
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

/*	Parse Input file and get global vaiables */
	while (fscanf(fp,"%s",&str) != EOF) {
		if(strcmp(str, "*PART")==0){
			//printf("%s\n",str);
			flag=0;
			while (flag!=1){
				fscanf(fp,"%s",&str);
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
				fscanf(fp,"%s",&str);
				//printf("%s\n",str);
			}
			max_nnodes_per_element = 8;
			flag=0;
			while ( flag != 1 ){
				fscanf(fp,"%s", &str);
				//printf("%s\n",str);
				if (strcmp(str,"*NODE")==0) {
						flag=1;
				}
				else {
					nelements++;
					for(i=0;i<9;i++){
							dummyint=0;
							fscanf(fp,"%d ",&dummyint);
							//printf("%d ",dummyint);
					}
				}
				//printf("\n");
			}
		 }
		if(strcmp(str, "*NODE")==0){
			//printf("%s\n",str);
			for(i=0;i<7;i++){
				fscanf(fp,"%s",&str);
				//printf("sq:%s\n",str);
			}
			ndim=3;
			flag=0;
			while( flag != 1 ){
				fscanf(fp,"%s", &str);
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
	printf("nparts=%d\n",nparts);

	/* initalize arrays */
	coordinates = (double*)malloc((ndim*nnodes)* sizeof(double));
	pid = (int*)malloc((nelements)* sizeof(int));
	mid = (int*)malloc((nelements)* sizeof(int));
	connectivity = (int*)malloc((max_nnodes_per_element*nelements)* sizeof(int));
	ElementType = (char**)malloc((MAX_ELEMENT_TYPE_SIZE*nelements)* sizeof(char));

	//initalize coordinates
	for (i=0;i<nnodes;i++){
		for (j=0;j<ndim;j++){
			coordinates[ndim*i+j]=0.0;
		}
	}
	//initalize connectivity
	for (i=0;i<nelements;i++){
		mid[i]=0;
		pid[i]=0;
		for (j=0;j<max_nnodes_per_element;j++){
			connectivity[max_nnodes_per_element*i+j]=0;
		}
	}

	/*rewind input file */
	rewind(fp);

	/* Read input and place values into arrays */
	while (fscanf(fp,"%s",&str) != EOF) {
		if(strcmp(str, "*PART")==0){
			//printf("%s\n",str);
			flag=0;
			while (flag!=1){
				fscanf(fp,"%s",&str);
				//printf("i:%s\n",str);
				if (strcmp(str,"tmid")==0) {
					flag=1;
				}
			}
		}
		if(strcmp(str, "*ELEMENT_SOLID")==0){
			//printf("%s\n",str);
			for(i=0;i<11;i++){
				fscanf(fp,"%s",&str);
				//printf("%s\n",str);
			}
			flag=0;
			counter=0;
			while ( flag != 1 ){
				fscanf(fp,"%s", &str);
				//printf("%s\n",str);
				if (strcmp(str,"*NODE")==0) {
						flag=1;
				}
				else {
					// read eid
					dummyint=0;
					//fscanf(fp,"%d",&dummyint);
					dummyint = atoi(str);
					//printf("counter:%d, e:%d ",counter,dummyint);
					// read pid
					fscanf(fp,"%d",&pid[counter]);
					//printf("p:%d ",pid[counter]);
					// read nodes
					for(i=0;i<max_nnodes_per_element;i++){
							fscanf(fp,"%d",&connectivity[max_nnodes_per_element*counter+i]);
							//printf("n%d:%d ",i+1,connectivity[max_nnodes_per_element*counter+i]);
					}

				}
				//printf("\n");
				counter++;
			}
			}
		if(strcmp(str, "*NODE")==0){
			//printf("%s\n",str);
			for(i=0;i<7;i++){
				fscanf(fp,"%s",&str);
				//printf("sq:%s\n",str);
			}
			flag=0;
			counter=0;
			while( flag != 1 ){
				fscanf(fp,"%s", &str);
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
		for(j=0;j<max_nnodes_per_element;j++){
		 connectivity[max_nnodes_per_element*i+j]--;
		}
	}

	//assign element type
	int nodecounter = 0;
	for(i=0;i<nelements;i++){
		nodecounter = 1;
		//printf("number of nodes in this element %d: \n",i);
		for (j = 1; j < max_nnodes_per_element; j++) {
			//printf("%d %d\n", connectivity[max_nnodes_per_element*i + (j - 1)], connectivity[max_nnodes_per_element*i + j]);
			if (connectivity[max_nnodes_per_element*i + (j-1)] == connectivity[max_nnodes_per_element*i + j]) {
				//printf("duplicate\n");
			}
			else {
				nodecounter++;
			}
		}
		//printf("nodecounter = %d\n",nodecounter);
		if (nodecounter == 8) {
			ElementType[i] = "C3D8";
		}
		if (nodecounter == 4) {
			ElementType[i] = "C3D4";
		}
		//printf("%s \n",ElementType[i]);
	}

	/*close input file */
	fclose(fp);
	return;
}
