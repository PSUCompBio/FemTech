#include "digitalbrain.h"

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
   //fopen(fp,"r");
   //InitializeGlobalVariables();
   //return 0;

	int flag;
	int dummyint;
	double dummyfloat;

	int nparts=0;
	int nelements=0;
	int nnodes=0;
	int ndim=3;

	while (fscanf(fp,"%s",&str) != EOF) {
		if(strcmp(str, "*PART")==0){
				printf("%s\n",str);
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
				printf("%s\n",str);
				for(i=0;i<11;i++){
					fscanf(fp,"%s",&str);
					//printf("%s\n",str);
				}
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
				printf("%s\n",str);
				for(i=0;i<7;i++){
				fscanf(fp,"%s",&str);
				//printf("sq:%s\n",str);
			}
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
	printf("nparts=%d\n",nparts);
	printf("nelements=%d\n",nelements);
	printf("nnodes=%d\n",nnodes);
  fclose(fp);
  return;
}
