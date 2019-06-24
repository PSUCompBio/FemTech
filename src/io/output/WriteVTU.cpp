#include "FemTech.h"

static void strip_ext(char *);

void WriteVTU(const char* FileName, int step,double time){
    static const int ARR_SIZE = 1000;

	FILE *fp;
	int i,j;
	char s[ARR_SIZE];
  char outfile[ARR_SIZE] = {0};
	char paths[ARR_SIZE] = {0};
	char paths2[ARR_SIZE] = {0};
	char outfileP[ARR_SIZE] = {0};
  char outfileP2[ARR_SIZE] = {0};

	if (strlen(FileName) < ARR_SIZE)
	{
	    strcpy(outfile, FileName);
	}
	strip_ext(outfile);

	//printf("\nwrite_VTU partition: %d\n", world_rank);
	sprintf(s, ".vtu.%04d.%04d",step,world_rank);
  sprintf(paths,"./results/vtu/");
	sprintf(paths2,"./results/");
  strcpy(outfileP, outfile);
  strcpy(outfileP2, outfile);
 	strcat(paths2,outfileP2);
  strcat(paths,outfile);
  strcat(paths,s);


	//printf("\nnew name: %s\n",outfile);
	fp=fopen(paths,"w");

	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(fp,"\t<UnstructuredGrid>\n");
	fprintf(fp,"\t\t<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",nnodes,nelements);
	// write coordinates

	// Temporary soution for ndim
	//const int NDim = ndim == 2 ? 3 : ndim;
	//

	fprintf(fp,"\t\t\t<Points>\n");
	fprintf(fp,"\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"%d\" format=\"ascii\">\n", ndim);
	for(i=0;i<nnodes;i++){
		fprintf(fp,"\t\t\t\t\t");
		for(j=0;j<ndim;j++){
			if(fabs(displacements[ndim*i+j]-0.0)< 1e-20){
				displacements[ndim*i+j]=0.0;
			}
			fprintf(fp,"%10.8e ",coordinates[ndim*i+j] + displacements[ndim*i+j]);
		}
		// Temporary solution for ndim
		if(ndim == 2){fprintf(fp,"%10.8e", 0.0);}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\t\t\t\t</DataArray>\n");
	fprintf(fp,"\t\t\t</Points>\n");

	//element connectivity
	fprintf(fp,"\t\t\t<Cells>\n");
	fprintf(fp,"\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
	for (i = 0; i < nelements; i++) {
		fprintf(fp, "\t\t\t\t\t");
		for (j = eptr[i]; j < eptr[i + 1]; j++) {
			fprintf(fp, "%d ", connectivity[j]);
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
		else if (strcmp(ElementType[i], "T3D2") == 0){
			fprintf(fp, "\t\t\t\t\t%d\n", 3);
		}
	}
	fprintf(fp, "\t\t\t\t</DataArray>\n");
	fprintf(fp, "\t\t\t</Cells>\n");



	/*-----------POINT DATA -----------------*/
	fprintf(fp, "\t\t\t<PointData>\n");
	// write displacements
	fprintf(fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"Displacements\" "
									"NumberOfComponents=\"%d\" ComponentName0=\"X\" "
									"ComponentName1=\"Y\" ComponentName2=\"Z\" "
									"format=\"ascii\">\n",ndim);
	for(i=0;i<nnodes;i++){
			fprintf(fp,"\t\t\t\t\t");
			for(j=0;j<ndim;j++){
				fprintf(fp,"%10.8e ",displacements[ndim*i+j]);
			}
			// Temporary solution for ndim
			if(ndim == 2){fprintf(fp,"%10.8e", 0.0);}
			fprintf(fp,"\n");
	}
	fprintf(fp,"\t\t\t\t</DataArray>\n");

  // write acceleration
  fprintf(fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"Accelerations\" "
                  "NumberOfComponents=\"%d\" ComponentName0=\"X\" "
                  "ComponentName1=\"Y\" ComponentName2=\"Z\" "
                  "format=\"ascii\">\n",ndim);
  for(i=0;i<nnodes;i++){
      fprintf(fp,"\t\t\t\t\t");
      for(j=0;j<ndim;j++){
        fprintf(fp,"%10.8e ",accelerations[ndim*i+j]);
      }
      // Temporary solution for ndim
      if(ndim == 2){fprintf(fp,"%10.8e", 0.0);}
      fprintf(fp,"\n");
  }
  fprintf(fp,"\t\t\t\t</DataArray>\n");

	// write boundary
	fprintf(fp,"\t\t\t\t<DataArray type=\"Int32\" Name=\"Boundary\" "
									"NumberOfComponents=\"%d\" ComponentName0=\"X\" "
									"ComponentName1=\"Y\" ComponentName2=\"Z\" "
									"format=\"ascii\">\n",ndim);
	for(i=0;i<nnodes;i++){
			fprintf(fp,"\t\t\t\t\t");
			for(j=0;j<ndim;j++){
				fprintf(fp,"%d ",boundary[ndim*i+j]);
			}
			// Temporary solution for ndim
			if(ndim == 1){fprintf(fp,"0.0  0.0");}
			if(ndim == 2){fprintf(fp," 0.0");}

			fprintf(fp,"\n");
	}
	fprintf(fp,"\t\t\t\t</DataArray>\n");


	fprintf(fp, "\t\t\t</PointData>\n");

	/*-----------CELL DATA -----------------*/
	fprintf(fp, "\t\t\t<CellData>\n");
	// write part ID
	fprintf(fp,"\t\t\t\t<DataArray type=\"Int32\" Name=\"PartID\" format=\"ascii\">\n");
	for(i=0;i<nelements;i++){
		 fprintf(fp,"\t\t\t\t\t%d\n",pid[i]);
	}
	fprintf(fp,"\t\t\t\t</DataArray>\n");

  fprintf(fp,"\t\t\t\t<DataArray type=\"Float64\" Name=\"AvgStrain\" "
                  "NumberOfComponents=\"%d\" ComponentName0=\"E11\" "
                  "ComponentName1=\"E12\" ComponentName2=\"E13\" "
                  "ComponentName3=\"E12\" ComponentName4=\"E13\" "
                  "ComponentName5=\"E12\" ComponentName6=\"E13\" "
                  "ComponentName7=\"E12\" ComponentName8=\"E33\" "
                  "format=\"ascii\">\n",ndim*ndim);
  for(i=0;i<nelements;i++){
      fprintf(fp,"\t\t\t\t\t");
        for(j=0;j<ndim*ndim;j++){
           fprintf(fp,"%10.8e ",Eavg[i*ndim*ndim+j]);
      }
      // Temporary solution for ndim
      if(ndim == 2){fprintf(fp,"%10.8e", 0.0);}
      fprintf(fp,"\n");
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
  if (world_rank == 0) {
	  //printf("\nRank 0 Writing PVTU file\n");
		sprintf(s, ".%.04d.pvtu",step);
		strcat(paths2, s);
		fp=fopen(paths2, "w");
		fprintf(fp,"<?xml version=\"1.0\"?>\n");
    fprintf(fp,"<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp,"\t<PUnstructuredGrid GhostLevel=\"0\">\n");
    fprintf(fp,"\t\t<PPoints>\n");
		fprintf(fp,"\t\t\t<PDataArray type=\"Float64\" NumberOfComponents=\"%d\" format=\"ascii\"/>\n",ndim);
    fprintf(fp,"\t\t</PPoints>\n");
    fprintf(fp,"\t\t<PCells>\n");
		fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>\n");
		fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\"/>\n");
		fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"types\" NumberOfComponents=\"1\"/>\n");
    fprintf(fp,"\t\t</PCells>\n");

		/*-----------POINT DATA -----------------*/
		fprintf(fp,"\t\t<PPointData  Vectors=\"Displacements Accelerations\" >\n");
		fprintf(fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Displacements\" "
               "NumberOfComponents=\"%d\" ComponentName0=\"X\" ComponentName1=\"Y\" "
               "ComponentName2=\"Z\" format=\"ascii\" />\n",ndim);
   fprintf(fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"Accelerations\" "
              "NumberOfComponents=\"%d\" ComponentName0=\"X\" ComponentName1=\"Y\" "
              "ComponentName2=\"Z\" format=\"ascii\" />\n",ndim);
		fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"Boundary\" "
							 "NumberOfComponents=\"%d\" ComponentName0=\"X\" ComponentName1=\"Y\" "
							 "ComponentName2=\"Z\" format=\"ascii\" />\n",ndim);
    fprintf(fp,"\t\t</PPointData>\n");

		/*-----------CELL DATA -----------------*/
  	fprintf(fp,"\t\t<PCellData Scalars=\"PartID\" Tensors=\"AvgStrain\">\n");
		fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"PartID\"/>\n");
    fprintf(fp,"\t\t\t<PDataArray type=\"Float64\" Name=\"AvgStrain\" "
                    "NumberOfComponents=\"%d\" ComponentName0=\"E11\" "
                    "ComponentName1=\"E12\" ComponentName2=\"E13\" "
                    "ComponentName3=\"E12\" ComponentName4=\"E13\" "
                    "ComponentName5=\"E12\" ComponentName6=\"E13\" "
                    "ComponentName7=\"E12\" ComponentName8=\"E33\" "
                    "format=\"ascii\"/>\n",ndim*ndim);
		fprintf(fp,"\t\t\t<PDataArray type=\"Int32\" Name=\"ProcID\"/>\n");
    fprintf(fp,"\t\t</PCellData>\n");
    for (int i = 0; i < world_size; ++i) {
      fprintf(fp,"\t\t<Piece Source=\"vtu/%s.vtu.%04d.%.4d\"/>\n", outfileP2,step, i);
    }
    fprintf(fp,"\t</PUnstructuredGrid>\n");
    fprintf(fp,"</VTKFile>\n");
		fclose(fp);
  }
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
