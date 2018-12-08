#include "digitalbrain.h"


static void strip_ext(char *);
//-------------------------------------------------------------------------------------------
bool PrepareVTUData(const char *FileName, const int PartIdx, const int NParts, const PARTITION Partition)
{
    // Global variables that change for each partition
    // These wouldn't be used if WriteUTV() function is called inside CreatePartitions() function
    coordinates = NULL;
    connectivity = NULL;
    eptr = NULL;
    pid = NULL;
    ElementType = NULL;
    nelements = 0;
    nnodes = 0;
    
    static NODE *Nodes = NULL;
    static int NodeIDsCount;
    
    if (Nodes == NULL)
    {
        FILE *File = fopen(FileName, "rt");
        if (File == NULL)
        {
            printf("\nERROR( proc %d ): Couldn't open input file.\n", world_rank);
            return false;
        }

        ndim = 0;  // Global variable
        
        NodeIDsCount = 0;
        char Line[MAX_FILE_LINE];
        long NodesSectionPos = 0;
        while (fgets(Line, sizeof(Line), File) != NULL)
        {
            if (strcmp(Line, "*NODE\n") == 0)
            {
                NodesSectionPos = ftell(File);
                while (fgets(Line, sizeof(Line), File) != NULL)
                {
                    if (ndim == 0 && strstr(Line, "nid") != NULL)
                    {
                        if (strstr(Line, "x") != NULL)
                        {
                            ndim++;
                        }
                        if (strstr(Line, "y") != NULL)
                        {
                            ndim++;
                        }
                        if (strstr(Line, "z") != NULL)
                        {
                            ndim++;
                        }
                        continue;
                    }
                    if ((ndim == 2 || ndim == 3) && LineToArray(false, false, 1, ndim + 1, Line, NULL) == (ndim + 1))
                    {
                        NodeIDsCount++;
                    }
                    if (strcmp(Line, "*END\n") == 0)
                    {
                        break;
                    }
                }
                break;
            }
        }

        if (NodeIDsCount == 0 || fseek(File, NodesSectionPos, SEEK_SET) != 0)
        {
            fclose(File);
            if (NodeIDsCount == 0)
            {
                printf("\nERROR( proc %d ): No nodes found. This means input file is empty or contains invalid data.\n", world_rank);
            }
            else
            {
                printf("\nERROR( proc %d ): 'fseek()' call failed.\n", world_rank);
            }
            return false;
        }

        int i = 0;
        Nodes = (NODE *)calloc(NodeIDsCount, sizeof(NODE));
        while (fgets(Line, sizeof(Line), File) != NULL)
        {
            double *NodeData;
            const int n = LineToArray(false, false, 1, ndim + 1, Line, (void**)&NodeData);
            if (n == 0)
            {
                continue;
            }
            else if (n == (ndim + 1) && i < NodeIDsCount)
            {
                Nodes[i].Id = NodeData[0];
                for (int j = 0; j < ndim; j++)
                {
                    Nodes[i].XYZ[j] = NodeData[j + 1];
                }
            }
            free(NodeData);
            if (strcmp(Line, "*END\n") == 0)
            {
                break;
            }
            i++;
        }
        
        fclose(File);
    }
    
    nelements = Partition.Count;
    for (int i = 0; i < Partition.Count; i++)
    {
        nnodes += Partition.Elements[i].Count;
    }
    
    if (nnodes == 0 || nelements == 0)
    {
        if (Nodes != NULL)
        {
            free(Nodes);
        }
        printf("\nERROR( proc %d ): Failed to initialize 'nnodes' and/or 'nelements' variables.\n", world_rank);
        return false;
    }
    
    const int CoordSize = nnodes * ndim;
    coordinates = (double *)calloc(CoordSize, sizeof(double));
    connectivity = (int *)calloc(nnodes, sizeof(int));
    eptr = (int *)calloc(nelements + 1, sizeof(int));
    pid = (int *)calloc(nelements, sizeof(int));
    ElementType = (char **)calloc(nelements, sizeof(char *));
    for (int i = 0; i < nelements; i++)
    {
        ElementType[i] = (char *)calloc(MAX_ELEMENT_TYPE_SIZE, sizeof(char));
    }
    
    memcpy(eptr, Partition.Eptr, (nelements + 1) * sizeof(int));
    
    for (int i = 0, ci = 0, ni = 0; i < nelements; i++)
    {
        pid[i] = Partition.Elements[i].PId;
        
        if (Partition.Elements[i].Count == 8)
        {
            strcpy(ElementType[i], "C3D8");
        }
        else if (Partition.Elements[i].Count == 4)
        {
            strcpy(ElementType[i], "C3D4");
        }
        
        for (int j = 0; j < Partition.Elements[i].Count; j++)
        {
            connectivity[ni++] = Partition.Elements[i].Nodes[j];

            for (int k = 0; k < NodeIDsCount; k++)
            {
                if (Partition.Elements[i].Nodes[j] == Nodes[k].Id)
                {
                    for (int x = 0; x < ndim; x++)
                    {
                        if ((ci + x) < CoordSize)
                        {
                            coordinates[ci + x] = Nodes[k].XYZ[x];
                        }
                    }
                    ci += ndim;
                    break;
                }
            }
        }
    }
        
    if (Nodes != NULL && PartIdx == (NParts - 1))
    {
        free(Nodes);
    }
    
    return true;
}
//-------------------------------------------------------------------------------------------
void WriteVTU(const char* FileName, const int PartIdx){
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
	printf("\nwrite_VTU partition: %d\n", PartIdx);
	sprintf(s, ".vtu.%04d", PartIdx);
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
