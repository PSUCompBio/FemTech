#include "digitalbrain.h"

/* Global Variables */
int nparts=0;
int nelements=0;        // number of elements in processor
int nallelements=0;     // number of all elements in mesh
int nnodes=0;
int ndim=0;
double *coordinates;
int *connectivity;
int *eptr;
int *pid;
int *mid;
char **ElementType;
int world_rank;
int world_size;


bool ReadInputFile(const char *FileName){
    // Checking MPI variables validity
    if (world_size < 1 || world_rank < 0 || world_rank >= world_size) {
        printf("\nERROR( proc %d ): 'world_size' and/or 'world_rank' variable is not valid.\n", world_rank);
        return false;
    }

    // Checking if mesh file can be opened or not
    FILE *File;
    if (FileName == NULL || strlen(FileName) == 0 || (File = fopen(FileName, "rt")) == NULL) {
        printf("\nERROR( proc %d ): Failed to open input file.\n", world_rank);
        return false;
    }

    // Counting number of elements in elements section of mesh file
    // And fixing elements and nodes sections' positions in the file
    char Line[MAX_FILE_LINE];
    long ElementsSectionPos = 0xFFFFFFFF, NodesSectionPos = 0xFFFFFFFF;
    while (fgets(Line, sizeof(Line), File) != NULL) {
        if (strcmp(Line, "*ELEMENT_SOLID\n") == 0) {
            ElementsSectionPos = ftell(File);
            while (fgets(Line, sizeof(Line), File) != NULL) {
                if (LineToArray(true, false, 2, 1, Line, NULL) > 0 && LineToArray(true, true, 3, 0, Line, NULL) > 0) {
                    nallelements++;
                }
                if (strcmp(Line, "*NODE\n") == 0) {
                    NodesSectionPos = ftell(File);
                    break;
                }
            }
            break;
        }
    }

    // Checking if elements were found in mesh file
    // And checking if we can be back to elements section in the file
    if (nallelements == 0 || fseek(File, ElementsSectionPos, SEEK_SET) != 0) {
        fclose(File);
        if (nallelements == 0) {
            printf("\nERROR( proc %d ): No elemens found. This means input file is empty or contains invalid data.\n", world_rank);
        }
        else {
            printf("\nERROR( proc %d ): 'fseek()' call failed.\n", world_rank);
        }
        return false;
    }
    
    // Determining number/count of elements processor will hold
    nelements = nallelements / world_size;
    if (world_rank == world_size - 1) {
        nelements += (nallelements % world_size);
    }

    // Fixing range of lines (in elements section of mesh file) from which processor will read its own elements
    const int From = world_rank * (nallelements / world_size);
    const int To = From + nelements - 1;
    
    // Creating and initializing "eptr" array, and determining size of "connectivity" array of processor
    eptr = (int *)calloc(nelements + 1, sizeof(int));
    int i = 0, j = 0, ConnectivitySize = 0;
    while (fgets(Line, sizeof(Line), File) != NULL) {
        const int np = LineToArray(true, false, 2, 1, Line, NULL);
        const int nn = LineToArray(true, true, 3, 0, Line, NULL);
        if (np == 0 || nn == 0) {
            continue;
        }
        else if (i >= From && i <= To) {
            j = j + 1;
            eptr[j] = eptr[j - 1] + nn;
            ConnectivitySize = ConnectivitySize + nn;
        }
        if (strcmp(Line, "*NODE\n") == 0) {
            break;
        }
        i = i + 1;
    }

    // Checking if calculated size of connectivity array is valid
    // And checking if we can be back to elements section of mesh file
    if (ConnectivitySize != eptr[nelements] || fseek(File, ElementsSectionPos, SEEK_SET) != 0) {
        fclose(File);
        FreeArrays();
        if (ConnectivitySize != eptr[nelements]) {
            printf("\nERROR( proc %d ): Size of 'eind' array is not vald.\n", world_rank);
        }
        else {
            printf("\nERROR( proc %d ): 'fseek()' call failed.\n", world_rank);
        }
        return false;
    }

    // Creating and initializing "connectivity", "pid" and "ElementType" arrays for processor
    connectivity = (int *)calloc(ConnectivitySize, sizeof(int));
    pid = (int *)calloc(nelements, sizeof(int));
    ElementType = (char **)calloc(nelements, sizeof(char *));
    for (int i = 0; i < nelements; i++) {
        ElementType[i] = (char *)calloc(MAX_ELEMENT_TYPE_SIZE, sizeof(char));
    }
    i = 0;
    int ci = 0, pi = 0;
    while (fgets(Line, sizeof(Line), File) != NULL) {
        int *PIDs, *Nodes;
        const int np = LineToArray(true, false, 2, 1, Line, (void**)&PIDs);
        const int nn = LineToArray(true, true, 3, 0, Line, (void**)&Nodes);
        if (np == 0 || nn == 0) {
            if (np > 0) {
                free(PIDs);
            }
            if (nn > 0) {
                free(Nodes);
            }
            continue;
        }
        else if (i >= From && i <= To) {
            pid[pi] = PIDs[0];
            
            if (nn == 8) {
                strcpy(ElementType[pi], "C3D8");
            }
            else if (nn == 4) {
                strcpy(ElementType[pi], "C3D4");
            }
            pi = pi + 1;
             
            for (int j = 0; j < nn; j++) {
                connectivity[ci] = Nodes[j] - 1;
                ci = ci + 1;
            }
        }
        free(PIDs);
        free(Nodes);
        if (strcmp(Line, "*NODE\n") == 0) {
            break;
        }
        i = i + 1;
    }

    // Checking if we can go to nodes section of mesh file
    if (fseek(File, NodesSectionPos, SEEK_SET) != 0) {
        fclose(File);
        FreeArrays();
        printf("\nERROR( proc %d ): 'fseek()' call failed.\n", world_rank);
        return false;
    }
    
    // "nnodes" is not unique node numbers for now
    // "coordinates" array will contain "ndim" coordinates for each element of "connectivity" array
    nnodes = 0;
    ndim = 0;
    coordinates = NULL;
    
    // Initializing "coordinates" array
    while (fgets(Line, sizeof(Line), File) != NULL) {
        if (ndim == 0 && strstr(Line, "nid") != NULL) {
            if (strstr(Line, "x") != NULL) {
                ndim = ndim + 1;
            }
            if (strstr(Line, "y") != NULL) {
                ndim = ndim + 1;
            }
            if (strstr(Line, "z") != NULL) {
                ndim = ndim + 1;
            }
            if (ndim == 2 || ndim == 3) {
                coordinates = (double *)calloc(ConnectivitySize * ndim, sizeof(double));
            }
            continue;
        }
        
        if (coordinates != NULL) {
            double *NodeData;
            if (LineToArray(false, false, 1, ndim + 1, Line, (void**)&NodeData) == (ndim + 1)) {
                const int NodeID = int(NodeData[0]) - 1;
                // Searching node id in "connectivity" array
                for (int i = 0; i < ConnectivitySize; i++) {
                    if (connectivity[i] == NodeID) {
                        // Node found. Getting coordinates of it
                        memcpy(&coordinates[i * ndim], &NodeData[1], ndim * sizeof(double));
                        nnodes = nnodes + 1;
                    }
                }
                free(NodeData);
            }
        }
        
        if (strcmp(Line, "*END\n") == 0) {
            break;
        }
    }
    
    // Checking if "coordinates" array is OK
    if (coordinates == NULL || nnodes != ConnectivitySize) {
        fclose(File);
        FreeArrays();
        printf("\nERROR( proc %d ): Failed to initialize 'coordinates' array.\n", world_rank);
        if (coordinates == NULL) {
            printf("\ncoordinates == NULL\n");
        }
        printf("\nnnodes = %d, ConnectivitySize = %d, ndim = %d\n", nnodes, ConnectivitySize, ndim);
        return false;
    }
    
    // Printing local arrays of processor (this section can be removed)
    printf("\neptr array in processor %d before partitioning = ", world_rank);
    for (int i = 0; i <= nelements ; i++) {
        printf("%d ", eptr[i]);
    }
    printf("\n");
    printf("\neind array in processor %d before partitioning =", world_rank);
    for (int i = 0; i < nelements; i++) {
        printf(" (%d)  ", i);
        for (int j = eptr[i]; j < eptr[i + 1]; j++) {
            printf("%d ", connectivity[j]);
        }
    }
    printf("\n");
    printf("\nType/PartID of element in processor %d before partitioning = ", world_rank);
    for (int i = 0; i < nelements ; i++) {
        printf("%s/%d  ", ElementType[i], pid[i]);
    }
    printf("\n");
    printf("\nSize of coordinates array in processor %d before partitioning = %d\n", world_rank, nnodes * ndim);
    printf("\nCoordinates array in processor %d before partitioning =", world_rank);
    for (int i = 0; i < nnodes; i++) {
        printf(" (%d)  ", i);
        for (int j = 0; j < ndim; j++) {
            printf("%.*f ", 1, coordinates[ndim * i + j]);
        }
    }
    printf("\n");

    // Everything is OK for now. Closing the mesh file, returning TRUE.
    fclose(File);
    return true;
}
//-------------------------------------------------------------------------------------------
/*
   Converts line into array of int or double if line contains
   valid int/double separated by whitespase or tab.
*/
int LineToArray(
    const bool IntOrFloat,      // If true, will be converted to int else to double
    const bool CheckLastVal,    // If true, if next column's value is the same as the last column's one, next coulmns won't be fetched
    const int ColumnToStart,    // One-based index of column to start to fetch
    const int ColumnCount,      // Number of columns to fetch. Will be used if greater than zero
    const char *ConstLine,      // Line to convert
    void **Array) {             // Output array, int or double

    int Result = 0;
    char Line[MAX_FILE_LINE];
    strcpy(Line, ConstLine);
    double LastVal = INT_MAX;
    int Column = 0;
    char *T = strtok(Line, " \t");
    while (T != NULL) {
        char *EndPtr;
        const double Val = strtod(T, &EndPtr);
        if (EndPtr == T) {
            // Error. Token in a line cannot be converted to number.
            return 0;
        }
        T = strtok(NULL, " \t");
        Column = Column + 1;
        if (Column < ColumnToStart) {
            // Skip columns to go to required ones.
            continue;
        }
        if ((CheckLastVal && Val == LastVal) || (ColumnCount > 0 && Result >= ColumnCount)) {
            // Stop fetching columns.
            break;
        }
        LastVal = Val;
        Result = Result + 1;
    }

    if (Result == 0 || Array == NULL) {
        return Result;
    }

    strcpy(Line, ConstLine);
    int *IntArray;
    double *DoubleArray;
    if (IntOrFloat) {
        IntArray = (int *)calloc(Result, sizeof(int));
        *Array = (void *)IntArray;
    }
    else {
        DoubleArray = (double *)calloc(Result, sizeof(double));
        *Array = (void *)DoubleArray;
    }
    
    int i = 0;
    Column = 0;
    T = strtok(Line, " \t");
    while (T != NULL && i < Result) {
        const double Val = strtod(T, NULL);
        T = strtok(NULL, " \t");
        Column = Column + 1;
        if (Column < ColumnToStart) {
            continue;
        }
        if (IntOrFloat) {
            IntArray[i] = int(Val);
        }
        else {
            DoubleArray[i] = Val;
        }
        i = i + 1;
    }
    return Result;
}
