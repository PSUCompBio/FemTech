#include "FemTech.h"

#include <assert.h>

#define TITLE_ELEMENT_SOLID  "*ELEMENT_SOLID"
#define TITLE_ELEMENT_SHELL  "*ELEMENT_SHELL"
#define TITLE_ELEMENT_BEAM   "*ELEMENT_BEAM"
#define TITLE_NODE           "*NODE"
#define TITLE_END            "*END"

//-------------------------------------------------------------------------------------------
static bool IsElementSection(char *Line, FILE *File, int *NDimExceptRuleExists = NULL) {
    const bool ElementSection =
        strncmp(Line, TITLE_ELEMENT_SOLID, strlen(TITLE_ELEMENT_SOLID)) == 0 ||
        strncmp(Line, TITLE_ELEMENT_SHELL, strlen(TITLE_ELEMENT_SHELL)) == 0 ||
        strncmp(Line, TITLE_ELEMENT_BEAM, strlen(TITLE_ELEMENT_BEAM)) == 0;

    if (ElementSection) {
        if (NDimExceptRuleExists != NULL && strncmp(Line, TITLE_ELEMENT_BEAM, strlen(TITLE_ELEMENT_BEAM)) == 0) {
            *NDimExceptRuleExists = 1;
        }
        char* read = fgets(Line, MAX_FILE_LINE * sizeof(char), File);
        assert(read != NULL);
    }
    return ElementSection;
}
//-------------------------------------------------------------------------------------------
static bool IsNodeSection(char *Line) {
    return strncmp(Line, TITLE_NODE, strlen(TITLE_NODE)) == 0;
}
//-------------------------------------------------------------------------------------------
static int CalcUniqueValues(const int *Array, const int Size) {
    int Result = 1;
    for (int i = 1; i < Size; i++) {
        int j = 0;
        for (j = 0; j < i; j++) {
            if (Array[i] == Array[j]) {
                break;
            }
        }
        if (i == j) {
            Result++;
        }
    }
    return Result;
}
//-------------------------------------------------------------------------------------------
bool ReadLsDyna(const char *FileName) {
    // Checking if mesh file can be opened or not
    FILE *File;
    if ((File = fopen(FileName, "rb")) == NULL) {
        FILE_LOG_SINGLE(ERROR, "Cannot open input mesh file");
        return false;
    }

    // Counting number of elements in elements section of mesh file
    // And fixing elements and nodes sections' positions in the file
    const char *Delim = " \t";
    char Line[MAX_FILE_LINE];
    long ElementsSectionPos = -1, NodesSectionPos = -1, CurrentPos = ftell(File);;
    int NDimExceptRuleExists = 0;
    ndim = 0;
    while (fgets(Line, sizeof(Line), File) != NULL) {
        if (IsElementSection(Line, File, &NDimExceptRuleExists)) {
            if (ElementsSectionPos == -1) {
                ElementsSectionPos = CurrentPos;
            }
            long LastValidElemDataLinePos = -1;
            while (fgets(Line, sizeof(Line), File) != NULL) {
                if (LineToArray(true, false, 2, 1, Line) == 0 || LineToArray(true, true, 3, 0, Line) == 0) {
                    break;
                }
                nallelements++;
                LastValidElemDataLinePos = ftell(File);
            }
            if (fseek(File, LastValidElemDataLinePos, SEEK_SET) != 0) {
                FILE_LOG_SINGLE(ERROR, "'fseek()' call for LastValidElemDataLinePos failed");
            }
        }
        else if (NodesSectionPos == -1 && IsNodeSection(Line)) {
            bool Dim3Exists = false;
            while (!Dim3Exists && fgets(Line, sizeof(Line), File) != NULL) {
                if (ndim == 0 && strstr(Line, "nid") != NULL) {
                    // Defining ndim using x, y, z column names
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
                        NodesSectionPos = ftell(File);
                    }
                }
                else if (ndim == 3) {
                    // if ndim is 3, let's check values of rows of Z column
                    // Dim3Exists gets TRUE if there is a non-zero row
                    double *ZColumnValue;
                    if (LineToArray(false, false, 4, 1, Line, Delim, (void**)&ZColumnValue) > 0) {
                        Dim3Exists = ZColumnValue[0] != 0;
                        free(ZColumnValue);
                    }
                }
            }
            // Even each row in Z column of NODE section is 0, ndim will be 3 IF input file contains ELEMENT_BEAM section
            // For ELEMENT_SHELL and ELEMENT_SOLID files, ndim will be 2 IF there are not non-zero rows in Z column
            if (ndim == 2 || ndim == 3) {
                if (!Dim3Exists && NDimExceptRuleExists == 0) {
                    ndim = 2;
                }
            }
        }

        if (ElementsSectionPos == -1) {
            CurrentPos = ftell(File);
        }
    }

    // Checking if elements and nodes were found in mesh file
    // And checking if we can be back to elements section in the file
    int fseeked = 0;
    if (nallelements == 0 || NodesSectionPos == -1 || (fseeked = fseek(File, ElementsSectionPos, SEEK_SET)) != 0) {
        fclose(File);
        if (nallelements == 0) {
            FILE_LOG_SINGLE(ERROR, "No element found. This means input file is empty or contains invalid data");
        }
        if (NodesSectionPos == -1) {
            FILE_LOG_SINGLE(ERROR, "Node section not found. This means input file is empty or contains invalid data");
        }
        if (fseeked != 0) {
            FILE_LOG_SINGLE(ERROR, "'fseek()' call for ElementsSectionPos failed");
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
        if (IsElementSection(Line, File)) {

            long LastValidElemDataLinePos = -1;
            while (fgets(Line, sizeof(Line), File) != NULL) {
                const int np = LineToArray(true, false, 2, 1, Line);
                const int nn = LineToArray(true, true, 3, 0, Line);
                if (np == 0 || nn == 0) {
                    break;
                }
                if (i >= From && i <= To) {
                    j = j + 1;
                    eptr[j] = eptr[j - 1] + nn;
                    ConnectivitySize = ConnectivitySize + nn;
                }
                i = i + 1;
                LastValidElemDataLinePos = ftell(File);
            }

            if (fseek(File, LastValidElemDataLinePos, SEEK_SET) != 0) {
                FILE_LOG_SINGLE(ERROR, "'fseek()' call for LastValidElemDataLinePos failed");
            }
        }
    }

    // Checking if calculated size of connectivity array is valid
    // And checking if we can be back to elements section of mesh file
    if (ConnectivitySize != eptr[nelements] || fseek(File, ElementsSectionPos, SEEK_SET) != 0) {
        fclose(File);
        FreeArrays();
        if (ConnectivitySize != eptr[nelements]) {
            FILE_LOG_SINGLE(ERROR, "Size of 'eind' array is not vald");
        }
        else {
            FILE_LOG_SINGLE(ERROR, "'fseek()' call failed");
        }
        return false;
    }

    // Creating and initializing "connectivity", "pid" and "ElementType" arrays for processor
    connectivity = (int *)calloc(ConnectivitySize, sizeof(int));
    pid = (int *)calloc(nelements, sizeof(int));
    global_eid = (int *)malloc(nelements*sizeof(int));
    ElementType = (char **)calloc(nelements, sizeof(char *));
    for (int i1 = 0; i1 < nelements; i1++) {
        ElementType[i1] = (char *)calloc(MAX_ELEMENT_TYPE_SIZE, sizeof(char));
    }
    i = 0;
    int ci = 0, pi = 0;
    while (fgets(Line, sizeof(Line), File) != NULL) {
        if (IsElementSection(Line, File)) {

            long LastValidElemDataLinePos = -1;
            while (fgets(Line, sizeof(Line), File) != NULL) {
                int *PIDs, *Nodes, *eids, *sID;
                const int ne = LineToArray(true, false, 1, 1, Line, Delim, (void**)&eids);
                const int np = LineToArray(true, false, 2, 1, Line, Delim, (void**)&PIDs);
                const int nn = LineToArray(true, true, 3, 0, Line, Delim, (void**)&Nodes);
		int secID = LineToArray(true, true, 3, 0, Line, Delim, (void**)&sID);
                if (np == 0 || nn == 0 || ne == 0) {
                    if (np > 0) {
                        free(PIDs);
                    }
                    if (nn > 0) {
                        free(Nodes);
                    }
                    if (ne > 0) {
                        free(eids);
                    }
                    break;
                }
                if (i >= From && i <= To) {
                    pid[pi] = PIDs[0]-1;
                    global_eid[pi] = eids[0];

                    const int N = nn == 8 ? CalcUniqueValues(Nodes, nn) : nn;
                    if (N > 4) {
			if(!RI)
                        strcpy(ElementType[pi], "C3D8");
			else strcpy(ElementType[pi], "C3D8R");
                    }
                    else if (N == 4) {
                        strcpy(ElementType[pi], "C3D4");
                    }
                    else if (N == 2) {
                        strcpy(ElementType[pi], "T3D2");
                    }
                    else {
                        FILE_LOG_SINGLE(ERROR, "Unknown element type \"%d\" in line %s", N, Line);
                    }
                    
                    pi = pi + 1;

                    for (int j1 = 0; j1 < nn; j1++) {
                        connectivity[ci] = Nodes[j1] - 1;
                        ci = ci + 1;
                    }
                }
                free(eids);
                free(PIDs);
                free(Nodes);
                i = i + 1;
                LastValidElemDataLinePos = ftell(File);
            }

            if (fseek(File, LastValidElemDataLinePos, SEEK_SET) != 0) {
                FILE_LOG_SINGLE(ERROR, "'fseek()' call for LastValidElemDataLinePos failed");
            }
        }
    }
    assert(pi == nelements);

    // Checking if we can go to nodes section of mesh file
    if (fseek(File, NodesSectionPos, SEEK_SET) != 0) {
        fclose(File);
        FreeArrays();
        FILE_LOG_SINGLE(ERROR, "'fseek()' call failed", world_rank);
        return false;
    }

    // Initializing "coordinates" array
    int *UniqueConnectivity = (int *)malloc(ConnectivitySize * sizeof(int));
    memcpy(UniqueConnectivity, connectivity, ConnectivitySize * sizeof(int));
    qsort(UniqueConnectivity, ConnectivitySize, sizeof(int), compare);
    const int UniqueConnectivitySize = unique(UniqueConnectivity, ConnectivitySize);
    const int csize = ndim * sizeof(double);
    double *UniqueConnCoordinates = (double *)calloc(UniqueConnectivitySize, csize);
    const int nCols = ndim + 1;
    int nCopied = 0;
    while (nCopied < UniqueConnectivitySize && fgets(Line, sizeof(Line), File) != NULL) {
        double *NodeData;
        const int n = LineToArray(false, false, 1, nCols, Line, Delim, (void**)&NodeData);
        if (n == nCols) {
            const int NodeID = int(NodeData[0]) - 1;
            const int *p = (const int *)bsearch(&NodeID, UniqueConnectivity, UniqueConnectivitySize, sizeof(int), compare);
            if (p != NULL) {
                const int I = p - UniqueConnectivity;
                memcpy(&UniqueConnCoordinates[I * ndim], &NodeData[1], csize);
                nCopied = nCopied + 1;
            }
            free(NodeData);
        }
        else if (n > 0) {
            free(NodeData);
        }
        if (strncmp(Line, TITLE_END, strlen(TITLE_END)) == 0) {
            break;
        }
    }
    nNodes = 0;
    coordinates = (double *)malloc(ConnectivitySize * csize);
    for (int i1 = 0; i1 < ConnectivitySize; i1++) {
        const int *p = (const int *)bsearch(&connectivity[i1], UniqueConnectivity, UniqueConnectivitySize, sizeof(int), compare);
        if (p != NULL) {
            const int I = p - UniqueConnectivity;
            memcpy(&coordinates[i1 * ndim], &UniqueConnCoordinates[I * ndim], csize);
            nNodes = nNodes + 1;
        }
    }
    free(UniqueConnCoordinates);
    free(UniqueConnectivity);

    // Checking if "coordinates" array is OK
    if (nNodes != ConnectivitySize) {
        FreeArrays();
        FILE_LOG_SINGLE(ERROR, "Failed to initialize 'coordinates' array");
        FILE_LOG_SINGLE(ERROR, "nnodes = %d, ConnectivitySize = %d, ndim = %d", nNodes, ConnectivitySize, ndim);
    }

    // Everything is OK for now. Closing the mesh file, returning TRUE.
    fclose(File);
    return nNodes == ConnectivitySize;
}
