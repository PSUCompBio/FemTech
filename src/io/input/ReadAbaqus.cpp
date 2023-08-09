#include "FemTech.h"

#include <assert.h>

#ifdef _WIN32
#include <string.h>

int StrCmpCI(const char *str1, const char *str2) {
    return stricmp(str1, str2);
}

#else
#include <strings.h>

int StrCmpCI(const char *str1, const char *str2) {
    return strcasecmp(str1, str2);
}

#endif

//-------------------------------------------------------------------------------------------
static void Free2DimArray(void **Array, const int Size);
static bool IsNodeSection(char *Line, FILE *File);
static bool IsElementSection(char *Line, FILE *File, char *Type = NULL, char *ElSetName = NULL);
//-------------------------------------------------------------------------------------------
bool ReadAbaqus(const char *FileName) {
    // Checking if mesh file can be opened or not
    FILE *File;
    if ((File = fopen(FileName, "rb")) == NULL) {
        FILE_LOG_SINGLE(ERROR, "Cannot open input mesh file");
        return false;
    }
    
    // Counting number of elements in elements section, number of nodes in nodes section of mesh file
    // And fixing elements and nodes sections' positions in the file
    const char *Delim = ", \t";
    char Line[MAX_FILE_LINE];
    int ElSetsCount = 0;
    long ElementsSectionPos = -1, NodesSectionPos = -1, CurrentPos = ftell(File);
    while (fgets(Line, sizeof(Line), File) != NULL) {
        if (NodesSectionPos == -1 && IsNodeSection(Line, File)) {
            const long P = ftell(File);
            if (fgets(Line, sizeof(Line), File) != NULL) {
                const int n = LineToArray(false, false, 1, 0, Line, Delim) - 1;
                if (n >= 2) {
                    ndim = n > 3 ? 3 : n;
                    NodesSectionPos = P;
                }
            }
        }
        
        if (IsElementSection(Line, File)) {
            if (ElementsSectionPos == -1) {
                ElementsSectionPos = CurrentPos;
            }
            long LastValidElemDataLinePos = -1;
            ElSetsCount++;
            while (fgets(Line, sizeof(Line), File) != NULL) {
                if (LineToArray(true, false, 2, 0, Line, Delim) == 0) {
                    break;
                }
                nallelements++;
                LastValidElemDataLinePos = ftell(File);
            }
            if (fseek(File, LastValidElemDataLinePos, SEEK_SET) != 0) {
                FILE_LOG_SINGLE(ERROR, "'fseek()' call for LastValidElemDataLinePos failed");
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
    nallelementsnoembed = nallelements - nembedel; //only partitioning on host elements
    nelements = nallelementsnoembed / world_size;
    if (world_rank == world_size - 1) {
        nelements += (nallelementsnoembed % world_size);
    }

    // Fixing range of lines (in elements section of mesh file) from which processor will read its own elements
    const int From = world_rank * (nallelementsnoembed / world_size);
    const int To = From + nelements - 1;
    
 /*  for(int i = 0; i<nembedel; i++){
	if(embedinfo[i]>=From && embedinfo[i]<=To){
	nelements = nelements + 1;
	embedproc[i] = 1;
	}
   }*/

    // Creating and initializing "eptr" and "ElementType" array, and determining size of "connectivity" array
    const int ELSETNAME_SIZE = 80;
    ElementType = (char **)malloc(nelements * sizeof(char *));
    char **ElSetNames = (char **)malloc(nelements * sizeof(char *));
    for (int i = 0; i < nelements; i++) {
        ElementType[i] = (char *)calloc(MAX_ELEMENT_TYPE_SIZE, sizeof(char));
        ElSetNames[i] = (char *)calloc(ELSETNAME_SIZE, sizeof(char));
    }
    char **UniqueElSetNames = (char **)malloc(ElSetsCount * sizeof(char *));
    for (int i = 0; i < ElSetsCount; i++) {
        UniqueElSetNames[i] = (char *)calloc(ELSETNAME_SIZE, sizeof(char));
    }
    eptr = (int *)calloc(nelements + 1, sizeof(int));
    int i = 0, j = 0, pi = 0, ConnectivitySize = 0, UniqueElSetsCount = 0;
    char Type[MAX_FILE_LINE] = {0}, ElSetName[MAX_FILE_LINE] = {0};
    while (fgets(Line, sizeof(Line), File) != NULL) {
        if (IsElementSection(Line, File, Type, ElSetName)) {

            ElSetName[ELSETNAME_SIZE - 1] = 0;
            bool Found = false;
            for (int k = 0; !Found && k < ElSetsCount; k++) {
                Found = strcmp(ElSetName, UniqueElSetNames[k]) == 0;
            }
            if (!Found && UniqueElSetsCount < ElSetsCount) {
                strcpy(UniqueElSetNames[UniqueElSetsCount], ElSetName);
                UniqueElSetsCount = UniqueElSetsCount + 1;
            }

            long LastValidElemDataLinePos = -1;
            while (fgets(Line, sizeof(Line), File) != NULL) {
                const int n = LineToArray(true, false, 2, 0, Line, Delim);
                if (n == 0) {
                    break;
                }
           //     if ((i >= From && i <= To)||embedproc[i-nallelementsnoembed]==1) {
		if (i >= From && i <= To) {
                    j = j + 1;
                    eptr[j] = eptr[j - 1] + n;
                    ConnectivitySize = ConnectivitySize + n;
                    Type[MAX_ELEMENT_TYPE_SIZE - 1] = 0;
                    if (!strcmp(Type, "C3D8") && reducedIntegration) {
                      strncpy(Type, "C3D8R", MAX_ELEMENT_TYPE_SIZE);
                    }
                    strcpy(ElementType[pi], Type);
                    strcpy(ElSetNames[pi], ElSetName);
                    pi = pi + 1;
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
        Free2DimArray((void **)ElSetNames, nelements);
        Free2DimArray((void **)UniqueElSetNames, ElSetsCount);
        if (ConnectivitySize != eptr[nelements]) {
            FILE_LOG_SINGLE(ERROR, "Size of 'eind' array is not valid");
        }
        else {
            FILE_LOG_SINGLE(ERROR, "'fseek()' call for ElementsSectionPos failed");
        }
        return false;
    }

    // Creating and initializing "pid" array
    pid = (int *)calloc(nelements, sizeof(int));
    for (int i1 = 0; i1 < UniqueElSetsCount; i1++) {
        for (int j1 = 0; j1 < nelements; j1++) {
            if (strcmp(UniqueElSetNames[i1], ElSetNames[j1]) == 0) {
                pid[j1] = i1;
            }
        }
    }
    Free2DimArray((void **)ElSetNames, nelements);
    Free2DimArray((void **)UniqueElSetNames, ElSetsCount);

    global_eid = (int *)malloc(nelements*sizeof(int));

    // Creating and initializing "connectivity" array for processor
    connectivity = (int *)calloc(ConnectivitySize, sizeof(int));
    i = 0, j = 0;
    int eIndex = 0;
    while (fgets(Line, sizeof(Line), File) != NULL) {
        if (IsElementSection(Line, File)) {
            long LastValidElemDataLinePos = -1;
            while (fgets(Line, sizeof(Line), File) != NULL) {
                int *Nodes;
                const int n = LineToArray(true, false, 1, 0, Line, Delim, (void**)&Nodes);
                if (n == 0) {
                    break;
                }
          //      if ((i >= From && i <= To||embedproc[i-nallelementsnoembed]==1)) {
                if (i >= From && i <= To) {
                  global_eid[eIndex] = Nodes[0];
                  eIndex = eIndex + 1;
                    for (int k = 1; k < n; k++) {
                        connectivity[j] = Nodes[k] - 1;
                        j = j + 1;
                    }
                }
                free(Nodes);
                i = i + 1;
                LastValidElemDataLinePos = ftell(File);
            }
            
            if (fseek(File, LastValidElemDataLinePos, SEEK_SET) != 0) {
                FILE_LOG_SINGLE(ERROR, "'fseek()' call for LastValidElemDataLinePos failed");
            }
        }
    }
    assert(eIndex == nelements);
    
    // Checking if we can go to nodes section of mesh file
    if (fseek(File, NodesSectionPos, SEEK_SET) != 0) {
        fclose(File);
        FreeArrays();
        FILE_LOG_SINGLE(ERROR, "'fseek()' call for NodesSectionPos failed");
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
        else {
            if (n > 0) {
                free(NodeData);
            }
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

    fclose(File);
    return nNodes == ConnectivitySize;
}
//-------------------------------------------------------------------------------------------
static char **LineToTokens(const char *Line, int *Size) {
    char **Result = NULL;
    char StrToken[MAX_FILE_LINE];
    int Count;
    
    for (int n = 0; n < 2; n++) {
        StrToken[0] = 0;
        Count = 0;
        for (int i = 0; i < strlen(Line); i++) {
            const char Ch = Line[i];
            if (Ch == ' ' || Ch == '\t' || Ch == '\r' || Ch == '\n') {
                if (strlen(StrToken) != 0) {
                    if (Result != NULL) {
                        Result[Count] = (char *)malloc((strlen(StrToken) + 1) * sizeof(char));
                        strcpy(Result[Count], StrToken);
                    }
                    StrToken[0] = 0;
                    Count++;
                }
            }
            else if ((Ch >= 'a' && Ch <= 'z') || (Ch >= 'A' && Ch <= 'Z') || (Ch >= '0' && Ch <= '9') || Ch == '_' || Ch == '-') {
                strncat(StrToken, &Line[i], 1);
                if (i == (strlen(Line) - 1)) {
                    if (Result != NULL) {
                        Result[Count] = (char *)malloc((strlen(StrToken) + 1) * sizeof(char));
                        strcpy(Result[Count], StrToken);
                    }
                    Count++;
                }
            }
            else {
                if (strlen(StrToken) != 0) {
                    if (Result != NULL) {
                        Result[Count] = (char *)malloc((strlen(StrToken) + 1) * sizeof(char));
                        strcpy(Result[Count], StrToken);
                    }
                    StrToken[0] = 0;
                    Count++;
                }
                if (Result != NULL) {
                    Result[Count] = (char *)calloc(2, sizeof(char));
                    strncat(Result[Count], &Line[i], 1);
                }
                Count++;
            }
        }
        *Size = Count;
        if (Count == 0) {
            return Result;
        }
        else if (Result != NULL) {
            continue;
        }
        Result = (char **)malloc(Count * sizeof(char *));
    }
    return Result;
}
//-------------------------------------------------------------------------------------------
static void Free2DimArray(void **Array, const int Size) {
    if (Array == NULL) {
        return;
    }
    for (int i = 0; i < Size; i++) {
        free(Array[i]);
    }
    free(Array);
}
//-------------------------------------------------------------------------------------------
static bool IsLineComplete(const char *Line) {
    for (int i = strlen(Line) - 1; i >= 0; i--) {
        if (Line[i] == '\r' || Line[i] == '\n' || Line[i] == ' ' || Line[i] == '\t') {
            continue;
        }
        else {
            return Line[i] != ',';
        }
    }
    return true;
}
//-------------------------------------------------------------------------------------------
static bool IsNodeSection(char *Line, FILE *File) {
    bool Result = false;
    int Size;
    char **Tokens = LineToTokens(Line, &Size);
    if (Size >= 2) {
        Result = strcmp(Tokens[0], "*") == 0 && StrCmpCI(Tokens[1], "NODE") == 0;
    }
    if (Result && Size > 2) {
        Result = strcmp(Tokens[2], ",") == 0;
        if (Result) {
            while (!IsLineComplete(Line) && fgets(Line, MAX_FILE_LINE * sizeof(char), File) != NULL) {
                // Going to beginning of data line
            }
        }
    }
    Free2DimArray((void **)Tokens, Size);
    return Result;
}
//-------------------------------------------------------------------------------------------
static bool IsElementSection(char *Line, FILE *File, char *Type, char *ElSetName) {
    bool Result = false;
    int Size;
    char **Tokens = LineToTokens(Line, &Size);
    if (Size > 2) {
        Result = strcmp(Tokens[0], "*") == 0 && StrCmpCI(Tokens[1], "ELEMENT") == 0 && strcmp(Tokens[2], ",") == 0;
    }
    if (Result) {
        Result = false;
        do {
            for (int i = 0; i < Size; i++) {
                if (!Result && StrCmpCI(Tokens[i], "TYPE") == 0) {
                    Result = (i + 2) < Size && strcmp(Tokens[i + 1], "=") == 0 && strlen(Tokens[i + 2]) > 0;
                    if (Result && Type != NULL) {
                        strcpy(Type, Tokens[i + 2]);
                    }
                }
                if (StrCmpCI(Tokens[i], "ELSET") == 0 && (i + 2) < Size && strcmp(Tokens[i + 1], "=") == 0) {
                    if (ElSetName != NULL) {
                        strcpy(ElSetName, Tokens[i + 2]);
                    }
                }
            }
            if (IsLineComplete(Line)) {
                break;
            }
            if (fgets(Line, MAX_FILE_LINE * sizeof(char), File) != NULL && !Result) {
                Free2DimArray((void **)Tokens, Size);
                Tokens = LineToTokens(Line, &Size);
            }
        } while (true);
    }
    Free2DimArray((void **)Tokens, Size);
    return Result;
}
