#include "FemTech.h"

#include <assert.h>

#ifdef _WIN32
#include <string.h>

int StrCmpCIembed(const char *str1, const char *str2) {
    return stricmp(str1, str2);
}

#else
#include <strings.h>

int StrCmpCIembed(const char *str1, const char *str2) {
    return strcasecmp(str1, str2);
}

#endif

//-------------------------------------------------------------------------------------------
static void Free2DimArray(void **Array, const int Size);
static bool IsNodeSection(char *Line, FILE *File);
static bool IsElementSection(char *Line, FILE *File, char *Type = NULL, char *ElSetName = NULL);
//-------------------------------------------------------------------------------------------
void ReadEmbed(const char *FileName) {
    // Checking if mesh file can be opened or not
    FILE *File;
    if ((File = fopen(FileName, "rb")) == NULL) {
        FILE_LOG_SINGLE(ERROR, "Cannot open input mesh file");
        return;
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
        return;
    }
  int nelementsold = nelements;  
   for(int i = 0; i<nembedel; i++){
	for(int j = 0; j<nelementsold; j++){
	    if(embedinfo[i]==global_eid[j]-1){
		nelements = nelements + 1;
		embedproc[i] = 1;
		}
	}
   }

    // Creating and initializing "eptr" and "ElementType" array, and determining size of "connectivity" array
    const int ELSETNAME_SIZE = 80;
    char **ElementType_temp = (char **)malloc(nelements * sizeof(char *));
    char **ElSetNames_temp = (char **)malloc(nelements * sizeof(char *));
    for (int i = 0; i < nelements; i++) {
        ElementType_temp[i] = (char *)calloc(MAX_ELEMENT_TYPE_SIZE, sizeof(char));
        ElSetNames_temp[i] = (char *)calloc(ELSETNAME_SIZE, sizeof(char));
    }
    char **UniqueElSetNames = (char **)malloc(ElSetsCount * sizeof(char *));
    for (int i = 0; i < ElSetsCount; i++) {
        UniqueElSetNames[i] = (char *)calloc(ELSETNAME_SIZE, sizeof(char));
    }
    int *eptr_temp = (int *)calloc(nelements + 1, sizeof(int));

    memcpy(eptr_temp, eptr, (nelementsold + 1)*sizeof(int));

    for (int i = 0; i < nelementsold; i++) {
	memcpy(ElementType_temp[i], ElementType[i], MAX_ELEMENT_TYPE_SIZE*sizeof(char));
    }
    

    int i = 0, j = nelementsold, pi = nelementsold, ConnectivitySize = nelementsold*8, UniqueElSetsCount = 0;
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
		if ((n==2)&&(embedproc[i-nallelementsnoembed]==1)) {
                    j = j + 1;
                    eptr_temp[j] = eptr_temp[j - 1] + n;
                    ConnectivitySize = ConnectivitySize + n;
                    Type[MAX_ELEMENT_TYPE_SIZE - 1] = 0;
                    strcpy(ElementType_temp[pi], Type);
                    strcpy(ElSetNames_temp[pi], ElSetName);
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

//	printf("%d %d\n", eptr[nelements], world_rank);

    // Checking if calculated size of connectivity array is valid
    // And checking if we can be back to elements section of mesh file
    if (ConnectivitySize != eptr_temp[nelements] || fseek(File, ElementsSectionPos, SEEK_SET) != 0) {
        fclose(File);
        FreeArrays();
        Free2DimArray((void **)ElSetNames_temp, nelements);
        Free2DimArray((void **)UniqueElSetNames, ElSetsCount);
        if (ConnectivitySize != eptr[nelements]) {
            FILE_LOG_SINGLE(ERROR, "Size of 'eind' array is not valid");
        }
        else {
            FILE_LOG_SINGLE(ERROR, "'fseek()' call for ElementsSectionPos failed");
        }
        return;
    }

    // Creating and initializing "pid" array
    int *pid_temp = (int *)calloc(nelements, sizeof(int));
    memcpy(pid_temp, pid, nelementsold*sizeof(int));
    for (int i1 = 0; i1 < UniqueElSetsCount; i1++) {
        for (int j1 = nelementsold; j1 < nelements; j1++) {
            if (strcmp(UniqueElSetNames[i1], ElSetNames_temp[j1]) == 0) {
                pid_temp[j1] = i1;
            }
        }
    }
    Free2DimArray((void **)ElSetNames_temp, nelements);
    Free2DimArray((void **)UniqueElSetNames, ElSetsCount);

    int *global_eid_temp = (int *)malloc(nelements*sizeof(int));
    memcpy(global_eid_temp, global_eid, (nelementsold)*sizeof(int));

    // Creating and initializing "connectivity" array for processor
    int *connectivity_temp = (int *)calloc(ConnectivitySize, sizeof(int));
    memcpy(connectivity_temp, connectivity, (nelementsold*8)*sizeof(int));
    i = 0, j = nelementsold*8;
    int eIndex = nelementsold;
    while (fgets(Line, sizeof(Line), File) != NULL) {
        if (IsElementSection(Line, File)) {
            long LastValidElemDataLinePos = -1;
            while (fgets(Line, sizeof(Line), File) != NULL) {
                int *Nodes;
                const int n = LineToArray(true, false, 1, 0, Line, Delim, (void**)&Nodes);
                if (n == 0) {
                    break;
                }
                if ((n-1==2)&&(embedproc[i-nallelementsnoembed]==1)) {
                  global_eid_temp[eIndex] = Nodes[0];
                  eIndex = eIndex + 1;
                    for (int k = 1; k < n; k++) {
                        connectivity_temp[j] = Nodes[k] - 1;
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
        return;
    }

    // Initializing "coordinates" array
    int *UniqueConnectivity = (int *)malloc(ConnectivitySize * sizeof(int));
    memcpy(UniqueConnectivity, connectivity_temp, ConnectivitySize * sizeof(int));
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

    double *coordinates_temp = (double *)malloc(UniqueConnectivitySize * csize);
    memcpy(coordinates_temp, coordinates, nDOF*sizeof(double));
    for (int i1 = nDOF/3; i1 < UniqueConnectivitySize; i1++) {
        const int *p = (const int *)bsearch(&connectivity_temp[i1], UniqueConnectivity, UniqueConnectivitySize, sizeof(int), compare);
            memcpy(&coordinates_temp[i1 * ndim], &UniqueConnCoordinates[i1 * ndim], csize);
    }
	int tempsize = UniqueConnectivitySize-nDOF/3;
	int *renumberConnectivity = (int *)malloc(tempsize * sizeof(int));
	for(int i2 = nDOF/3, j2=1; i2<UniqueConnectivitySize; i2++,j2++){
		renumberConnectivity[j2-1] = UniqueConnectivity[nDOF/3-1]+j2;
	}
	for(int i2=0; i2<ConnectivitySize; i2++){
		for(int j2=0; j2<tempsize; j2++){
			if(connectivity_temp[i2]==UniqueConnectivity[nDOF/3+j2]){
		  	  	connectivity_temp[i2]=renumberConnectivity[j2];	
			}
		}
	}

//delete original arrays, recreate them and copy new info into arrays
  if (ElementType != NULL) {
    for (int i3 = 0; i3 < nelementsold; i3++){
        free(ElementType[i3]);
    }
  }
    free(eptr);
    free(global_eid);
    free(connectivity);
    free(coordinates);
    free(pid);


  ElementType = (char **)malloc(nelements * sizeof(char *));
  for (int i3 = 0; i3 < nelements; i3++) {
        ElementType[i3] = (char *)calloc(MAX_ELEMENT_TYPE_SIZE, sizeof(char));
  }
  eptr = (int *)calloc(nelements + 1, sizeof(int));
  global_eid = (int *)malloc(nelements*sizeof(int));
  connectivity = (int *)calloc(ConnectivitySize, sizeof(int));
  coordinates = (double *)malloc(UniqueConnectivitySize * csize);
  pid = (int *)malloc(nelements*sizeof(int));
  for (int i3 = 0; i3 < nelements; i3++) {
	memcpy(ElementType[i3], ElementType_temp[i3], MAX_ELEMENT_TYPE_SIZE*sizeof(char));
  } 
  memcpy(eptr, eptr_temp, (nelements+1)* sizeof(int));
  memcpy(global_eid, global_eid_temp, nelements* sizeof(int));
  memcpy(connectivity, connectivity_temp, ConnectivitySize* sizeof(int));
  memcpy(coordinates, coordinates_temp, UniqueConnectivitySize* csize);
  memcpy(pid, pid_temp, nelements* sizeof(int));
  nDOF = UniqueConnectivitySize*3;
  nNodes = UniqueConnectivitySize;
  free(UniqueConnCoordinates);
  free(UniqueConnectivity);
  if (ElementType != NULL) {
    for (int i3 = 0; i3 < nelements; i3++){
        free(ElementType_temp[i3]);
    }
  }	
  free(eptr_temp);
  free(global_eid_temp);
  free(connectivity_temp);
  free(coordinates_temp);
  free(pid_temp);
  fclose(File);
  return;
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
        Result = strcmp(Tokens[0], "*") == 0 && StrCmpCIembed(Tokens[1], "NODE") == 0;
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
        Result = strcmp(Tokens[0], "*") == 0 && StrCmpCIembed(Tokens[1], "ELEMENT") == 0 && strcmp(Tokens[2], ",") == 0;
    }
    if (Result) {
        Result = false;
        do {
            for (int i = 0; i < Size; i++) {
                if (!Result && StrCmpCIembed(Tokens[i], "TYPE") == 0) {
                    Result = (i + 2) < Size && strcmp(Tokens[i + 1], "=") == 0 && strlen(Tokens[i + 2]) > 0;
                    if (Result && Type != NULL) {
                        strcpy(Type, Tokens[i + 2]);
                    }
                }
                if (StrCmpCIembed(Tokens[i], "ELSET") == 0 && (i + 2) < Size && strcmp(Tokens[i + 1], "=") == 0) {
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
