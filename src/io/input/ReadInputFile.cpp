#include "digitalbrain.h"

/* Global Variables */
int nparts=0;
int nelements=0;
int nnodes=0;
int nCoordinates=0;
int ndim=0;
double *coordinates;
int *connectivity;
int *eptr;
int *pid;
int *mid;
char **ElementType;
int world_rank;
int world_size;


//-------------------------------------------------------------------------------------------
/*
   elmdist, processor's eptr and eind arrays are created in this function.
*/
bool ReadInputFile(const char *FileName){
    if (world_size < 1 || world_rank < 0 || world_rank >= world_size) {
        printf("\nERROR( proc %d ): 'world_size' and/or 'world_rank' variable is not valid.\n", world_rank);
        return false;
    }

    FILE *File;
    if (FileName == NULL || strlen(FileName) == 0 || (File = fopen(FileName, "rt")) == NULL) {
        printf("\nERROR( proc %d ): Failed to open input file.\n", world_rank);
        return false;
    }

    char Line[MAX_FILE_LINE];
    int AllElementsCount = 0;
    long ElementsSectionPos = 0;
    while (fgets(Line, sizeof(Line), File) != NULL) {
        if (strcmp(Line, "*ELEMENT_SOLID\n") == 0) {
            ElementsSectionPos = ftell(File);
            while (fgets(Line, sizeof(Line), File) != NULL) {
                if (LineToArray(true, true, 3, 0, Line, NULL) > 0) {
                    AllElementsCount++;
                }
                if (strcmp(Line, "*NODE\n") == 0) {
                    break;
                }
            }
            break;
        }
    }

    if (AllElementsCount == 0 || fseek(File, ElementsSectionPos, SEEK_SET) != 0) {
        fclose(File);
        if (AllElementsCount == 0) {
            printf("\nERROR( proc %d ): No elemens found. This means input file is empty or contains invalid data.\n", world_rank);
        }
        else {
            printf("\nERROR( proc %d ): 'fseek()' call failed.\n", world_rank);
        }
        return false;
    }
    
    nelements = AllElementsCount / world_size;
    if (world_rank == world_size - 1) {
        nelements += (AllElementsCount % world_size);
    }

    int *NodesCountPerElement = (int *)calloc(AllElementsCount, sizeof(int));
    int i = 0;
    while (fgets(Line, sizeof(Line), File) != NULL) {
        const int n = LineToArray(true, true, 3, 0, Line, NULL);
        if (n == 0) {
            continue;
        }
        else if (i < AllElementsCount) {
            NodesCountPerElement[i] = n;
        }
        if (strcmp(Line, "*NODE\n") == 0) {
            break;
        }
        i++;
    }
    int ConnectivitySize = 0;
    const int From = world_rank * AllElementsCount / world_size;
    const int To = From + nelements - 1;
    eptr = (int *)calloc(nelements + 1, sizeof(int));
    for (int i = 0, j = 0, n = 0; i < AllElementsCount; i++) {
        if (i >= From && i <= To) {
            n += NodesCountPerElement[i];
            eptr[++j] = n;
            ConnectivitySize += NodesCountPerElement[i];
        }
    }
    free(NodesCountPerElement);

    if (ConnectivitySize != eptr[nelements] || fseek(File, ElementsSectionPos, SEEK_SET) != 0) {
        fclose(File);
        free(eptr);
        eptr = NULL;
        if (ConnectivitySize != eptr[nelements]) {
            printf("\nERROR( proc %d ): Size of 'eind' array is not vald.\n", world_rank);
        }
        else {
            printf("\nERROR( proc %d ): 'fseek()' call failed.\n", world_rank);
        }
        return false;
    }

    i = 0;
    int ei = 0;
    connectivity = (int *)calloc(ConnectivitySize, sizeof(int));
    while (fgets(Line, sizeof(Line), File) != NULL) {
        int *Nodes;
        const int n = LineToArray(true, true, 3, 0, Line, (void**)&Nodes);
        if (n == 0) {
            continue;
        }
        else if (i >= From && i <= To) {
            for (int j = 0; j < n; j++) {
                connectivity[ei++] = Nodes[j] - 1;
            }
        }
        free(Nodes);
        if (strcmp(Line, "*NODE\n") == 0) {
            break;
        }
        i++;
    }

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
        if (++Column < ColumnToStart) {
            // Skip columns to go to required ones.
            continue;
        }
        if ((CheckLastVal && Val == LastVal) || (ColumnCount > 0 && Result >= ColumnCount)) {
            // Stop fetching columns.
            break;
        }
        LastVal = Val;
        Result += 1;
    }

    if (Result == 0 || Array == NULL) {
        return Result;
    }

    strcpy(Line, ConstLine);
    *Array = calloc(Result, IntOrFloat ? sizeof(int) : sizeof(double));
    int i = 0;
    Column = 0;
    T = strtok(Line, " \t");
    while (T != NULL && i < Result) {
        const double Val = strtod(T, NULL);
        T = strtok(NULL, " \t");
        if (++Column < ColumnToStart) {
            continue;
        }
        if (IntOrFloat) {
            ((int *)*Array)[i] = Val;
        }
        else {
            ((double *)*Array)[i] = Val;
        }
        i++;
    }
    return Result;
}

