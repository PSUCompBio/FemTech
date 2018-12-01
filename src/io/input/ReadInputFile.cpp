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
int world_rank;
int world_size;

#define MAX_LINE 128

//-------------------------------------------------------------------------------------------
/*
   Converts line into array of int or double if line contains
   valid int/double separated by whitespase or tab.
*/
static size_t LineToArray(
    const bool IntOrFloat,      // If true, will be converted to int else to double
    const bool CheckLastVal,    // If true, if next column's value is the same as the last column's one, next coulmns won't be fetched
    const int ColumnToStart,    // One-based index of column to start to fetch
    const int ColumnCount,      // Number of columns to fetch. Will be used if greater than zero
    const char *ConstLine,      // Line to convert
    void **Array);              // Output array, int or double
//-------------------------------------------------------------------------------------------
/*
   elmdist, processor's eptr and eind arrays are created in this function.
   So the caller MUST FREE them when/if the function returns true.
*/
bool ReadInputFile(
    const char *FileName,
    idx_t **elmdist,
    idx_t **MyEptr,
    idx_t **MyEind,
    size_t *partArraySize)
{
    if (world_size < 1 || world_rank < 0 || world_rank >= world_size)
    {
        return false;
    }

    FILE *File;
    if (FileName == NULL || strlen(FileName) == 0 || (File = fopen(FileName, "rt")) == NULL)
    {
        return false;
    }

    char Line[MAX_LINE];
    size_t AllElementsCount = 0;
    long ElementsSectionPos = 0;
    while (fgets(Line, sizeof(Line), File) != NULL)
    {
        if (strcmp(Line, "*ELEMENT_SOLID\n") == 0)
        {
            ElementsSectionPos = ftell(File);
            while (fgets(Line, sizeof(Line), File) != NULL)
            {
                if (LineToArray(true, true, 3, 0, Line, NULL) > 0)
                {
                    AllElementsCount++;
					//printf("P%d: Line: %s", world_rank, Line);
                }
                if (strcmp(Line, "*NODE\n") == 0)
                {
                    break;
                }
            }
        }
    }

    if (AllElementsCount == 0 || fseek(File, ElementsSectionPos, SEEK_SET) != 0)
    {
        fclose(File);
        return false;
    }

    //ELMDIST: THIS ARRAY DESCRIBES HOW THE ELEMENTS OF THE MESH ARE DISTRIBUTED AMONG THE PROCESSORS.
    //         IT IS ANALOGOUS TO THE VTXDIST ARRAY. ITS CONTENTS ARE IDENTICAL FOR EVERY PROCESSOR.
    //         Size of this array equals to p + 1, where p is count of processors       
    *elmdist = (idx_t *)calloc(world_size + 1, sizeof(idx_t));
    for (int i = 0, n = 0; i <= world_size; i++)
    {
        (*elmdist)[i] = n;
        n += (AllElementsCount / world_size);
    }
    (*elmdist)[world_size] += (AllElementsCount % world_size);

    const size_t MyElementsCount = (*elmdist)[world_rank + 1] - (*elmdist)[world_rank];
    const size_t MyEptrSize = MyElementsCount + 1;
    *partArraySize = MyElementsCount;
    size_t MyEindSize = 0;
	//printf("P%d: MyElementsCount: %d\n", world_rank, MyElementsCount);
    size_t *NodesCountPerElement = (size_t *)calloc(AllElementsCount, sizeof(size_t));
    size_t i = 0;
    while (fgets(Line, sizeof(Line), File) != NULL)
    {
        const size_t n = LineToArray(true, true, 3, 0, Line, NULL);
        if (n == 0)
        {
            continue;
        }
        else if (i < AllElementsCount)
        {
            NodesCountPerElement[i] = n;
        }
        if (strcmp(Line, "*NODE\n") == 0)
        {
            break;
        }
        i++;
    }
    const size_t From = world_rank * AllElementsCount / world_size;
    const size_t To = From + MyElementsCount - 1;
    *MyEptr = (idx_t *)calloc(MyEptrSize, sizeof(idx_t));
    for (size_t i = 0, j = 0, n = 0; i < AllElementsCount; i++)
    {
        if (i >= From && i <= To)
        {
            n += NodesCountPerElement[i];
            (*MyEptr)[++j] = n;
            MyEindSize += NodesCountPerElement[i];
        }
    }
    free(NodesCountPerElement);

    if (MyEindSize != ((*MyEptr)[MyElementsCount]) || fseek(File, ElementsSectionPos, SEEK_SET) != 0)
    {
        fclose(File);
        free(*MyEptr);
        free(*elmdist);
        return false;
    }

    i = 0;
    size_t ei = 0;
    *MyEind = (idx_t *)calloc(MyEindSize, sizeof(idx_t));
    while (fgets(Line, sizeof(Line), File) != NULL)
    {
        int *Nodes;
        const size_t n = LineToArray(true, true, 3, 0, Line, (void**)&Nodes);
        if (n == 0)
        {
            continue;
        }
        else if (i >= From && i <= To)
        {
            for (size_t j = 0; j < n; j++)
            {
                (*MyEind)[ei++] = Nodes[j];
            }
        }
        free(Nodes);
        if (strcmp(Line, "*NODE\n") == 0)
        {
            break;
        }
        i++;
    }

	idx_t *eptr = *MyEptr;
	idx_t *eind = *MyEind;
	printf("\neptr in processor %d = ", world_rank);
	for (size_t i = 0; i < MyEptrSize; i++){
		printf("%d ", eptr[i]);
	}
	printf("\n");

	printf("\neind in processor %d = \n", world_rank);
	for (size_t i = 0; i < MyEptrSize - 1; i++){
		printf("p%d, element %d: ",world_rank, i);
		for (size_t j = eptr[i]; j < eptr[i + 1]; j++){
			printf("%d ", eind[j]);
		}
		printf("\n");
	}
	printf("\n");
	
    fclose(File);
    return true;
}
//-------------------------------------------------------------------------------------------
static size_t LineToArray(
    const bool IntOrFloat, const bool CheckLastVal, const int ColumnToStart,
    const int ColumnCount, const char *ConstLine, void **Array)
{
    size_t Result = 0;
    char Line[MAX_LINE];
    strcpy(Line, ConstLine);
    double LastVal = INT_MAX;
    int Column = 0;
    char *T = strtok(Line, " \t");
    while (T != NULL)
    {
        char *EndPtr;
        const double Val = strtod(T, &EndPtr);
        if (EndPtr == T)
        {
            // Error. Token in a line cannot be converted to number.
            return 0;
        }
        T = strtok(NULL, " \t");
        if (++Column < ColumnToStart)
        {
            // Skip columns to go to required ones.
            continue;
        }
        if ((CheckLastVal && Val == LastVal) || (ColumnCount > 0 && Result >= ColumnCount))
        {
            // Stop fetching columns.
            break;
        }
        LastVal = Val;
        Result += 1;
    }

    if (Result == 0 || Array == NULL)
    {
        return Result;
    }

    strcpy(Line, ConstLine);
    *Array = calloc(Result, IntOrFloat ? sizeof(int) : sizeof(double));
    size_t i = 0;
    Column = 0;
    T = strtok(Line, " \t");
    while (T != NULL && i < Result)
    {
        const double Val = strtod(T, NULL);
        T = strtok(NULL, " \t");
        if (++Column < ColumnToStart)
        {
            continue;
        }
        if (IntOrFloat)
        {
            ((int *)*Array)[i] = Val;
        }
        else
        {
            ((double *)*Array)[i] = Val;
        }
        i++;
    }
    return Result;
}

