#include "digitalbrain.h"

static ELEMENT *ReadAllElements(const char *FileName, int *Count);
static void FreeParMETISOutput(PARTOUTPUT *Parts);
//-------------------------------------------------------------------------------------------
int CreatePartitions(const char *FileName, int NParts, PARTOUTPUT *Parts, PARTITION **PartitionsP)
{
    // Reading elements for temporary use
    int AllElementsCount;
    ELEMENT *AllElements = ReadAllElements(FileName, &AllElementsCount);
    if (AllElements == NULL)
    {
        printf("\nERROR( proc %d ): Couldn't read elements from input file.\n", world_rank);
        FreeParMETISOutput(Parts);
        return 0;
    }
    
    // Printing all elements
    printf("\nAll elements = ");
    for (int i = 0; i < AllElementsCount; i++)
    {
        printf(" (%d)  ", i);
        for (int j = 0; j < AllElements[i].Count; j++)
        {
            printf("%d ", AllElements[i].Nodes[j]);
        }
    }
    printf("\n");

    // Calculating number of parts using part arrays and assigning part number/id for each element
    int nparts = 0;
    for (int i = 0, ei = 0; i < world_size; i++)
    {
        for (int j = 0; j < Parts[i].part_size; j++)
        {
            if (Parts[i].part[j] > nparts)
            {
                nparts = Parts[i].part[j];
            }
            AllElements[ei++].PartNumber = Parts[i].part[j];
        }
    }
    
    // Checking if part arrays received from processors are valid, and return if not
    if (nparts >= NParts)
    {
        printf("\nERROR( proc %d ): One or more 'part' arrays received from processors may be not valid: %d / %d\n", world_rank, NParts, nparts);
        for (int i = 0; i < AllElementsCount; i++)
        {
            free(AllElements[i].Nodes);
        }
        free(AllElements);
        FreeParMETISOutput(Parts);
        return 0;
    }

    NParts = nparts + 1;
    
    // Creating partitions
    PARTITION *Partitions = (PARTITION *)calloc(NParts, sizeof(PARTITION));
    *PartitionsP = Partitions;
    for (int i = 0; i < NParts; i++)
    {
        for (int j = 0; j < world_size; j++)
        {
            for (int k = 0; k < Parts[j].part_size; k++)
            {
                if (Parts[j].part[k] == i)
                {
                    Partitions[i].Count++;
                }
            }
        }
        Partitions[i].Elements = (ELEMENT *)calloc(Partitions[i].Count, sizeof(ELEMENT));
    }

    // Assigning elements to partitions
    for (int i = 0; i < NParts; i++)
    {
        for (int j = 0; j < Partitions[i].Count; j++)
        {
            if (Partitions[i].Elements[j].Nodes == NULL)
            {
                for (int k = 0; k < AllElementsCount; k++)
                {
                    if (AllElements[k].PartNumber == i)
                    {
                        Partitions[i].Elements[j] = AllElements[k];
                        Partitions[i].Elements[j].Nodes = (int *)calloc(AllElements[k].Count, sizeof(int));
                        memcpy(Partitions[i].Elements[j].Nodes, AllElements[k].Nodes, AllElements[k].Count * sizeof(int));
                        free(AllElements[k].Nodes);
                        AllElements[k].Nodes = NULL;
                        AllElements[k].PartNumber = -1;
                        break;
                    }
                }
            }
        }
    }
    
    // Creating eptr array for each partitiion
    for (int i = 0; i < NParts; i++)
    {
        Partitions[i].Eptr = (idx_t *)calloc(Partitions[i].Count + 1, sizeof(idx_t));
        for (int j = 0, n = 0; j <= Partitions[i].Count; j++)
        {
            Partitions[i].Eptr[j] = n;
            if (j < Partitions[i].Count)
            {
                n += Partitions[i].Elements[j].Count;
            }
        }
    }
    
    
    // Printing partitions
    for (int i = 0; i < NParts; i++)
    {
        printf("\nPartition %d:\n", i);
        printf(" eptr = ");
        for (int j = 0; j <= Partitions[i].Count; j++)
        {
            printf("%d ", Partitions[i].Eptr[j]);
        }
        printf("\n eind =");
        for (int j = 0; j < Partitions[i].Count; j++)
        {
            printf(" (%d)  ", j);
            for (int k = 0; k < Partitions[i].Elements[j].Count; k++)
            {
                printf("%d ", Partitions[i].Elements[j].Nodes[k]);
            }
        }
        printf("\n");
    }
    printf("\n");

    // Freeing arrays
    FreeParMETISOutput(Parts);
    for (int i = 0; i < AllElementsCount; i++)
    {
        if (AllElements[i].Nodes != NULL)
        {
            free(AllElements[i].Nodes);
        }
    }
    free(AllElements);
    
    return NParts;
}

//-------------------------------------------------------------------------------------------
static ELEMENT *ReadAllElements(const char *FileName, int *Count)
{
    ELEMENT *Result = NULL;
    
    FILE *File = fopen(FileName, "rt");
    if (File == NULL)
    {
        return Result;
    }
    
    char Line[MAX_FILE_LINE];
    int AllElementsCount = 0;
    long ElementsSectionPos = 0;
    while (fgets(Line, sizeof(Line), File) != NULL)
    {
        if (strcmp(Line, "*ELEMENT_SOLID\n") == 0)
        {
            ElementsSectionPos = ftell(File);
            while (fgets(Line, sizeof(Line), File) != NULL)
            {
                if (LineToArray(true, true, 1, 0, Line, NULL) > 0)
                {
                    AllElementsCount++;
                }
                if (strcmp(Line, "*NODE\n") == 0)
                {
                    break;
                }
            }
            break;
        }
    }

    if (AllElementsCount == 0 || fseek(File, ElementsSectionPos, SEEK_SET) != 0)
    {
        fclose(File);
        return Result;
    }
    
    *Count = AllElementsCount;
    Result = (ELEMENT *)calloc(AllElementsCount, sizeof(ELEMENT));
    int i = 0;
    while (fgets(Line, sizeof(Line), File) != NULL)
    {
        int *IDs, *Nodes;
        const int ni = LineToArray(true, false, 1, 2, Line, (void**)&IDs);
        const int nn = LineToArray(true, true, 3, 0, Line, (void**)&Nodes);
        if (ni != 2 || nn == 0)
        {
            continue;
        }
        else if (i < AllElementsCount)
        {
            Result[i].Id = IDs[0];
            Result[i].PId = IDs[1];
            Result[i].Count = nn;
            Result[i].Nodes = (int *)malloc(nn * sizeof(int));
            Result[i].PartNumber = 0;
            memcpy(Result[i].Nodes, Nodes, nn * sizeof(int));
        }
        free(IDs);
        free(Nodes);
        if (strcmp(Line, "*NODE\n") == 0)
        {
            break;
        }
        i++;
    }

    fclose(File);
    return Result;
}
//-------------------------------------------------------------------------------------------
static void FreeParMETISOutput(PARTOUTPUT *Parts)
{
    for (int i = 0; i < world_size; i++)
    {
        free(Parts[i].part);
    }
    free(Parts);
}
