#include "digitalbrain.h"

static ELEMENT *ReadAllElements(const char *FileName, int *Count);
//-------------------------------------------------------------------------------------------
bool CreatePartitions(const char *FileName, int NParts, const PARTOUTPUT *Parts)
{
    // Reading elements for temporary use
    int AllElementsCount;
    ELEMENT *AllElements = ReadAllElements(FileName, &AllElementsCount);
    if (AllElements == NULL)
    {
        printf("\nERROR( proc %d ): Couldn't read elements from input file.\n", world_rank);
        return false;
    }
    
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
        return false;
    }

    NParts = nparts + 1;
    
    // Creating partitions
    PARTITION *Partitions = (PARTITION *)calloc(NParts, sizeof(PARTITION));
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
        Partitions[i].Elements = (PELEMENT *)calloc(Partitions[i].Count, sizeof(PELEMENT));
    }

    // Assigning elements to partitions
    for (int i = 0; i < NParts; i++)
    {
        for (int j = 0; j < Partitions[i].Count; j++)
        {
            if (Partitions[i].Elements[j] == NULL)
            {
                for (int k = 0; k < AllElementsCount; k++)
                {
                    if (AllElements[k].PartNumber == i)
                    {
                        Partitions[i].Elements[j] = &AllElements[k];
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
                n += Partitions[i].Elements[j]->Count;
            }
        }
    }
    
    
    // Printing partitions
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
            for (int k = 0; k < Partitions[i].Elements[j]->Count; k++)
            {
                printf("%d ", Partitions[i].Elements[j]->Nodes[k]);
            }
        }
        printf("\n");
    }
    printf("\n");

    // Freeing arrays
    for (int i = 0; i < NParts; i++)
    {
        free(Partitions[i].Eptr);
        free(Partitions[i].Elements);
    }
    free(Partitions);
    for (int i = 0; i < AllElementsCount; i++)
    {
        free(AllElements[i].Nodes);
    }
    free(AllElements);

    return true;
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
                if (LineToArray(true, true, 3, 0, Line, NULL) > 0)
                {
                    AllElementsCount++;
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
        return Result;
    }
    
    *Count = AllElementsCount;
    Result = (ELEMENT *)calloc(AllElementsCount, sizeof(ELEMENT));
    int i = 0;
    while (fgets(Line, sizeof(Line), File) != NULL)
    {
        int *Nodes;
        const int n = LineToArray(true, true, 3, 0, Line, (void**)&Nodes);
        if (n == 0)
        {
            continue;
        }
        else if (i < AllElementsCount)
        {
            Result[i].Count = n;
            Result[i].Nodes = (int *)malloc(n * sizeof(int));
            Result[i].PartNumber = 0;
            for (int j = 0; j < n; j++)
            {
                Nodes[j] -= 1; // Subtracting 1 from node numbers
            }
            memcpy(Result[i].Nodes, Nodes, n * sizeof(int));
        }
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

