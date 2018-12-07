#include "digitalbrain.h"

void FreeArrays()
{
    if (coordinates != NULL)
    {
        free(coordinates);
    }
    
    if (connectivity != NULL)
    {
        free(connectivity);
    }
    
    //free(mid);
    
    if (pid != NULL)
    {
        free(pid);
    }
    
    if (eptr != NULL)
    {
        free(eptr);
    }
    
    if (ElementType != NULL)
    {
        for (int i = 0; i < nelements; i++)
        {
            free(ElementType[i]);
        }
        free(ElementType);
    }
}
