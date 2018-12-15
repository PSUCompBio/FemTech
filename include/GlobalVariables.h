extern int nparts;
extern int nelements;
extern int nnodes;
extern int nCoordinates;
extern int ndim;
extern int world_rank;
extern int world_size;


extern double *coordinates;
extern int *connectivity;
extern int *eptr; /*number of nodes per element pointer*/
extern int *pid; /* part ID */
extern int *mid; /*material ID */
extern char **ElementType; /* element type, e.g. C3D8 */
