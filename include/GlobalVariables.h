
static double HUGE = 1e20;

extern int nparts;
extern int nelements;
extern int nallelements;
extern int nnodes;
extern int ndim;
extern int world_rank;
extern int world_size;

extern double ExplicitTimeStepReduction;
extern double FailureTimeStep;

extern double *coordinates;
extern int *connectivity;
extern int *GaussPoints;	/*holds how many guass points per element*/
extern int *gptr;			/*gauss point pointer - helps step through shp array*/
extern int *dsptr;			/*deriviative of shp functions pointer array - helps step through dshp array*/
extern int *gpPtr;			/*gauss point pointer - helps step through detJ and gaussWeights*/
extern int *fptr; /*deformation gradient pointer - helps step through F array */
extern int *cptr; /*counter array for iterating through cauchy stress */
extern int *nShapeFunctions;/*number of shp functions per element */
extern double *shp;			/*shape functions*/
extern double *dshp;		/*derivatives of shape functions*/
extern int *eptr;			/*number of nodes per element pointer*/
extern int *detFptr; //pointer for iterating through detF array.
extern int *pid;			/* part ID */
extern int *mid;			/*material ID */
extern double *mass;        /*mass matrix*/
extern double *stiffness;        /*stiffness matrix*/
extern double *rhs;              /*to store rhs matrix equation (implicit) */
extern char **ElementType;	/* element type, e.g. C3D8 */

extern double *C; /*Stores C matrix for isotropic elastic material */
extern double rho; /*Stores material density for isotropic material */
extern double *gaussWeights;
extern double *detJacobian;

/* For BC */
extern int nSpecifiedDispBC;

extern double *displacements;
extern int *boundary;
/* Variables required for unsteady problem */
extern double *velocities;
extern double * velocities_half;
extern double *accelerations;
extern bool ImplicitStatic;
extern bool ImplicitDynamic;
extern bool ExplicitDynamic;

extern double *fe;
extern double *fi;
extern double *fi_prev;
extern double *f_net;
extern double *fr_prev;
extern double *fr_curr;
extern double *f_damp_prev;
extern double *f_damp_curr;
extern double *displacements_prev;

extern double *F; /* deformation graident tensor array */
extern double *detF; /*determinate of F for all gauss points */
extern double *invF; /*Inverse of F for all gauss points */
extern double *b; /*Left Cauchy  Greeen tensor */
extern double *E; /*Green Lagrange Strain Tensor for each Gauss point */
extern double *cauchy; /*Cauchy Stress */

extern int *materialID; /* material id for each element */
extern double *properties; /* holds material parameters for each element */
static int MAXMATPARAMS = 10; /* maximum number of material parameters stored for each element */



static int debug = 1;

// Variables to keep track of time ans step count
// Initial configuration nStep = 0 and time = 0
// For steady, solution is stored at nStep = 1 and time = 1
// For unsteady nStep and time are set by solver
extern int nStep;
extern double Time;
