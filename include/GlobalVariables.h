#include <map>
#include <string>

/* Currently Ogden model uses 9 properties + 13 properties for
 * viscoelastic properties*/
static const int MAXMATPARAMS = 22; /* maximum number of material parameters stored for each element */
static const int MAXINTERNALVARS = 5; /* max internal variables per gauss point */
static const int MAXPLOTSTEPS = 1000; /* maximum plot output steps */
static const int MAXOGDEN = 3; /* maximum terms for Ogden material model */
/* Weights for parmetis generated by running Benchmark in the order of material
   ID */
static const int materialWeight[] = {1, 10, 10, 10, 15, 20, 10, 20, 20};
/* Weights based on element type 0 = C3D8, 1 = C3D4, 2 = T3D2 */
static const std::map<std::string, int> elemID = {{"C3D8", 0}, {"C3D4", 1}, {"T3D2", 2}, {"C3D8R", 3}};
static const int elemWeight[] = {2, 1, 1, 1};

const double huge = 1e20;

extern std::string uid; /* Simulation unique id */
extern int nparts;
extern int nelements;
extern int nallelements;
extern int nNodes;
extern int ndim;
extern int nDOF;
extern int nPID, nPIDglobal; /* Local and global number of unique parts*/
extern int world_rank;
extern int world_size;

extern int ndim2; // To store ndim*ndim, a widely used constant

extern double ExplicitTimeStepReduction;
extern double FailureTimeStep;

extern double *coordinates;
extern int *connectivity;
extern int *globalNodeID; /*Stores global node ID for computing axis of rotation */
extern int *GaussPoints;	/*holds how many guass points per element*/
extern int *gptr;			/*gauss point pointer - helps step through shp array*/
extern int *dsptr;			/*deriviative of shp functions pointer array - helps step through dshp array*/
extern int *gpPtr;			/*gauss point pointer - helps step through detJ and gaussWeights*/
extern int *fptr; /*deformation gradient pointer - helps step through F array */
extern int *pk2ptr; /*counter array for iterating through PK2 stress */
extern int *nShapeFunctions;/*number of shp functions per element */
extern double *shp;			/*shape functions*/
extern double *dshp;		/*derivatives of shape functions*/
extern int *eptr;			/*number of nodes per element pointer*/
extern int *detFptr; //pointer for iterating through detF array.
extern int *pid;			/* part ID */
extern int *global_eid; /* Variable to store global element id for post processing */
extern double *mass;        /*mass matrix*/
extern double *stiffness;        /*stiffness matrix*/
extern double *rhs;              /*to store rhs matrix equation (implicit) */
extern char **ElementType;	/* element type, e.g. C3D8 */

extern double *C; /*Stores C matrix for isotropic elastic material */
extern double *gaussWeights;
extern double *detJacobian;

extern int countstress;

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

extern double *fe; /*external forces for current step - used in engery calculation*/
extern double *fe_prev; /*external forces for previous step - used in engery calculation*/
extern double *fi; /*internal forces for current step - used in engery calculation*/
extern double *fi_prev; /*internal forces for previous step - used in engery calculation*/
extern double *f_net;
extern double *fr_prev;
extern double *fr_curr;
extern double *f_damp_prev;
extern double *f_damp_curr;
extern double *displacements_prev; /*displacements for previous step - used in engery calculation*/
extern double *accelerations_prev; /*displacements for previous step - used in engery calculation*/

extern double *F; /* deformation graident tensor array */
extern double *detF; /*determinate of F for all gauss points */
extern double *invF; /*Inverse of F for all gauss points */
extern double *pk2; /*PK2 Stress */
extern double *Eavg; /*Green Lagrange Strain Tensor Averaged */

extern int *materialID; /* material id for each element */
extern double *properties; /* holds material parameters for each element */

extern double *internals; /* internal variables, typically used for damage or plasticity varialbes*/
extern int *InternalsPtr; /* pointer for iterating through internals array */

static int debug = 1; /* Setting it to 1 has no effect when compiled in non-debug mode */

// Variables to keep track of time ans step count
// Initial configuration nStep = 0 and time = 0
// For steady, solution is stored at nStep = 1 and time = 1
// For unsteady nStep and time are set by solver
extern double Time;
extern double dt;
extern double *stepTime;

// Variables to keep store the communication patterns between processes
// Implemented in PartitionMesh.cpp
extern int sendProcessCount;
extern int *sendProcessID;
extern int *sendNeighbourCount;
extern int *sendNeighbourCountCum;
extern int *sendNodeIndex;
extern double *sendNodeDisplacement;
extern double *recvNodeDisplacement;

// File to write energy
extern FILE *energyFile;

// Variables to be allocated based on material models used
extern double* mat1;
extern double* mat2;
extern double* mat3;
extern double* mat4;
// Internal Force Update
extern double *fintGQ;
extern double *B;
// Variables to store stress history in viscoelastic materials
extern double **Hn;
extern double **S0n;
extern int *nProny;

// Variables required for updated Lagrangian formulation
extern double *F_Xi, *F_XiInverse, J_Xi;
extern double *sigma_n; // To store stress in Voigot notation for updated Lagrangian
extern double *B0; // Store dN_I/dX_j. Used for H computation
// Temperory variables to test updated lagrangian
// extern double *F_Xi_0;
