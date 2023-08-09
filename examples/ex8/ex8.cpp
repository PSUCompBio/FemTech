#include "FemTech.h"

#include <string>

bool ImplicitDynamic = true;
bool ExplicitDynamic = false;

/*Drupal*/
bool reducedIntegration = true;
double tInitial = 0.000, Time = 0.0;
bool embed = false; /*Drupal*/

int main(int argc, char **argv){

	// Initialize the MPI environment
	MPI_Init(NULL, NULL);
	// Get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (argc < 1) {
    fprintf(stdout, "ERROR(%d): Please provide an input mesh file\n", world_rank);
    fprintf(stdout, "ERROR(%d): Usage %s meshfile\n", world_rank, argv[0]);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  std::string logFile = "femtech.log";

  // Initialize the output log file
  initLog(logFile.c_str());
	
	ReadInputFile(argv[1]);
  ReadMaterials();
	PartitionMesh();
	//test

#if 0
	// Printing local arrays of processor (this section can be removed)
	printf("\neptr array in processor %d after partitioning = ", world_rank);
	for (int i = 0; i <= nelements; i++) {
		printf("%d ", eptr[i]);
	}
	printf("\n");
	printf("\neind array in processor %d after partitioning =", world_rank);
	for (int i = 0; i < nelements; i++) {
		printf(" (%d)  ", i);
		for (int j = eptr[i]; j < eptr[i + 1]; j++) {
			printf("%d ", connectivity[j]);
		}
	}
	printf("\n");
	printf("\nType/PartID of element in processor %d after partitioning = ", world_rank);
	for (int i = 0; i < nelements; i++) {
		printf("%s/%d  ", ElementType[i], pid[i]);
	}
	printf("\n");
	printf("\nSize of coordinates array in processor %d after partitioning = %d\n", world_rank, nDOF);
	printf("\nCoordinates array in processor %d after partitioning =", world_rank);
	for (int i = 0; i < nNodes; i++) {
		printf(" (%d)  ", i);
		for (int j = 0; j < ndim; j++) {
			printf("%.*f ", 1, coordinates[ndim * i + j]);
		}
	}
	printf("\n");
#endif

  AllocateArrays();

  std::string meshFile(argv[1]);
  size_t lastindex = meshFile.find_last_of(".");
  std::string outputFileName = meshFile.substr(0, lastindex);

  WriteVTU(outputFileName.c_str(), 0);
  ShapeFunctions();

  finaliseLog();
  FreeArrays();
  MPI_Finalize();

  return 0;
}
