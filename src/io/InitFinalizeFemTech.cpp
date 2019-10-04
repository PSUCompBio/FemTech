#include "FemTech.h"

#include "mpi.h"

// int world_rank;
// int world_size;

void InitFemTech() {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  const char *logFile = "femtech.log";
  initLog(logFile);
}

void FinalizeFemTech() {
  finaliseLog();
  FreeArrays();
  MPI_Finalize();
}
