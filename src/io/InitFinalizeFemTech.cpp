#include "FemTech.h"
#include "gitbranch.h"

#include "jsonfuncs.h"

#include "mpi.h"
#include <string>
#include <ctime>
#include <stdio.h>

int world_rank;
int world_size;

Json::Value InitFemTech(int argc, char **argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (argc < 1) {
    fprintf(stdout, "ERROR(%d): Please provide an input file\n", world_rank);
    fprintf(stdout, "ERROR(%d): Usage %s inputFile\n", world_rank, argv[0]);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  // Read the input file
  Json::Value inputJson = getConfig(argv[1]);
  // create simulation unique id from time
  // Check if UID is present in input JSON
  std::string uid;
  if (inputJson["uid"].empty()) {
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    char buffer [80];
    strftime(buffer, 80, "%d_%m_%Y_%H_%M_%S", now);
    uid = buffer;
    inputJson["uid"] = uid;
  } else {
    uid = inputJson["uid"].asString();
  }
  std::string logFile = "femtech_"+uid+".log";

  // Initialize the output log file
  initLog(logFile.c_str());
  FILE_LOG_MASTER(INFO, "Code with commit hash : %s of branch %s", GIT_COMMIT_HASH,
           GIT_BRANCH);
  return inputJson;
}

void InitFemTechWoInput(int argc, char **argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  // Get the rank of the process
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (argc < 1) {
    fprintf(stdout, "ERROR(%d): Please provide an input file\n", world_rank);
    fprintf(stdout, "ERROR(%d): Usage %s inputFile\n", world_rank, argv[0]);
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  // create simulation unique id from time
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  char buffer [80];
  strftime(buffer, 80, "%d_%m_%Y_%H_%M_%S", now);
  std::string uid = buffer;
  std::string logFile = "femtech_"+uid+".log";

  // Initialize the output log file
  initLog(logFile.c_str());
  FILE_LOG_MASTER(INFO, "Code with commit hash : %s of branch %s", GIT_COMMIT_HASH,
           GIT_BRANCH);
}

void FinalizeFemTech() {
  FILE_LOG_MASTER(INFO, "Terminating Application Normally");
  finaliseLog();
  FreeArrays();
  MPI_Finalize();
}

void TerminateFemTech(int errorCode) {
  finaliseLog();
  FreeArrays();
  MPI_Abort(MPI_COMM_WORLD, errorCode);
}
