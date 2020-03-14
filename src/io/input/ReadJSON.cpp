#include <fstream>

#include "json/json.h"

#include "FemTech.h"
#include "parallel_log.h"

// Function to read JSON input file with confiuration settings
Json::Value getConfig(const char* inputFile) {
  Json::Value root;
  std::ifstream ifs;
  ifs.open(inputFile);
  if (!ifs.is_open()) {
    fprintf(stdout, "Failed to open configuration file");
    MPI_Abort(MPI_COMM_WORLD, 3);
  }

  Json::CharReaderBuilder builder;
  JSONCPP_STRING errs;
  if (!parseFromStream(builder, ifs, &root, &errs)) {
    fprintf(stdout, "%s", errs.c_str());
    MPI_Abort(MPI_COMM_WORLD, 3);
  }
  return root;
}

void jsonToArray(double* array, const Json::Value& jsonArray) {
  int size = jsonArray.size();
  for (int i = 0; i < size; ++i) {
    array[i] = jsonArray[i].asDouble();
  }
}
