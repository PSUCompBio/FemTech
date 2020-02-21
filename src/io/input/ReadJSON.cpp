#include <fstream>
#include <stdlib.h>

#include "json/json.h"

#include "parallel_log.h"

// Function to read JSON input file with confiuration settings
Json::Value getConfig(const char* inputFile) {
  Json::Value root;
  std::ifstream ifs;
  ifs.open(inputFile);
  if (!ifs.is_open()) {
    FILE_LOG_SINGLE(ERROR, "Failed to open configuration file");
    exit(1);
  }

  Json::CharReaderBuilder builder;
  JSONCPP_STRING errs;
  if (!parseFromStream(builder, ifs, &root, &errs)) {
    FILE_LOG_SINGLE(ERROR, "%s", errs.c_str());
    exit(1);
  }
  return root["simulation"];
}

void jsonToArray(double* array, const Json::Value& jsonArray) {
  int size = jsonArray.size();
  for (int i = 0; i < size; ++i) {
    array[i] = jsonArray[i].asDouble();
  }
}
