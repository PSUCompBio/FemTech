#include <fstream>
#include <stdlib.h>

#include "json/json.h"

// Function to read JSON input file with confiuration settings
Json::Value getConfig(const char* inputFile) {
  Json::Value root;
  std::ifstream ifs;
  ifs.open(inputFile);
  if (!ifs.is_open()) {
    printf("ERROR: Failed to open configuration file\n");
    exit(1);
  }

  Json::CharReaderBuilder builder;
  JSONCPP_STRING errs;
  if (!parseFromStream(builder, ifs, &root, &errs)) {
    printf("ERROR : %s\n", errs.c_str());
    exit(1);
  }
  return root["simulation"];
}
