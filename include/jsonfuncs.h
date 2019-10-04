#ifndef INCLUDE_JSONFUNCS_H_
#define INCLUDE_JSONFUNCS_H_

#include "json/json.h"

// Function to read JSON input file with confiuration settings
Json::Value getConfig(const char* inputFile);

#endif  // INCLUDE_JSONFUNCS_H_
