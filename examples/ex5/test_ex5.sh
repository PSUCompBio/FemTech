#! /bin/bash
set -eo pipefail

[ "$DEBUG" ] && set -x

# set current working directory to the directory of the script
cd "$(dirname "$0")"

if [ -e output_test.json ]; then
  rm output_test.json
fi
mpirun --oversubscribe -np 8 -mca btl_vader_single_copy_mechanism none ex5 input_test.json
simulationSuccess=$?
if [ $simulationSuccess -eq 0 ]; then
  if [ -e output_test.json ]; then
    echo "Output Generated ... Testing Results"
    max_shear_id=`cat output_test.json| jq -r .'["maximum-shear-strain"]["global-element-id"]'`
    min_principal_id=`cat output_test.json| jq -r .'["principal-min-strain"]["global-element-id"]'`
    max_principal_id=`cat output_test.json| jq -r .'["principal-max-strain"]["global-element-id"]'`
    if [ $max_shear_id != "14196" ] || [ $min_principal_id != "14196" ] || [ $max_principal_id != "14196" ]; then
      echo "Results do not match. Possible issue with code"
      exit 1
    else
      echo "TEST PASSED!"
    fi
  else
    echo "Output not generated...Issue with image" >&2
    exit 1
  fi
else
  echo "MPIRUN failed...Issue with image" >&2
  exit 1
fi
