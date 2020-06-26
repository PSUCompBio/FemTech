#! /bin/bash
set -eo pipefail

[ "$DEBUG" ] && set -x

# set current working directory to the directory of the script
cd "$(dirname "$0")"

# input_test_?.json array
testInputFiles=( p1 pt x s )
# Values of element id expected in the order, max-shear, min-principal,
# max-principal
value_p1=( 4185 4185 4185 )
value_pt=( 4749 16263 4749 )
value_s=( 5299 5299 5167 )
value_x=( 5439 5439 6963 )

# Run all tests and check output
for name in "${testInputFiles[@]}"
do
  input_filename="input_test_${name}.json"
	echo "Input file name : ${input_filename}"
  uid=`cat ${input_filename}|jq -r .'["uid"]'`
  output_filename="output_${uid}.json"
  eval expectedVal=( \${value_$name[@]} )
  # Delete output file if it exists
  if [ -e $output_filename ]; then
    rm $output_filename
  fi
  mpirun --oversubscribe -np 8 -mca btl_vader_single_copy_mechanism none ex5 $input_filename
  simulationSuccess=$?
  if [ $simulationSuccess -eq 0 ]; then
    if [ -e $output_filename ]; then
      echo "  Output Generated ... Testing Results"
      max_shear_id=`cat ${output_filename}| jq -r .'["maximum-shear-strain"]["global-element-id"]'`
      min_principal_id=`cat ${output_filename}| jq -r .'["principal-min-strain"]["global-element-id"]'`
      max_principal_id=`cat ${output_filename}| jq -r .'["principal-max-strain"]["global-element-id"]'`
      if [ $max_shear_id != ${expectedVal[0]} ] || [ $min_principal_id != ${expectedVal[1]} ] || [ $max_principal_id != ${expectedVal[2]} ]; then
        echo "  Results do not match. Possible issue with code"
        exit 1
      else
        echo "  TEST PASSED!"
      fi
    else
      echo "  Output not generated...Issue with image" >&2
      exit 1
    fi
  else
    echo "  MPIRUN failed...Issue with image" >&2
    exit 1
  fi
done
exit 0
