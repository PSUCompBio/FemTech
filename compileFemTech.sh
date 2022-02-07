#!/bin/bash

# Script to build for different architectures on docker
# If no argument provided default to amd64
if [ -z "$1" ]; then
  PARCH="amd64"
else
  PARCH=$1
fi

if [ "$PARCH" == "amd64" ]; then
  mkdir FemTech/build_c6i FemTech/build_c5 FemTech/build_c5a FemTech/build_native FemTech/build_all
  # Build for c6i, icelake 
  cd FemTech/build_c6i
  cmake .. -DPROC_ARCH=icelake-server -DEXAMPLES=ON -DEXAMPLE5=ON -DEXAMPLE21=ON && make -j 8
  cd ../build_c5
  # Build for c5, cascadelake 
  cmake .. -DPROC_ARCH=cascadelake -DEXAMPLES=ON -DEXAMPLE5=ON -DEXAMPLE21=ON && make -j 8
  cd ../build_c5a
  # Build for c5a, amd epyc 2 
  cmake .. -DPROC_ARCH=znver2 -DEXAMPLES=ON -DEXAMPLE5=ON -DEXAMPLE21=ON && make -j 8
  cd ../build_native
  # Build for travis to test container
  cmake .. -DPROC_ARCH=native -DEXAMPLES=ON -DEXAMPLE5=ON -DEXAMPLE21=ON && make -j 8
  ### Copy all executables and required input files to a single folder #####
  cd ../build_all
  ## Copy files required for displacement BC ##
  cp ../build_native/examples/ex5/{input.json,materials.dat,coarse_brain.inp,simulationMovie.py,mps95Movie.py,addGraph.py,fine_cellcentres.txt,coarse_cellcentres.txt,updateOutputJson.py,fine_brain.inp} .
  cp ../build_c6i/examples/ex5/ex5 ex5_c6i
  cp ../build_c5/examples/ex5/ex5 ex5_c5
  cp ../build_c5a/examples/ex5/ex5 ex5_c5a
  ## Copy files required for pressure/traction BC ##
  cp ../build_c6i/examples/ex21/materials.dat materialsPressure.dat
  cp ../build_c6i/examples/ex21/ex21 ex21_c6i
  cp ../build_c5/examples/ex21/ex21 ex21_c5
  cp ../build_c5a/examples/ex21/ex21 ex21_c5a
fi

if [ "$PARCH" == "arm64" ]; then
  mkdir FemTech/build_c6g FemTech/build_all
  cd FemTech/build_c6g
  # Build for c6g, aws graviton 2 
  cmake .. -DPROC_ARCH=arm64 -DEXAMPLES=ON -DEXAMPLE5=ON -DEXAMPLE21=ON && make -j 8
  ### Copy all input and executables to required folder ###
  cd ../build_all
  cp ../build_c6g/examples/ex5/{input.json,materials.dat,coarse_brain.inp,simulationMovie.py,mps95Movie.py,addGraph.py,fine_cellcentres.txt,coarse_cellcentres.txt,updateOutputJson.py,fine_brain.inp} .
  cp ../build_c6g/examples/ex5/ex5 ex5_c6g
  ### Copy pressure/traction bc executable ##
  cp ../build_c6g/examples/ex21/ex21 ex21_c6g
  cp ../build_c6g/examples/ex21/materials.dat materialsPressure.dat
fi
