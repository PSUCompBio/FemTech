#!/bin/bash
#
#UNAMEX="ubuntu"
#cd /home/$UNAMEX
#sudo apt-get update
#sudo apt-get install -y build-essential
#sudo apt-get install -y cmake-curses-gui
#sudo apt-get install -y libblas-dev liblapack-dev
## sudo apt-get install -y openmpi-bin openmpi-common openssh-client openssh-server libopenmpi1.10
#sudo apt-get install -y openmpi-bin openmpi-common libopenmpi-dev
#git clone https://github.com/PSUCompBio/FemTech
#cd FemTech
mkdir build
cd build
cmake .. -DEXAMPLES=ON -DEXAMPLE12=ON -DEXAMPLE11=ON
make -j8

#make -j8 ex12 In the examples/ex12/ there is also a makefile. use it

#mpirun -n 1 ex11 1-elt-cube.k > debug_log.txt
#mpirun -n 1 ex12 1-elt-truss.k > debug_log.txt
#paraview 1-elt-cube.pvd  here the *.pvd file
