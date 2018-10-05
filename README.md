# Cross-Compiling Test
This is a minimal example of using CMake and C programming that should
be able to be compiled on Windows, Mac and Linux with little or no modification.
##  Windows
#### Software required
- CMake-GUI (https://cmake.org/download/)
- Visual Studio 2017 (community version, free with Microsoft MPI option, see https://www.microsoft.com/en-us/download/details.aspx?id=56727)


#### To compile
- download zip file of project and extract 
- use CMake-GUIpoint CMake to digitalbrain source code and also specify the build to be in the same location as source but in build (i.e., digitalbrain/build)
- select "configure"
- For parallel Enable MPI and Examples then select example to turn on, e.g. ex1 
- configure again and generate makefile and select open project
- once Visual Studio opens up, right click on ex1 "solution"(on right hand side) and select "Set as Startup Project"
- then right click again the ex1 solution again and select Properties
- Under the Configuration Properties-> Debugging, set "Command" to C:\Program Files\Microsoft MPI\Bin\mpiexec.exe
- Under the Configuration Properties-> Debugging, set "Commnad Arguments" to -n 2 "$(TargetPath)" 1-elt-cube.k (where the .k file is the input file in the ex1 directory)
- Then go to top menu bar and select Debug->Start Without Debugging

## Linux
#### To compile
- open terminal (navigate to desired directory)
```
git clone https://github.com/PSUCompBio/cross-compiling-test
cd cross-compiling-test
mkdir build
cd bulid
ccmake ..
make -j8
```

## Mac OS
#### To compile
- open terminal (navigate to desired directory)
```
git clone https://github.com/PSUCompBio/cross-compiling-test
cd cross-compiling-test
mkdir build
cd bulid
ccmake ..
make -j8
```
