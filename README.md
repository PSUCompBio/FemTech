# Cross-Compiling Test
This is a minimal example of using CMake and C programming that should
be able to be compiled on Windows, Mac and Linux with little or no modification.
##  Windows
#### Software required
- CMAKE gui (https://cmake.org/download/)
- Visual Studio 2017 (community version, free)
- (suggested) atom for ediing source code

#### To compile
- download zip file of project and extract
- use cmake gui
  - point cmake to source code and also specify the build to
    be in the same location as source but in build (i.e., source/build).
  - once configured once specify "Debug" for the CMAKE_CONFIGUATION_TYPES
  - configure and generate makefile and select open project


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
