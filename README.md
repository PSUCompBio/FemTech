# Cross-Compiling Test
This is a minimal example of using CMake and C programming that should
be able to be compiled on Windows, Mac and Linux with little or no modification.
##  Windows
#### Software required
- CMake-GUI (https://cmake.org/download/)
- Visual Studio 2017 (community version, free with Microsoft MPI option, see https://www.microsoft.com/en-us/download/details.aspx?id=56727)


#### To compile
- download zip file of project and extract
- use CMake-GUI point CMake to digitalbrain source code and also specify the build to be in the same location as source but in build (i.e., digitalbrain/build)
- select "configure"
- For parallel:
	- Enable MPI and Examples then select example to turn on, e.g. ex1
	- configure again and generate makefile and select open project
	- Once Visual Studio opens up, right click on ex1 "solution" (on right hand side) and select "Set as Startup Project"
	- Then right click again the ex1 solution again and select Properties
	- Under the Configuration Properties-> Debugging, set "Command" to C:\Program Files\Microsoft MPI\Bin\mpiexec.exe
	- Under the Configuration Properties-> Debugging, set "Command Arguments" to -n 2 "$(TargetPath)" 1-elt-cube.k (where the .k file is the input file in the ex1 directory)
	- Then go to top menu bar and select Debug->Start Without Debugging
- For Serial:
	- Configure using defaults, but turn on examples (for at least example 1)
	- Configure again and generate makefile and select open project
	- Once Visual Studio opens up, right click on ex1 "solution" (on right hand side) and select "Set as Startup Project"
	- Then right click again the ex1 solution again and select Properties
	- Under the Configuration Properties-> Debugging, set "Command Arguments" to 1-elt-cube.k (where the .k file is the input file in the ex1 directory)
	- Then go to top menu bar and select Debug->Start Without Debugging

## Linux
#### To compile
- open terminal (navigate to desired directory)
```
git clone https://github.com/PSUCompBio/digitalbrain
cd digitalbrain
mkdir build
(on PSU ACI-I) module load cmake gcc/5.3.1 openmpi/1.10.1
(on Amazon E2C) ...TBD
cd bulid
ccmake ..
make -j8
(for parallel on ACI-B) qsub -A open -l walltime=1:00:00 -l nodes=1:ppn=2 -I
```

## Mac OS
#### Software required
- Xcode with developer command line tools (https://www.embarcadero.com/starthere/berlin/mobdevsetup/ios/en/installing_the_xcode_command_line_tools_on_a_mac.html)
- CMake (https://cmake.org/download/)
#### To compile
- open terminal (navigate to desired directory)
```
git clone https://github.com/PSUCompBio/digitalbrain
cd cross-compiling-test
mkdir build
cd bulid
ccmake ..
make -j8
```
