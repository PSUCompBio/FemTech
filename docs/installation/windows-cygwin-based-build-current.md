---
description: Explains how to build using Cygwin...
---

# Windows - Cygwin Based Build \(Current\)

The current suggested setup to work on Windows is to use Atom for editing source code and Cygwin terminal for building \(a.k.a compiling\). 

## Step-by-step Instructions

1. Download and install Atom: [https://atom.io/](https://atom.io/)
2. Download and install Cygwin: [https://cygwin.com/](https://cygwin.com/) 
   1. Here you will need to use the installer to select a few different packages. I am able to provide a full list of the package I currently have - it is long and it is possible to go down the list one-by-one to check that you are installing it but, a good first step might be to install a few big ones I remember and then go from there. Below is a full list \(text file\), but for now just try the ones I can remember:  
      1. MPI, including: 
         1. vim 8.0
         2. vim-common 8.0
         3. fontconfig \(??\)
         4. make 4.2.1-2
         5. cmake 3.6.2-1
         6. mpfr 4.0.2-1
         7. libmpfr4 3.1.6
         8. gcc-core 7.4
         9. gcc-fortran 7.4
         10. gcc-g++ 7.4
         11. libopenmpi12 1.10.7-1        
         12. libopenmpi-devel 3.1.1-2  
         13. openmpi 1.10.7-1  


             additional there are:  


             libopenmpicxx1  -   C++ library  
      2. cmake 3.6.2 
      3. gcc \(including gfortran\). 
      4. lapack and blas 
3. Download and install Paraview: [https://www.paraview.org/](https://www.paraview.org/)
4. add entry
5. add entry
6. add entry 

{% file src="../.gitbook/assets/cygwin-configuration.bin" caption="List of Cygwin install packages" %}

