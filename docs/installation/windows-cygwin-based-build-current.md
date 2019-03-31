---
description: Explains how to build using Cygwin...
---

# Windows - Cygwin Based Build \(Current\)

The current suggested setup to work on Windows is to use Atom \(or Vi for advanced linux users\) for editing source code and Cygwin terminal for building \(a.k.a compiling\). 

## Step-by-step Instructions

1. Download and install Atom: [https://atom.io/](https://atom.io/)
2. Download and install Cygwin: [https://cygwin.com/](https://cygwin.com/). Use Cygwin installer to install the following packages:
   1. vim 8.0
   2. vim-common 8.0
   3. fontconfig \(??\)
   4. make 4.2.1-2
   5. cmake 3.6.2-1
   6. mpfr
   7. libmpfr4 3.1.6
   8. gcc-core 7.4
   9. gcc-fortran 7.4
   10. gcc-g++ 7.4
   11. libopenmpi12 1.10.7-1        
   12. libopenmpi-devel 3.1.1-2  
   13. openmpi 3.1.2-2 
   14. liblapack-devel 3.8.0-1
   15. libopenblas 0.3.5-1 
3. Download and install Paraview: [https://www.paraview.org/](https://www.paraview.org/)
4. **Build FemTech**:  Open cygwin terminal Navigate to a directory where you would like to install FemTech or make one, e.g.  mkdir code cd code git clone [https://github.com/PSUCompBio/FemTech](https://github.com/PSUCompBio/FemTech) cd FemTech mkdir build cd build ccmake .. In the ccmake gui, add "Enable MPI" and "Examples", then configure In the ccmake gui, select the examples you would like to turn on \(see list here\) In the ccmake gui, configure again, then generate \(g\) to create the Makefile. make -j 8 \(where 8 is the number of processors you have\)



