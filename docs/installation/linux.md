---
description: How to build FemTech on Linux.
---

# Linux

 Open terminal \(navigate to desired directory\) and run commands outlined below.

The flow below was completed on ACI-I. 

```text
> git clone https://github.com/PSUCompBio/FemTech
> cd FemTech
> mkdir build
(on PSU ACI-I) > module load cmake gcc/5.3.1 openmpi/1.10.1 blas/3.6.0 lapack/3.6.0
(on Amazon E2C) ...TBD
> cd bulid
> ccmake ..
> make -j8
(for interactive parallel on ACI-B) qsub -A open -l walltime=1:00:00 -l nodes=1:ppn=2 -I
```

Here is video of a build on Penn State's Advanced Cyber Infrastructure system \(ACI-I\)

{% embed url="https://youtu.be/\_y0PvQxQSak" %}



