---
description: How to build FemTech on Linux.
---

# Linux

 open terminal \(navigate to desired directory\)

```text
git clone https://github.com/PSUCompBio/FemTech
cd FemTech
mkdir build
(on PSU ACI-I) module load cmake gcc/5.3.1 openmpi/1.10.1
(on Amazon E2C) ...TBD
cd bulid
ccmake ..
make -j8
(for interactive parallel on ACI-B) qsub -A open -l walltime=1:00:00 -l nodes=1:ppn=2 -I
```

