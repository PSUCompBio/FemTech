FROM ubuntu:18.04 AS buildFemTech

# using mpich because of bug with openmpi
# https://bugs.launchpad.net/ubuntu/+source/openmpi/+bug/1838684
RUN apt-get update && \
  apt-get install -y  --no-install-recommends \
  cmake git g++ mpich libmpich-dev ca-certificates\
  libopenblas-base libopenblas-dev vim make\
  && rm -rf /var/lib/apt/lists/*

# Setup home environment
RUN useradd -rm -d /home/ubuntu -s /bin/bash -g root -G sudo -u 1000 ubuntu
USER ubuntu
WORKDIR /home/ubuntu

# Setup FemTech
ADD https://api.github.com/repos/PSUCompBio/FemTech/git/refs/heads/develop version.json
RUN git clone --single-branch --branch develop https://github.com/PSUCompBio/FemTech
RUN mkdir FemTech/build
RUN cd FemTech/build;cmake .. -DEXAMPLES=ON -DEXAMPLE5=ON; make -j 8;
# OpenMPI : RUN cd FemTech/build/examples/ex5; mpirun -mca btl_vader_single_copy_mechanism none ex5 input.json
# Mpich   : RUN cd FemTech/build/examples/ex5; mpirun -np 2 ./ex5 input.json

# To create image : docker build --pull --cache-from nsfcareer/femtech:develop --tag nsfcareer/femtech:develop --target buildFemTech  -f Dockerfile .

FROM ubuntu:18.04

RUN apt-get update && \
  apt-get install -y  --no-install-recommends \
  mpich libopenblas-base \
  && rm -rf /var/lib/apt/lists/*

# Setup home environment
RUN useradd -rm -d /home/ubuntu -s /bin/bash -g root -G sudo -u 1000 ubuntu
USER ubuntu
WORKDIR /home/ubuntu

# Setup FemTech
RUN mkdir FemTechRun FemTechRun/results FemTechRun/results/vtu

COPY --from=buildFemTech ["/home/ubuntu/FemTech/build/examples/ex5/ex5", \
  "/home/ubuntu/FemTech/build/examples/ex5/input.json", \
  "/home/ubuntu/FemTech/build/examples/ex5/materials.dat", \
  "/home/ubuntu/FemTech/build/examples/ex5/coarse_brain.inp", \
  "/home/ubuntu/FemTechRun/"]

# To create image : docker build --pull --cache-from nsfcareer/femtech:develop --cache-from nsfcareer/femtech:production --tag nsfcareer/femtech:production -f Dockerfile .
