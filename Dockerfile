FROM ubuntu:18.04 AS buildFemTech

ARG BRANCH=develop
ARG PARCH=native

RUN apt-get update && \
  apt-get install -y  --no-install-recommends \
  cmake git g++ openmpi-bin openmpi-common ca-certificates \
  libopenmpi-dev libopenblas-base libopenblas-dev vim make \
  openssh-client jq \
  && rm -rf /var/lib/apt/lists/*

# Setup home environment
RUN useradd -rm -d /home/ubuntu -s /bin/bash -g root -G sudo -u 1000 ubuntu
USER ubuntu
WORKDIR /home/ubuntu

# Setup FemTech
RUN git clone --single-branch --branch $BRANCH https://github.com/PSUCompBio/FemTech
RUN mkdir FemTech/build_native
RUN cd FemTech/build_native;cmake .. -DEXAMPLES=ON -DEXAMPLE5=ON; make -j 8;
RUN if [ "$PARCH" != "native" ]; then \
    mkdir FemTech/build_$PARCH && \
    cd FemTech/build_$PARCH && \
    cmake .. -DPROC_ARCH=$PARCH -DEXAMPLES=ON -DEXAMPLE5=ON && \
    make -j 8; \ 
  fi

# OpenMPI : RUN cd FemTech/build/examples/ex5; mpirun -mca btl_vader_single_copy_mechanism none ex5 input.json
# Mpich   : RUN cd FemTech/build/examples/ex5; mpirun -np 2 ./ex5 input.json

# To create image : docker build --build-arg BRANCH=CI --pull --cache-from nsfcareer/femtech:develop --tag nsfcareer/femtech:develop --target buildFemTech  -f Dockerfile .

FROM ubuntu:18.04

ARG PARCH=native

RUN apt-get update && \
  apt-get install -y  --no-install-recommends \
  openmpi-bin libopenblas-base openssh-client \
  && rm -rf /var/lib/apt/lists/*

# Setup home environment
RUN useradd -rm -d /home/ubuntu -s /bin/bash -g root -G sudo -u 1000 ubuntu
USER ubuntu
WORKDIR /home/ubuntu

# Setup FemTech
RUN mkdir FemTechRun FemTechRun/results FemTechRun/results/vtu

COPY --from=buildFemTech ["/home/ubuntu/FemTech/build_${PARCH}/examples/ex5/ex5", \
  "/home/ubuntu/FemTech/build_${PARCH}/examples/ex5/input.json", \
  "/home/ubuntu/FemTech/build_${PARCH}/examples/ex5/materials.dat", \
  "/home/ubuntu/FemTech/build_${PARCH}/examples/ex5/coarse_brain.inp", \
  "/home/ubuntu/FemTech/build_${PARCH}/examples/ex5/simulationMovie.py", \
  "/home/ubuntu/FemTech/build_${PARCH}/examples/ex5/mps95Movie.py", \
  "/home/ubuntu/FemTech/build_${PARCH}/examples/ex5/addGraph.py", \
  "/home/ubuntu/FemTech/build_${PARCH}/examples/ex5/fine_cellcentres.txt", \
  "/home/ubuntu/FemTech/build_${PARCH}/examples/ex5/coarse_cellcentres.txt", \
  "/home/ubuntu/FemTech/build_${PARCH}/examples/ex5/updateOutputJson.py", \
  "/home/ubuntu/FemTechRun/"]

# To create image : docker build --pull --cache-from nsfcareer/femtech:develop --cache-from nsfcareer/femtech:production --tag nsfcareer/femtech:production -f Dockerfile .
