FROM ubuntu:focal AS buildFemTech
ARG DEBIAN_FRONTEND=noninteractive

ARG TARGETARCH
ARG BRANCH=feature/multiarch

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
COPY --chown=ubuntu:ubuntu ./compileFemTech.sh .
RUN chmod +x compileFemTech.sh
RUN git clone --single-branch --branch $BRANCH https://github.com/PSUCompBio/FemTech
RUN bash compileFemTech.sh $TARGETARCH 
# OpenMPI : RUN cd FemTech/build/examples/ex5; mpirun -mca btl_vader_single_copy_mechanism none ex5 input.json
# Mpich   : RUN cd FemTech/build/examples/ex5; mpirun -np 2 ./ex5 input.json

# To create image : docker build --build-arg BRANCH=CI --pull --cache-from nsfcareer/femtech:develop --tag nsfcareer/femtech:develop --target buildFemTech  -f Dockerfile .

FROM ubuntu:focal
ARG DEBIAN_FRONTEND=noninteractive

ARG TARGETARCH

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

COPY --from=buildFemTech ["/home/ubuntu/FemTech/build_all/.", \
  "/home/ubuntu/FemTechRun/"]

# To create image : docker build --pull --cache-from nsfcareer/femtech:develop --cache-from nsfcareer/femtech:production --tag nsfcareer/femtech:production -f Dockerfile .
# To create image : sudo docker buildx build --platform linux/arm64,linux/amd64 --push -t nsfcareer/femtech:multiarch -f Dockerfile .
