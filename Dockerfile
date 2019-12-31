FROM ubuntu:18.04

RUN apt-get update -y
RUN apt-get install -y cmake git g++ openmpi-bin openmpi-common libopenmpi-dev libopenblas-base libopenblas-dev vim

# Setup home environment
RUN useradd -rm -d /home/ubuntu -s /bin/bash -g root -G sudo -u 1000 ubuntu
USER ubuntu
WORKDIR /home/ubuntu

# Setup FemTech
RUN git clone --single-branch --branch develop https://github.com/PSUCompBio/FemTech
RUN mkdir FemTech/build
RUN cd FemTech/build;cmake .. -DEXAMPLES=ON -DEXAMPLE5=ON; make -j 8;
# RUN cd FemTech/build/examples/ex5; mpirun -mca btl_vader_single_copy_mechanism none ex5 input.json
