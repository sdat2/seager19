FROM ubuntu:16.04

# https://github.com/tianon/docker-brew-ubuntu-core/blob/c5bc8f61f0e0a8aa3780a8dc3a09ae6558693117/trusty/Dockerfile

# sudo apt-get update
# sudo apt-get install --fix-missing libnetcdf-dev libnetcdff-dev
# sudo apt-get install cloc
ENV USERNAME=backus

# install sudo
RUN apt-get -yq update && apt-get -yq install sudo

# create and switch to a user
RUN echo "backus ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
RUN useradd --no-log-init --home-dir /home/$USERNAME --create-home --shell /bin/bash $USERNAME
RUN adduser $USERNAME sudo
USER $USERNAME
WORKDIR /home/$USERNAME

# install packages
RUN sudo apt-get install -yq git curl && \
    sudo apt-get install --no-install-recommends -yq make cmake gfortran libcoarrays-dev libopenmpi-dev open-coarrays-bin && \
    sudo apt-get clean -q

# get modern-fortran code
RUN git clone https://github.com/modern-fortran/tsunami
RUN git clone https://github.com/modern-fortran/stock-prices
RUN git clone https://github.com/modern-fortran/weather-buoys
RUN git clone https://github.com/modern-fortran/generic-procedures
RUN git clone https://github.com/modern-fortran/countdown
RUN git clone https://github.com/modern-fortran/tcp-client-server
RUN git clone https://github.com/modern-fortran/listings

# extras
RUN git clone https://github.com/modern-fortran/neural-fortran
RUN git clone https://github.com/wavebitscientific/datetime-fortran
RUN git clone https://github.com/wavebitscientific/functional-fortran

# https://urban-institute.medium.com/fortran-and-docker-how-to-combine-legacy-code-with-cutting-edge-components-35e934b15023

# JAK HelloWorld_Fortran test

# start by building the basic container
FROM centos:latest
# MAINTAINER Jessica Kelly <jkelly@urban.org>
RUN yum update -y
# add gfortran, debugging tools and make
RUN yum install -y gcc-gfortran gdb make 

# build the hello world code
COPY Makefile HelloWorld.f90 /fortran/
WORKDIR /fortran/
RUN make HelloWorld

# configure the container to run the hello world executable by default
CMD ["./HelloWorld"]
