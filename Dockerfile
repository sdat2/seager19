FROM ubuntu:16.04 

# first half also works with ubuntu:14.04
# I pieced varioud different docker files together so that they could reproduce the environment.

# https://github.com/tianon/docker-brew-ubuntu-core/blob/c5bc8f61f0e0a8aa3780a8dc3a09ae6558693117/trusty/Dockerfile

# sudo apt-get update
# sudo apt-get install --fix-missing libnetcdf-dev libnetcdff-dev
# sudo apt-get install cloc

# install sudo
RUN apt-get update && apt-get -yq install sudo && \
    sudo apt-get -yq install git curl && \
    sudo apt-get -yq --fix-missing install make cmake gfortran-4.8 gcc-4.8 cloc  && \
    sudo apt-get -yq install --fix-missing libnetcdf-dev libnetcdff-dev && \
    sudo apt-get clean -q &&\
    sudo apt-get -yq install ncurses-dev && \
    cd /tmp \ 
    && curl -O https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf36_4/linux/cdf36_4-dist-all.tar.gz \ 
    && tar xzf cdf36_4-dist-all.tar.gz \ 
    && cd cdf36_4-dist \ 
    && make OS=linux ENV=gnu CURSES=yes FORTRAN=no UCOPTIONS=-O2 SHARED=yes all \
    && sudo make INSTALLDIR=/usr/local/cdf install

# x86_64-linux-gnu-gcc-4.8
# x86_64-linux-gnu-gfortran-4.8

# libcoarrays-dev libopenmpi-dev open-coarrays-bin

# https://urban-institute.medium.com/fortran-and-docker-how-to-combine-legacy-code-with-cutting-edge-components-35e934b15023

# JAK HelloWorld_Fortran test

# start by building the basic container
# FROM centos:latest
# MAINTAINER Jessica Kelly <jkelly@urban.org>
# RUN yum update -y
# add gfortran, debugging tools and make
# RUN yum install -y gcc-gfortran gdb make

# build the hello world code
# COPY Makefile HelloWorld.f90 /fortran/
# WORKDIR /fortran/
# RUN make HelloWorld

# configure the container to run the hello world executable by default
# CMD ["./HelloWorld"]

################ Anaconda section. ####################
# https://github.com/ContinuumIO/docker-images/blob/master/anaconda3/debian/Dockerfile


ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

# hadolint ignore=DL3008
RUN apt-get update --fix-missing && \
    sudo apt-get install -y --no-install-recommends wget bzip2 ca-certificates && \
    # build-essential
    # libglib2.0-0 libxext6 libsm6 libxrender1
    # apt-get install -y --no-install-recommends wget build-essential bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion sudo && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    wget --quiet https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

# ARG PYTHON_VERSION=3.8

CMD [ "/bin/bash" ]
