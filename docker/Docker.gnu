FROM ubuntu:22.04

# Install System Dependencies
ENV DEBIAN_FRONTEND noninteractive
RUN apt update -y && \
    apt install -y gcc gfortran mpich && \
    apt install -y mpich git vim libxml2-dev curl unzip wget && \
    apt install -y python3 python-is-python3 python3-pip python3-pybind11 && \
    apt install -y m4 libcurl4-openssl-dev zlib1g-dev

ENV CC=gcc
ENV CXX=g++
ENV FC=gfortran

RUN echo "export CPPFLAGS="-I/opt/hdf5/include -I/usr/lib/aarch64-linux-gnu/mpich/include/"" > /root/.bashrc
RUN echo "export LDFLAGS=-L/opt/hdf5/lib" >> /root/.bashrc


# Install HDF5
RUN wget https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-1_14_3.tar.gz && \
    tar xzf hdf5-1_14_3.tar.gz && \
    cd hdf5-hdf5-1_14_3 && \
    CC=mpicc FC=mpifort FCFLAGS=-fPIC CFLAGS=-fPIC \
    ./configure --enable-parallel \
                --enable-fortran \
                --prefix=/opt/hdf5 && \
    make -j$(nproc) install && \
    cd .. && \
    rm -rf hdf5-hdf5-1_14_3 hdf5-1_14_3.tar.gz
ENV LD_LIBRARY_PATH /opt/hdf5/lib:$LD_LIBRARY_PATH
ENV HDF5_ROOT=/opt/hdf5

RUN echo "export CPPFLAGS="-I/opt/hdf5/include -I/usr/lib/aarch64-linux-gnu/mpich/include/"" > /root/.bashrc
RUN echo "export LDFLAGS=-L/opt/hdf5/lib" >> /root/.bashrc
RUN wget -q https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz && \
    tar xf v4.9.2.tar.gz && \
    cd netcdf-c-4.9.2 && \
    CPPFLAGS="-I/opt/zlib/include -I/opt/hdf5/include" LDFLAGS="-L/opt/zlib/lib -L/opt/hdf5/lib"  CC=mpicc ./configure --prefix=/opt/hdf5 --enable-parallel-tests && \
    make -j 2 && make install && \
    cd .. && \
    rm -rf netcdf-c-4.9.2
ENV LD_LIBRARY_PATH /opt/hdf5/lib:$LD_LIBRARY_PATH

RUN wget -q https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.2.tar.gz && \
    tar xf v4.6.2.tar.gz && \
    cd netcdf-fortran-4.6.2 && \
    CPPFLAGS="-I/opt/hdf5/include" LDFLAGS="-L/opt/hdf5/lib"  F77=mpif77 F90=mpif90 FC=mpifort FCFLAGS="-fPIC" FFLAGS="-O3" ./configure --prefix=/opt/hdf5 --enable-parallel && \
    make -j 2 &&  make install && \
    cd .. && \
    rm -rf netcdf-fortran-4.6.2


ENTRYPOINT bash

