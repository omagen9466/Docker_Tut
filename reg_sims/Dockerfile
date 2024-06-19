# Import the base image
FROM alpine:latest

# Install necessary packages: python3, g++, make, and openmpi (instead of mpich)
RUN apk update && \
    apk add --no-cache python3 py3-pip g++ make openmpi openmpi-dev vim openssh-client

# Ensure the mpi environment is set up
ENV PATH /usr/lib64/openmpi/bin:$PATH
ENV LD_LIBRARY_PATH /usr/lib64/openmpi/lib:$LD_LIBRARY_PATH

ENV HDF5_MAJOR_REL       hdf5-1.12
ENV HDF5_MINOR_REL       hdf5-1.12.0
ENV HDF5_SRC_URL   http://support.hdfgroup.org/ftp/HDF5/releases
RUN mkdir -p /usr/local/src ;\
    cd /tmp                                                                        ; \
    echo "Getting: ${HDF5_SRC_URL}/${HDF5_MAJOR_REL}/${HDF5_MINOR_REL}/src/${HDF5_MINOR_REL}.tar"                ; \
    wget ${HDF5_SRC_URL}/${HDF5_MAJOR_REL}/${HDF5_MINOR_REL}/src/${HDF5_MINOR_REL}.tar                           ; \
    tar -xvf ${HDF5_MINOR_REL}.tar --directory /usr/local/src                      ; \
    rm ${HDF5_MINOR_REL}.tar                                                       ; \
    cd /usr/local/src/${HDF5_MINOR_REL}                                            ; \
    ./configure --enable-parallel --enable-shared --prefix=/usr/local/hdf5-install/                                           ; \
    make                                                                           ; \
#    make check                                                                     ; \
    make install  
# creating athena folder in the contained and cd into it
WORKDIR /usr/app/Sim/athena
# Set up work directory
WORKDIR /usr/app/Sim/athena

# Copy the athena source code from local machine to the container directory
ADD athena .

# If needed, set the default command to run when the container starts (optional)
# CMD ["python3", "your_script.py"]
