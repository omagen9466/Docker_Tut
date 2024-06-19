# import base image: as athena runs locally on ubuntu 22.04,
# im importing it as my base image, however other base images can be used. see: https://hub.docker.com/ 
FROM ubuntu:22.04 
# updating, installing python3 for configuring and g++ for compiling the athena code
RUN apt update
RUN apt install python3 -y
RUN apt-get install g++ -y
RUN apt-get install make -y
RUN apt install mpich -y
RUN apt-get install wget -y

ENV HDF5_MAJOR_REL       hdf5-1.12
ENV HDF5_MINOR_REL       hdf5-1.12.0
ENV HDF5_SRC_URL   http://support.hdfgroup.org/ftp/HDF5/releases
RUN cd /tmp                                                                        ; \
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

# cp the athena code from local machine to the container dir. NOTE: athena source is located in same dir as the Dockerfile 
ADD athena_version4 .

# # configuring and compiling. NOTE: "space" input translates into separate strings as seen bellow 
# CMD ["python3","configure.py","--prob","jet","-mpi"]
# CMD ["make","clean"]
# CMD ["make"]

# # switching back to where i want to run the sim
# WORKDIR /usr/app/Sim/


# #define volume mount
# VOLUME /usr/app/Sim/data

# # switching back to where I want the output to be in 
# WORKDIR /usr/app/Sim/data/MPI_docker
 
# # starting the sim
# # CMD ["../athena/bin/athena","-i","../athinput.jet",">","log.txt"]
# CMD ["mpiexec", "-n","4", "../../athena/bin/athena -i athinput_par_vtk.jet > log"]

