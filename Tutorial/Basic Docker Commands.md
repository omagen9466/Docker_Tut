- To see which docker version is currently installed 

```
docker version
```

- To verify that docker is running

```
sudo systemctl start docker
```

- To verify that docker is running

```
sudo systemctl status docker
```

- Build a Dockerfile

```
docker build -t <usr defined docker name> .
```

- Run the Docker image

```
docker run <usr defined docker name>
```

- Run the Docker image manually

```
docker run -it conf_test /bin/bash
```

- Run Docker image and allocate volume into local machine

```
docker run -v ~/Documents/Docker_Tut/data:/usr/app/Sim/data <usr defined docker name>
```

- configure 
```
python3 configure.py --prob jet -mpi -hdf5 --hdf5_path=/home/admina/Software
```
- export lib path
```
export LD_LIBRARY_PATH=/usr/local/hdf5-install/lib/:$LD_LIBRARY_PATH
```