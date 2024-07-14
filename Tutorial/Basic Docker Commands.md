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

FOR DAUDE:

```
python3 configure.py --prob water_air_shock_mixure  -mpi -hdf5 -h5double --hdf5_path /home/admina/Software  --nghost 6 --nscalars 5 --nprimitive_scalars 5 -mixture_model  -multi_materials 
```

```
python3 configure.py --prob water_air_shock_mixurev2  -mpi -hdf5 -h5double --hdf5_path /home/admina/Software  --nghost 6 --nscalars 5 --nprimitive_scalars 5  -multi_materials --eos general/ideal
```

