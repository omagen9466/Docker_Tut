
Download zip from:
https://www.hdfgroup.org/downloads/hdf5/source-code/?1717482886

Follow the steps:

1. Unzip .tar file 

```
tar -xvf <name>.tar --directory <unzip_to_this_path>
```

2. change directory to where you unzipped

4. configure hdf5 to enable parallel runs (if needed, change prefix to compile somewhere else)

```
./configure --enable-parallel --enable-shared --prefix=/usr/local/hdf5-install/
```

4. compile and install

```
make
make install
```

5. If needed, add to .bashrc:
```
export LD_LIBRARY_PATH=<your installation path>/lib:$LD_LIBRARY_PATH
```