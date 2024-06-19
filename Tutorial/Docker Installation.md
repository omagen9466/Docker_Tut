Full guide on: https://docs.docker.com/desktop/install/linux-install/#system-requirements



> [!faq] Prerequisites
> **Supported Platforms**
> - Ubuntu
> - Debian
> - Red Hat Enterprise Linux
> - Fedora
> 
> **Further requirements**
> - 64-bit CPU with Virtualization Support enabled.
> - At least 4GB RAM
> - A GUI desktop environment (Preferably GNOME, MATE, or KDE )
> - A [Sudo User](https://www.linuxtechi.com/how-to-create-sudo-user-on-ubuntu/) with admin rights


## KVM Virtulisation Support

The `kvm` module should load automatically if the host has virtualization support. To load the module manually, run:


```
madprobe kvm
```

Depending on the processor of the host machine, the corresponding module must be loaded:

```
modprobe kvm_intel  # Intel processors  
modprobe kvm_amd    # AMD processors
```


1. Gnome-terminal needs to be installed for situations without the Gnome Desktop: 

```
sudo apt install gnome-terminal
```

2. Uninstall the tech preview or beta version of Docker Desktop for Linux. Run:

```
sudo apt remove docker-desktop
```

## Installation 

Now install the Docker. But before that update package lists and install the requisite dependencies as follows.

```
sudo apt update  
sudo apt install software-properties-common curl apt-transport-https ca-certificates -y
```

Once the installation is complete, add Docker’s GPG signing key.

```
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/trusted.gpg.d/docker-archive-keyring.gpg
```

Next, add the official Docker’s repository to your system as follows.

```
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
```

Once the repository is in place, install Docker and other Docker tools as shown.

```
sudo apt install docker-ce docker-ce-cli containerd.io uidmap -y
```

To verify that docker is running, execute the following command:

```
sudo systemctl status docker
```


> [!tip] The following should be seen in the terminal:
> 
> (base) admina@om-pc:~/Documents/Docker_Tut$ sudo systemctl status docker
> [sudo] password for admina: 
> ● docker.service - Docker Application Container Engine
>      Loaded: loaded (/lib/systemd/system/docker.service; enabled; vendor preset: enabled)
>      Active: active (running) since Wed 2024-04-20 08:47:16 IDT; 1h 19min ago
> TriggeredBy: ● docker.socket
>        Docs: https://docs.docker.com
>    Main PID: 1652 (dockerd)
>       Tasks: 24
>      Memory: 108.3M
>         CPU: 1.796s
>      CGroup: /system.slice/docker.service
>              └─1652 /usr/bin/dockerd -H fd:// --containerd=/run/containerd/containerd.sock

You can also check your version:

```
docker version
```

**Install Docker Desktop on Ubuntu 22.04:**

Install the docker desktop using below wget command. The latest version of Docker Desktop is Docker Desktop version 4.19.0.

```
wget https://desktop.docker.com/linux/main/amd64/docker-desktop-4.19.0-amd64.deb
```

Also you can download the DEB package in this [link](https://docs.docker.com/desktop/install/ubuntu/) and use the below command for installation.Also you can download the DEB package in this link and use the below command for installation.

```
sudo apt install ./docker-desktop-*-amd64.deb
```

**You should now be able to start docker from your desktop!**

