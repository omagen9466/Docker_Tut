We start by installing Docker locally on our machine. This will allow to package our athena++ code with all its various dependencies (e.g: python, g++ etc.). The main advantage of this is that we will be able to run it directly on the google cloud server, without having to start building these dependencies/ change operating system on a Virtual Machine (VM) that we will be creating for our purpose. 

Generally speaking, the following steps are taken if we want to run our code on google's VM:

1. Create a Dockerfile locally 
2. Build a docker container 
3. Push the docker image onto google's container repository 
4. Create a VM with the chosen specifications and link the container from the repository to the VM
5. Enter the VM through SSH and run the container
6. Download the outputs to your local machine
Quick and easy!

In the following, I will take you through some necessary installations. Click the links to be sent to the respective tutorial file.
## Installation and other requirements

Lets begin with the [[Docker Installation]]




