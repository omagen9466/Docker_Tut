#!/bin/bash                         #must be at the start of every script
#SBATCH --job-name=jet_4_mcells     #name the job
#SBATCH --nodes=1                   #needed nodes amount
#SBATCH --ntasks-per-node=16        #tasks per node
#SBATCH --ntasks=16                 #total number of tasks
#SBATCH --cpus-per-task=1           #cpus per task
#SBATCH --partition=h3              #VM type (h3 recommended from google for CFD)
#SBATCH --time=00:01:00             #total time needed

##running sim commands:

module load openmpi
mpiexec -n 16 ../athena/bin/athena -i athinput.jet > log