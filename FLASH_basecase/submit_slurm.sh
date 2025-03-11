#!/bin/bash
#SBATCH --job-name=flashLasersJet    # Job name
#SBATCH --time=96:00:00               # Time limit hrs:min:sec
#SBATCH -p defq
#SBATCH -N 1                   # Run on a single CPU
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=a.hirsch@hzdr.de     # Where to send mail    	

echo 'start...'
mpirun -np 1 flash4
#mpiexec -n 4 flash4
echo 'done'
