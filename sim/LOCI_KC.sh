#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem 100G            # memory pool per process
#SBATCH -o LOCI_KC.out    # STDOUT
#SBATCH -t 12:00:00            # time (D-HH:MM)

source ~/py312/bin/activate
python 1_sim_ld.py --bfile LOCI --props 0.001 0.01 0.1 --flips -1.0 -0.5 0.0 0.5 --kc KC --ld LOCI.ld --nite 50
module load StdEnv/2020 r/4.0.0
Rscript 2_susie_ld.R LOCI-KC-0.0001-500000-z.txt LOCI-KC-0.0001-500000-ztrue.txt LOCI.ld