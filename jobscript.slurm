#!/bin/bash 
#SBATCH --partition=short-serial
#SBATCH --job-name=ENSOTRendGamma
#SBATCH -o %j.out 
#SBATCH -e %j.err
#SBATCH --time=12:00:00

# executable 
conda activate ./env

python src/main.pyatm.mem=EEEEE coup.c_d=0.00225 atm.eps_days=0.75 atm.e_frac=0.5 atm.vary_cloud_const=False
python src/main.pyatm.mem=EEEEE coup.c_d=0.00225 atm.eps_days=0.75 atm.e_frac=1.0 atm.vary_cloud_const=False
python src/main.pyatm.mem=EEEEE coup.c_d=0.00225 atm.eps_days=0.75 atm.e_frac=1.0 atm.vary_cloud_const=True
python src/main.pyatm.mem=EEEEE coup.c_d=0.00225 atm.eps_days=0.75 atm.e_frac=2.0 atm.vary_cloud_const=False
python src/main.pyatm.mem=EEEEE coup.c_d=0.00225 atm.eps_days=0.75 atm.e_frac=2.0 atm.vary_cloud_const=True
