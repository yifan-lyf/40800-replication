#!/bin/bash
#SBATCH --job-name=jags_stage1_0.5
#SBATCH --output=stage1_0.5out.txt
#SBATCH --error=stage1_0.5err.txt
#SBATCH --time=20:00:00
#SBATCH --partition=caslake 
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --account=macs30123

module load R/4.2.2
module load jags/4.3.1

Rscript stage1_0.5.R