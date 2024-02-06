#!/bin/bash

#SBATCH -p sched_mit_g4nier
#SBATCH -t 100:00:00
#SBATCH -N 2
#SBATCH --mem=0
#SBATCH -J DPANN_drep
#SBATCH -o DPANN_drep
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=izhang@mit.edu

dRep dereplicate DPANN_derep_99 -g all_DPANN_fastas/*.fna --S_algorithm fastANI -sa 0.99 --ignoreGenomeQuality --debug
