#!/bin/bash

#SBATCH -p sched_mit_g4nier
#SBATCH -t 100:00:00
#SBATCH -N 2
#SBATCH --mem=0
#SBATCH -J DPANN_gtotree
#SBATCH -o DPANN_gtotree
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=izhang@mit.edu

GToTree -f all_DPANN.txt -H Archaea -o DPANN_derep_GToTree_out -B -G 0.25 -t
