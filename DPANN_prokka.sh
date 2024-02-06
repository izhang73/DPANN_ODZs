#!/bin/bash

#SBATCH -p sched_mit_g4nier
#SBATCH -t 100:00:00
#SBATCH -N 2
#SBATCH --mem=0
#SBATCH -J DPANN_prokka
#SBATCH -o DPANN_prokka
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=izhang@mit.edu

for genome in $(cat all-DPANN.txt)
do
prokka --kingdom Archaea "$genome" --outdir $genome.prokka.output --metagenome --cpus 8
done
