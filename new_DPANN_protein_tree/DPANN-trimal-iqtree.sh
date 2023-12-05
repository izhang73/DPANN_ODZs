#!/bin/bash

#SBATCH -p sched_mit_g4nier
#SBATCH -t 100:00:00
#SBATCH -N 2
#SBATCH --mem=0
#SBATCH -J DPANN_tree_2
#SBATCH -o DPANN_tree_2
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=izhang@mit.edu

module add engaging/iqtree/1.6.3
wait
iqtree -nt AUTO -s all_pmo_Cox_nosZ_new_MAFFT_trimal.faa -alrt 1000 -bb 1000 -keep-ident
