#!/bin/bash
#SBATCH -p sched_mit_g4nier
#SBATCH -t 100:00:00
#SBATCH -N 2
#SBATCH --mem=0
#SBATCH -J DPANN_anvio2
#SBATCH -o DPANN_anvio2
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=izhang@mit.edu

for i in *.db
do
anvi-run-kegg-kofams -c $i --num-threads 2 --kegg-data-dir /home/izhang/KEGG_archive_unpacked/KEGG
done

for i in *fa.db
do
anvi-estimate-metabolism -c $i -O $i --kegg-data-dir /home/izhang/KEGG_archive_unpacked/KEGG
done
