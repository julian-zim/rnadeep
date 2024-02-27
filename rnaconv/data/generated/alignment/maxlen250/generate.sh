#!/usr/bin/bash

module load miniconda3 || true
eval "$(conda shell.bash hook)"
source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep

python ../../../../data_generator.py \
	   --sissi-path ../../../../sissi099 --type alignments --number 5 \
	   --tree-path ../../../rfam/maxlen250/seed_trees/rescaled \
	   --neigh-path ../../../rfam/maxlen250/seed_neighbourhoods/dbn \
	   --single-freq-path ../../../rfam/maxlen250/seed_frequencies/single \
	   --doublet-freq-path ../../../rfam/maxlen250/seed_frequencies/doublet \
	   --ali-path ../../../rfam/maxlen250/seed_alignments \
	   --out-path . \
	   #> ./generate_o.txt
