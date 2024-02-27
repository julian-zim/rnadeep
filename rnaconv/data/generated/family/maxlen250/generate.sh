#!/usr/bin/bash

module load miniconda3 || true
eval "$(conda shell.bash hook)"
source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep

python ../../../../data_generator.py \
	   --sissi-path ../../../../sissi099 --type families --number 5 --min-length 20 --max-length 250 \
	   --tree-path ../../../rfam/maxlen250/seed_trees/rescaled \
	   --single-freq-path ../../../rfam/maxlen250/seed_frequencies/single \
	   --doublet-freq-path ../../../rfam/maxlen250/seed_frequencies/doublet \
	   --out-path . \
	   #> ./generate_o.txt
