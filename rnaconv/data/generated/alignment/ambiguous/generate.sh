#!/usr/bin/bash

module load miniconda3 || true
eval "$(conda shell.bash hook)"
source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep

python ../../../../alignment_generator.py ../../../../sissi099 3 ../../../rfam/ambiguous . > ./generate_o.txt
