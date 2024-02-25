#!/usr/bin/bash

module load miniconda3 || true
eval "$(conda shell.bash hook)"
source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep

python ../../../../family_generator.py ../../../../sissi099 10 85 ../../../rfam/maxlen290 . > ./generate_o.txt