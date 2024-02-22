#!/usr/bin/bash

module load miniconda3 || true
eval "$(conda shell.bash hook)"
source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep

cd ../../..
python ./alignment_generator.py 3 data/rfam/ambiguous data/sissi/ambiguous > ./data/sissi/ambiguous/generate_o.txt
