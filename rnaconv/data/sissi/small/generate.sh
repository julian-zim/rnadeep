#!/usr/bin/bash

module load miniconda3 || true
eval "$(conda shell.bash hook)"
source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep

cd ../../..
python ./alignment_generator.py 3 data/rfam/small data/sissi/small > ./data/sissi/small/generate_o.txt
