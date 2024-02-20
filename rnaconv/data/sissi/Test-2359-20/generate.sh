#!/usr/bin/bash

module load miniconda3 || true
eval "$(conda shell.bash hook)"
source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep

cd ../../..
python ./alignment_generator.py 20 data/rfam/Test-2359-20 data/sissi/Test-2359-20 > ./data/sissi/Test-2359-20/generate_o.txt
