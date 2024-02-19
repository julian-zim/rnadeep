#!/usr/bin/bash

module load miniconda3 || true
eval "$(conda shell.bash hook)"
source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep

cd ../../..
python ./rfam_filter.py data/rfam/full > ./data/rfam/full/filter_o.txt
