#!/usr/bin/bash

module load miniconda3 || true
eval "$(conda shell.bash hook)"
source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep

cd ../../..
python ./sissi_filter.py data/rfam/Test-1483-10 data/sissi/Test-1483-10 > ./data/sissi/Test-1483-10/filter_o.txt
