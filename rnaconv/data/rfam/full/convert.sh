#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep || true

cd ../../..
python ./rfam_converter.py data/rfam/full > ./data/rfam/full/convert_o.txt
