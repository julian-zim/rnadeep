#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate rnadeep

cd ../../..
python ./rfam_converter.py data/rfam/full > ./data/rfam/full/output.txt
