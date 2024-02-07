#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate rnadeep

cd ../../..
python ./alignment_generator.py 1 data/rfam/full data/sissi/full > ./data/sissi/full/output.txt
