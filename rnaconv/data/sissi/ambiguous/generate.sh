#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep || true

cd ../../..
python ./alignment_generator.py 1 data/rfam/ambiguous data/sissi/ambiguous > ./data/sissi/ambiguous/output.txt
