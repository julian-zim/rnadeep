#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate rnadeep

cd ../../..
python ./alignment_generator.py 3 data/rfam/dummy data/sissi/dummy > ./data/sissi/dummy/output.txt
