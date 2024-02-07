#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate rnadeep

cd ../../..
python ./alignment_generator.py 3 data/rfam/small data/sissi/small > ./data/sissi/small/output.txt
