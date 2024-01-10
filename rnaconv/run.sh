#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate spb

python ./alignment_generator.py 5 ./data/full/rfam/ ./data/full/sissi/ > ./data/full/output.txt
