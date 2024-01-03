#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate spb

#python ./alignment_generator.py 2 ./data/rfam/ ./data/sissi/ > ./data/output.txt
python ./alignment_generator.py 1 ./data/dummy/rfam/ ./data/dummy/sissi/ > ./data/dummy/output.txt
