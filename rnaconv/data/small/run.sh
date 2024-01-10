#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate spb

python ../../alignment_generator.py 5 ./rfam/ ./sissi/ > ./output.txt
