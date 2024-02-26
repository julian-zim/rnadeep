#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep || true

python ../../../rfam_converter.py ./Rfam.seed ./seed_trees/original . > ./convert_o.txt
