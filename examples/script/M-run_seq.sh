#!/usr/bin/bash

module load miniconda3 || true
eval "$(conda shell.bash hook)"
source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep

python ../train.py --train-data-file ../data/uniform_len25-30_n10000.fa-train \
                   --validation-data-file ../data/uniform_len25-30_n2000.fa-test \
                   --model-log-dir ../models \
                   --data-tag sm3-l25-30 --smodel 3 --batch-size 5 --epochs 5
