#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate spb

python train.py --train-data-file ./data/uniform_len25-30_n10000.fa-train \
                --validation-data-file ./data/uniform_len25-30_n2000.fa-test \
                --model-log-dir ./models/seq/ \
                --data-tag l30 --smodel 3 --batch-size 50 --epochs 5
