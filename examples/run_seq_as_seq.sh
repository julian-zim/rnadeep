#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate spb

python train_seq_as_seq.py --train-data-file ./data/uniform_len25-30_n12000.fa \
                           --validation-data-file ./data/uniform_len25-30_n12000.fa \
                           --model-log-dir models/milestone \
                           --data-tag l30 --smodel 3 --batch-size 50 --epochs 20 \
                           2>> models/milestone/l30.err >> models/milestone/l30.out
#--train-data-file ./data/uniform_len25-30_n12000.fa --validation-data-file ./data/uniform_len25-30_n12000.fa --model-log-dir models/milestone --data-tag l30 --smodel 3 --batch-size 50 --epochs 20
