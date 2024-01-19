#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate spb

python train_seq_as_ali.py --ali-dir ./data/seq_as_ali/ali/ \
                           --dbn-dir ./data/seq_as_ali/dbn/ \
                           --model-log-dir models/milestone \
                           --data-tag l30ali --smodel 3 --batch-size 50 --epochs 20 \
                           2>> models/milestone/l30ali.err >> models/milestone/l30ali.out
#--ali-dir ./data/seq_as_ali/ali/ --dbn-dir ./data/seq_as_ali/dbn/ --model-log-dir models/milestone --data-tag l30ali --smodel 3 --batch-size 50 --epochs 20
