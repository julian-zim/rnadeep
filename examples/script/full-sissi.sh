#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate rnadeep

python ../train_ali.py --ali-dir ../../rnaconv/data/sissi/full/alignments/ \
                    --dbn-dir ../../rnaconv/data/rfam/full/seed_neighbourhoods/dbn/ \
                    --model-log-dir ../models \
                    --data-tag sm3-full-sissi --smodel 3 --batch-size 50 --epochs 20
