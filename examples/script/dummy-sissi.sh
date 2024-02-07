#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate rnadeep

python ../train_ali.py --ali-dir ../../rnaconv/data/sissi/dummy/alignments/ \
                    --dbn-dir ../../rnaconv/data/rfam/dummy/seed_neighbourhoods/dbn/ \
                    --model-log-dir ../models \
                    --data-tag sm3-dummy-sissi --smodel 3 --batch-size 50 --epochs 20
