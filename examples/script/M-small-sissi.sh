#!/usr/bin/bash

module load miniconda3 || true
eval "$(conda shell.bash hook)"
source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep

python ../train_ali.py --ali-dir ../../rnaconv/data/generated/alignment/small/alignments/ \
                    --dbn-dir ../../rnaconv/data/rfam/small/seed_neighbourhoods/dbn/ \
                    --model-log-dir ../models \
                    --data-tag sm3-small-sissi --smodel 3 --batch-size 5 --epochs 5
