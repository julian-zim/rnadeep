#!/usr/bin/bash

module load miniconda3 || true
eval "$(conda shell.bash hook)"
source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh || true
conda activate rnadeep

python ../train_ali.py --ali-dir ../../rnaconv/data/sissi/full/alignments/ \
                    --dbn-dir ../../rnaconv/data/rfam/full/seed_neighbourhoods/dbn/ \
                    --model-log-dir ../models \
                    --data-tag sm3-full-sissi --smodel 3 --batch-size 50 --epochs 20
