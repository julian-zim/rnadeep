#!/bin/bash

module load miniconda3
eval "$(conda shell.bash hook)"
conda activate rnadeep

DATASET="maxlen250"
TAG="$DATASET-genali"

python ../train_ali.py --ali-dir ../../rnaconv/data/generated/alignment/$DATASET/alignments/ \
                    --dbn-dir ../../rnaconv/data/generated/alignment/$DATASET/neighbourhoods/dbn/ \
                    --model-log-dir ../models \
                    --data-tag $TAG --smodel 3 --batch-size 4 --epochs 5

