#!/bin/bash

module load miniconda3
eval "$(conda shell.bash hook)"
conda activate rnadeep

DATASET="dummy"
TAG="$DATASET-rfam"

python ../train_ali.py --ali-dir ../../rnaconv/data/rfam/$DATASET/seed_alignments/ \
                       --dbn-dir ../../rnaconv/data/rfam/$DATASET/seed_neighbourhoods/dbn/ \
                       --model-log-dir ../models \
                       --data-tag $TAG --smodel 3 --batch-size 4 --epochs 5
