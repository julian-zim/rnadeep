#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate rnadeep

python train_ali.py --ali-dir ../rnaconv/data/sissi/small/alignments/ \
                    --dbn-dir ../rnaconv/data/rfam/small/seed_neighbourhoods/dbn/ \
                    --model-log-dir models \
                    --data-tag sm3-small-sissi --smodel 3 --batch-size 10 --epochs 20

#python train_ali.py --ali-dir ../rnaconv/data/rfam/small/seed_alignments/ \
                    #--dbn-dir ../rnaconv/data/rfam/small/seed_neighbourhoods/dbn/ \
                    #--model-log-dir models \
                    #--data-tag rfam-sissi --load-model ./models/sm0_rfam-sissi_002 --epoch0 10 --batch-size 50 --epochs 20
