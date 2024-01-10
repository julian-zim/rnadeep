#!/usr/bin/bash

source /home/julian-zim/Programs/anaconda3/etc/profile.d/conda.sh
conda activate spb

                    #--ali-dir ../rnaconv/data/ambiguous/sissi/ \
python train_ali.py --ali-dir ../rnaconv/data/ambiguous/rfam/seed_frequency/ali/ \
                    --dbn-dir ../rnaconv/data/ambiguous/rfam/seed_neighbourhood/dbn/ \
                    --model-log-dir models \
                    --data-tag rfam-sissi --smodel 0 --batch-size 10 --epochs 2

python train_ali.py --ali-dir ../rnaconv/data/ambiguous/rfam/seed_frequency/ali/ \
                    --dbn-dir ../rnaconv/data/ambiguous/rfam/seed_neighbourhood/dbn/ \
                    --model-log-dir models \
                    --data-tag rfam-sissi --load-model ./models/sm0_rfam-sissi_002 --epoch0 2 --batch-size 10 --epochs 4
