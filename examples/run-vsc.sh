#!/bin/sh
#SBATCH -J rnadeep
#SBATCH --partition=zen3_0512_a100x2
#SBATCH --qos zen3_0512_a100x2
#SBATCH --gres=gpu:2
#SBATCH --mail-user=a12144285@unet.univie.ac.at
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/home/fs71475/julianz123/rnadeep/examples/slurm_out/gcn_slurm.%A.out
#SBATCH --error=/home/fs71475/julianz123/rnadeep/examples/gcn_slurm.%A.error

spack load cuda
spack load miniconda3@4.12.0
source /home/fs71475/julianz123/.bashrc
conda activate /home/fs71475/julianz123/.conda/envs/jzspb

python train_ali.py --ali-dir ../rnaconv/data/small/sissi/ \
                    --dbn-dir ../rnaconv/data/small/rfam/seed_neighbourhoods/dbn/ \
                    --model-log-dir models \
                    --data-tag rfam-sissi --smodel 0 --batch-size 10 --epochs 2

python train_ali.py --ali-dir ../rnaconv/data/small/rfam/seed_alignments/ \
                    --dbn-dir ../rnaconv/data/small/rfam/seed_neighbourhoods/dbn/ \
                    --model-log-dir models \
                    --data-tag rfam-sissi --load-model ./models/sm0_rfam-sissi_002 --epoch0 2 --batch-size 10 --epochs 4
