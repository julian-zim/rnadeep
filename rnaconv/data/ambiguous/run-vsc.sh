#!/bin/sh

spack load miniconda3@4.12.0
source /home/fs71475/julianz123/.bashrc
conda activate /home/fs71475/julianz123/.conda/envs/jzspb

python ../../alignment_generator.py 1 ./rfam/ ./sissi/ > ./output.txt
python filter.py
