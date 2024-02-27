#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate rnadeep

mkdir html
cd ./html
python ../create.py
