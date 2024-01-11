#!/bin/bash

scp README.md install-vsc.sh setup.py julianz123@vsc5.vsc.ac.at:/home/fs71475/julianz123/rnadeep/
scp -r rnadeep julianz123@vsc5.vsc.ac.at:/home/fs71475/julianz123/rnadeep/
scp examples/batch.sh examples/run-vsc.slrm examples/train_ali.py examples/predict_ali.py julianz123@vsc5.vsc.ac.at:/home/fs71475/julianz123/rnadeep/examples/
scp rnaconv/run-vsc.sh rnaconv/sissi099 rnaconv/alignment_generator.py rnaconv/rfam_converter.py julianz123@vsc5.vsc.ac.at:/home/fs71475/julianz123/rnadeep/rnaconv/
scp -r rnaconv/data julianz123@vsc5.vsc.ac.at:/home/fs71475/julianz123/rnadeep/rnaconv/
