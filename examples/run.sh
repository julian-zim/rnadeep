#python train.py -t ./data/uniform_len25-30_n10000.fa-train -v ./data/uniform_len25-30_n2000.fa-test --data-tag l30 --smodel 3 --batch-size 50 --epochs 5
python train_ali.py --ali-dir ../rnaconv/data/sissi/ --dbrs-dir ../rnaconv/data/rfam/seed_neighbourhood/dbs/ --data-tag l30 --smodel 3 --batch-size 10 --epochs 5
