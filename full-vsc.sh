#./install-vsc.sh

cd ./rnaconv/data/dummy/
./run-vsc.sh

cd ../small
./run-vsc.sh

cd ../ambiguous
./run-vsc.sh

cd ../full
./run-vsc.sh

cd ../../../examples
./batch.sh
