## INSTALLATION
- Create a virtual environment: `conda env create -n rnadeep --file rnadeep.yaml`
- Run `pip install -r requirements` in root.

## DIRECTORIES
* *rnadeep*: The nessecary tools to sample and encode the learning data, the spotrna models, and important metrics.
* *rnaconv*: Contains the rfam database as well as multiple tools to generate artifical data in different ways using RNAfold<sup>3</sup> and SISSI<sup>4</sup>. Refer to `rnaconv/ReadMe.md`
* *examples*: Example python, bash and slurm scripts to train and predict data. Refer to `examples/ReadMe.md`

## REFS
* <sub><sup>3</sup>https://github.com/ViennaRNA</sub><br/>
* <sub><sup>4</sup>https://cibiv.github.io/software/sissi</sub><br/>

## DOCUMENTATION

[![Preview](preview.png)](doc/documentation.pdf)
