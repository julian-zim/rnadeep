## INSTALLATION

* Clone the repository
* Run `pip install .` in the root directory. It's recommended to create a virtual environment first, use `rnadeep.yaml` for that (e.g. `conda env create -n rnadeep --file rnadeep.yaml`).


## DIRECTORIES
* *rnadeep*: The nessecary tools to sample and encode the learning data, the spotrna models, and important metrics.
* *rnaconv*: Contains the rfam database as well as multiple tools to generate artifical data in different ways using RNAfold<sup>3</sup> and SISSI<sup>4</sup>. Refer to `rnaconv/ReadMe.md`
* *examples*: Example python, bash and slurm scripts to train and predict data. Refer to `examples/ReadMe.md`

* <sub><sup>3</sup>https://github.com/ViennaRNA</sub><br/>
* <sub><sup>4</sup>https://cibiv.github.io/software/sissi</sub><br/>


## DOCUMENTATION
* For a documentation document of the modules and functions refer to `doc/documentation.pdf`


# Caveats to deep learning approaches to RNA secondary structure prediction

Christoph Flamm<sup>1</sup>, Julia Wielach<sup>1</sup>, Michael T. Wolfinger<sup>1,2</sup>, Stefan Badelt<sup>1</sup>,  Ronny Lorenz<sup>1</sup>, Ivo L. Hofacker<sup>1,2</sup>

<sub><sup>1</sup>Department of Theoretical Chemistry, University of Vienna, Vienna, Austria</sub><br/>
<sub><sup>2</sup>Research Group Bioinformatics and Computational Biology, Faculty of Computer Science, University of Vienna, Vienna, Austria</sub><br/>

This repository contains additional resources that were used during preparation of the manuscript. These include custom code, Google Colab notebook and data files.
