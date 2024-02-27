## MODULES
* *rfam_converter.py*:<br>
&nbsp;&nbsp;&nbsp;&nbsp;Tool to convert the Rfam database into filetypes and formats that can be used by SISSI.<br>
&nbsp;&nbsp;&nbsp;&nbsp;Usage: ./rfam_converter.py <seed_file_path> <tree_directory_path> <out_path>

* *rfam_filter.py*:<br>
&nbsp;&nbsp;&nbsp;&nbsp;Tool to filter the converted rfam database by throwing out all alignments (and corresponding consensus structure, frequencies and tree) exceeding a certain length<br>
&nbsp;&nbsp;&nbsp;&nbsp;Usage: ./rfam_filter.py <rfam_directory_path> \[max_length\]

* *data_generator.py*:<br>
&nbsp;&nbsp;&nbsp;&nbsp;Tool to generate n RNA families of a certain length or alignments based on a certain family using SISSI, for given single and doublet equilibrium frequencies and phylogenetic trees.<br>
&nbsp;&nbsp;&nbsp;&nbsp;Arguments:<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--sissi-path:        Path to compiled sissi099 file<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--type:              Type of data to generate: can be "alignments" "families"<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--number:            Number of data points to generate<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--min-length:        Minimum length of sequences and alignments to generate. Only used of type=families<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--max-length:        Maximum length of sequences and alignments to generate. Only used of type=families<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--tree-path:         Path to the tree files used for generation<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--neigh-path:        Path to the neighbourhood files used for generation. Only used if type=alignments<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--single-freq-path:  Path to single equilibrium frequency files used for generation<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--doublet-freq-path: Path to doublet equilibrium frequency files used for generation<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--ali-path:          Path to the alignment files used to readd indels. Only used if type=alignments. Can be omitted if nothing should be readded<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--out-path:          Path to the directory in which to save the generated files<br>

* *data_filter.py*:<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Tool to filter the alignments created by the data_generator by removing sequences that deviate too much from the desired consensus structure of their alignment. Also contains code to extract the equilibrium frequencies from the generated alignments and form the differences between them and the original equilibrium frequencies of the rfam alignments that were used for generation.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Usage: ./data_filter.py <data_directory_path> \[<max_ss_deviation>\]<br>

* *data/copy.py*:<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Tool to copy selected data points from a given rfam directory into a new location for creating sub data sets.<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Usage: ./copy.py <path_to_copy_to> <rfam_path_to_copy_from> <filenames_to_copy> \[<additional_filenames_to_copy> ...\]<br>


## DATASETS
There are three types of datasets:
* real data (data/rfam/)
* generated alignments (data/generated/alignment/)
* generated families (data/generated/family/)
<br>
Each path contains a dummy folder with example scripts of how to create the types of datasets. Refer to them for more information.


## Typical Workflow
### Prepare the real data:
* run data/rfam/full/convert.sh to convert the Rfam database into a state usable by SISSI
* run data/rfam/full/filter.sh to do basic filtering to the converted database (removing all alignments longer than 700nt)

* Create subdatasets: Create a folder for your dataset and call the rfam/copy.py module to copy specified files from the converted Rfam database, similar to the example in data/rfam/dummy/copy.sh.
* Filter subdatasets: Call the rfam_filter.py module to do basic filtering to your subdataset, similar to the example in data/rfam/dummy/filter.sh. If you already filtered the converted Rfam database, this doesn't do anything.

### Generate artificial data:
* Call the data_generator.py module to either generate more alignments based on certain consensus structures and evolutionary trees, or to generate more families based on certain evolutionary trees, similar to the examples in data/generated/alignment/dummy/generate.sh and data/generated/family/dummy/generate.sh.
* Call the data_filter.py module to filter the generated alignments by throwing out sequences which predicted secondary structure deviate too much from the desired consensus structure of its alignment, similar to the examples in data/generated/alignment/dummy/filter.sh and data/generated/family/dummy/filter.sh.

### Train on real or artificial data:
* Use the ../rnadeep/sample_ali.py module to sample data from artificial or real datasets and ../examples/train_ali.py module to train a model on it, similar to the example scripts starting with "dummy-" in ../examples/script/ as well as their corresponding slurm scripts in ../examples/slurm/.
