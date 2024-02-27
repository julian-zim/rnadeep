## MODULES
*rfam_converter.py*:<br>
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

## Typical Workflow
### Prepare the real data:
* run data/rfam/full/convert.sh to convert the Rfam database into a state usable by SISSI
* run data/rfam/full/filter.sh to do basic filtering to the converted database

* Create subdatasets: Create a folder in which you run a script to copy specified files from the converted database, similar to the example in data/rfam/dummy/copy.sh. Alternatively, you can use a script to get filesnames to copy based on certain properties, similar to the examples in data/rfam/maxlen250/copy.sh.
* Filter subdatasets: Run a script to do basic filtering the created subdataset, similar to the example in data/rfam/dummy/filter.sh. If the converted Rfam database was already filtered, this doesn't do anything.

### Generate artificial data:
* Create a folder in which you run a script to generate either more alignments baased on certain families or more families based on certain evolutionary trees, similar to the examples in data/rfam/generated/alignment/dummy/generate.sh and data/rfam/generated/families/dummy/generate.sh.
* Filter your artificial data: Run a script to filter the generated data by throwing out data that deviates too much from the desired data, similar to the examples in data/rfam/generated/alignment/dummy/filter.sh and data/rfam/generated/families/dummy/filter.sh.

### Train on real or artificial data:
* Use the alignment sampler to sample data from the artificial or real datasets, similar to the examples in examples/script/dummy-rfam.sh, examples/script/dummy-genali.sh and examples/script/dummy-genfam.sh as well as their corresponding slurm scripts in examples/slurm.
